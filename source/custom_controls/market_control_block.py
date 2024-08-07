
from re import I
import sys

from market_system.Markets import Market
#from market_system.Assets import NondispatchableAsset_3ph
from market_system.EnergySystem import EnergySystem
from market_system.Network3phOpenDSS import Network_3ph
from market_system.DSS_info_extraction import dssConvert

import pickle
import json
import helics as h
import logging
logger = logging.getLogger(__name__)

from global_aux import OpenDSS_message_types

#####
# this controller requires a horizon load and market variables as inputs
# it optimizes EV and BTMS dispatching and outputs net site load
#####

class LPMarketController():
    def __init__(self, name='greedy', helics_config_path='', timestep_sec=60*5, feeder_name='ieee_34', dss_file_name='Main.dss', horizon_sec=24*60*60, evse_df=[], stationary_storage_df=[]):
        # add important params here
        self.name = name
        self.dss_file_name = dss_file_name
        self.feeder_name = feeder_name
        print(f'feeder_name in market_control_block: {feeder_name}')
        print(f'dss_file_name in market_control_block: {dss_file_name}')
        src_df, bus_df, load_df, line_df, solution_df, Ybus = dssConvert(feeder_name,dss_file_name)
        self.src_df = src_df # voltage sources on opendss model
        self.bus_df = bus_df # bus df of voltage base, name, number, connection type, P, Q, phases
        self.load_df = load_df # load df of name, number, bus name, connection type, P, Q, phases
        self.horizon_sec = horizon_sec
        self.timestep_sec = timestep_sec
        #self.line_df = line_df
        #self.solution_df = solution_df
        #self.Ybus = Ybus
        self.evse_df = evse_df # this should come from a csv and define the Pmin and max, as well as Emin and max, E0, Ehorizon, and busid
        self.stationary_storage_df = stationary_storage_df # this should come from a csv and define the Pmin and max, Emin and max, E0, Ehorizon, and busid
        # we can also assume that E0=Ehorizon if we want conservation of energy
        # this controller also requires market data
        self.bus_id_market = 0
        self.prices_export = []
        self.demand_charge = [] 
        self.nondispatch_assets = []
        #  this controller also requires a network



    def setup_market_controller(self, prices_export, demand_charge, Pmax_market=1000000, Pmin_market=-100000):       
        hs = self.horizon_sec
        ts = self.timestep_sec
        if not isinstance(Pmax_market, list):
            Pmax_market = [Pmax_market]*int(hs/ts)
        if not isinstance(Pmin_market, list):
            Pmin_market = [Pmin_market]*int(hs/ts)
        # setup market and energy management objects

        if self.name == "TOU_withoutV":            
            prices_import = [[0.2]]* int(24*3600/self.timestep_sec)

        elif self.name == "DA_withoutV":
            prices_import = [[0.15]]* int(24*3600/self.timestep_sec)

        else: #self.name == 'greedy':
            prices_import = [[0.1]]* int(24*3600/self.timestep_sec) #[]

        # pull lists of non-dispatchable and dispatchable assets from opendss model
        storage_assets = self.stationary_storage_df
        self.storage_assets = storage_assets
        evse_assets = self.evse_df
        # in opendss all nondispatchable assets are loadshapes, so just reformat the load_df a little
        nondispatch_assets = self.load_df
        nondispatch_assets['Pnet_pred'] = nondispatch_assets['P']
        nondispatch_assets['Qnet_pred'] = nondispatch_assets['Q']
        nondispatch_assets = nondispatch_assets.rename(columns={'phase':'phases', 'P':'Pnet', 'Q':'Qnet', 'bus_name':'bus_id'})
        #for load in self.load_df:
        #    # make a list of non-dispatchable asstes and their bus numbers
        #    # parse the opendss model
        #    ND_load_ph = NondispatchableAsset_3ph(Pnet=0, Qnet=0, bus_num=load.bus_id, ph_i=load.phases, dt=ts,
        #                                   T=hs, Pnet_pred = 0,
        #                                   Qnet_pred = 0)
        #    nondispatch_assets.append(ND_load_ph)
        self.nondispatch_assets = nondispatch_assets
        nda_list = nondispatch_assets.to_dict('records')
        for nda in nda_list:
            nda['number'] = self.bus_df[self.bus_df['name']==nda['bus_id']]['number'][0]
        bus_id_market = 0 # all buses have the same market, so this identifier can be 0

        self.network = Network_3ph(self.feeder_name)
        market = Market(bus_id_market, prices_export, prices_import, demand_charge, Pmax_market, Pmin_market, ts, hs)
        self.market = market
        energy_system = EnergySystem(storage_assets, nda_list, self.network, market, ts, hs, ts, hs)
        self.energy_system = energy_system


    def solve(self, DSS_state_info_dict):
        # this function solves the control parameters
        ev_control_setpoints = {}
        hs = self.horizon_sec
        ts = self.timestep_sec

        # first get the updated grid status
        #voltages = json.loads(h.helicsInputGetString(self.subscriptions[0]))
        #currents = json.loads(h.helicsInputGetString(self.subscriptions[1]))
        
        if OpenDSS_message_types.get_basenetloads in DSS_state_info_dict.keys():
            non_disp_loads = DSS_state_info_dict[OpenDSS_message_types.get_basenetloads]
        #non_disp_loads = json.loads(h.helicsInputGetString(self.subscriptions[2]))
        # update the non dispatchable asset loads
        self.nondispatch_assets['Pnet'] = self.nondispatch_assets['Pnet_pred']
        self.nondispatch_assets['Qnet'] = self.nondispatch_assets['Qnet_pred']
        for busname_i in non_disp_loads.keys():
            self.nondispatch_assets[self.nondispatch_assets['bus_id']==busname_i]['Pnet_pred'] = non_disp_loads[busname_i][0]
            self.nondispatch_assets[self.nondispatch_assets['bus_id']==busname_i]['Qnet_pred'] = non_disp_loads[busname_i][1]
        #for i_nda in range(len(self.nondispatch_assets)):
        #    busname_i = self.nondispatch_assets.loc[i_nda,'bus_id']
        #    self.nondispatch_assets.loc[i_nda, 'Pnet_pred'] = self.nondispatch_assets[i_nda, 'Pnet']
        #    self.nondispatch_assets[i_nda].Qnet_pred = self.nondispatch_assets[i_nda, 'Qnet']
        #    self.nondispatch_assets[i_nda].Pnet = non_disp_loads[busname_i][0]
        #    self.nondispatch_assets[i_nda].Qnet = non_disp_loads[busname_i][1]
        # update energy system with new asset loads
        nda_list = self.nondispatch_assets
        nda_list = nda_list.to_dict('records')
        #self.energy_system = EnergySystem(self.storage_assets, nda_list, self.network, self.market, ts, hs, ts, hs)

        network = self.network

        i_line_unconst_list = list(range(network.N_lines))   
        v_bus_unconst_list = list(range((network.N_phases)*(network.N_buses-1))) # no voltage constraints 
        EMS_output = self.energy_system.simulate_network_3phPF('3ph',\
                                        i_unconstrained_lines=\
                                        i_line_unconst_list,\
                                        v_unconstrained_buses=\
                                        v_bus_unconst_list)                                    
                                  
        #pickle.dump(EMS_output, open(join(EMS_path_string, normpath('Month' + str(Case_Month) + '_Day' + str(day) + '_EMS_output_' + str(x) + '.p')), "wb"))    

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EMS Result post-processing and Analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Step 1:  Extract the original EMS results  
        #PF_network_res = EMS_output['PF_network_res']
        P_import_ems = EMS_output['P_import_ems']
        P_export_ems = EMS_output['P_export_ems']
        #P_ES_ems = EMS_output['P_ES_ems']
        #P_demand_ems = EMS_output['P_demand_ems']  
        for load_name in self.load_df.keys():
            ev_control_setpoints[load_name] = P_import_ems[load_name] - P_export_ems[load_name]

        self.control_setpoints = ev_control_setpoints
        return ev_control_setpoints
