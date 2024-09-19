
from re import I
import sys
import os

from market_system.Markets import Market
#from market_system.Assets import NondispatchableAsset_3ph
from market_system.EnergySystem import EnergySystem
from market_system.Network3phOpenDSS import Network_3ph
from market_system.DSS_info_extraction import dssConvert

import opendssdirect as dss

import pickle
import json
import helics as h
import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from global_aux import OpenDSS_message_types

#####
# this controller requires a horizon load and market variables as inputs
# it optimizes EV and BTMS dispatching and outputs net site load
#####

class LPMarketController():
    def __init__(self, simulation_time_constraints, name='greedy', helics_config_path='', feeder_name='ieee_34', dss_file_name='Main.dss', evse_df=[], stationary_storage_df=[]):
        # add important params here
        self.name = name
        self.dss_file_name = dss_file_name
        self.feeder_name = feeder_name
        print(f'feeder_name in market_control_block: {feeder_name}')
        print(f'dss_file_name in market_control_block: {dss_file_name}')
        #src_df, bus_df, load_df, line_df, solution_df, Ybus = dssConvert(feeder_name,dss_file_name)
        #self.src_df = src_df # voltage sources on opendss model
        dss.Command(f'Redirect {dss_file_name}')
        LoadColumns = ["name","number","bus_name","connect","P","Q","phase_number","phase"]
        Load_names = dss.Loads.AllNames()
        #load_df = pd.DataFrame(data=np.zeros((nLoads,len(LoadColumns))), index=Load_names, columns=LoadColumns)
        load_df = pd.DataFrame(index=Load_names, columns=LoadColumns)
        bus_df = {}
        bus_names= dss.Circuit.AllBusNames()
        bus_df['name'] = bus_names
        bus_df['number'] = np.arange((len(bus_names)))

        load_df['name'] = Load_names
        load_df['number'] = np.arange((len(Load_names)))
        self.bus_df = bus_df # bus df of voltage base, name, number, connection type, P, Q, phases
        self.load_df = load_df # load df of name, number, bus name, connection type, P, Q, phases
        self.horizon_sec = simulation_time_constraints.end_simulation_unix_time
        self.timestep_sec = simulation_time_constraints.grid_timestep_sec
        self.start_time_sec = simulation_time_constraints.start_simulation_unix_time
        #self.line_df = line_df
        #self.solution_df = solution_df
        #self.Ybus = Ybus
        self.evse_df = evse_df # this should come from a csv and define the Pmin and max, as well as Emin and max, E0, Ehorizon, and busid
        self.stationary_storage_df = stationary_storage_df # this should come from a csv and define the Pmin and max, Emin and max, E0, Ehorizon, and busid
        # we can also assume that E0=Ehorizon if we want conservation of energy
        # only run the optimization once and then query the outputs
        # this controller also requires market data
        self.bus_id_market = 0
        self.prices_export = []
        self.demand_charge = [] 
        self.nondispatch_assets = []
        self.evse_assets = [] # filled in the setup function
        #  this controller also requires a network




    def setup_market_controller(self, prices_export, demand_charge, Pmax_market=1000000, Pmin_market=-100000, ce_file='inputs/CE_Sep_Shellbank_22700.csv', pev_file='inputs/EV_inputs.csv', evse_types_file='inputs/EVSE_inputs.csv'):   #TODO    
        hs = (self.horizon_sec-self.start_time_sec)
        ts = self.timestep_sec
        n_timesteps = int(hs/ts)
        if not isinstance(Pmax_market, list):
            Pmax_market = [Pmax_market]*int(hs/ts)
        if not isinstance(Pmin_market, list):
            Pmin_market = [Pmin_market]*int(hs/ts)
        # setup market and energy management objects
        market_based_control_names = ['TOU_withoutV', 'DA_withoutV', 'greedy', 'LMP_dayahead', 'LMP_realtime']
        emissions_based_control_names = ['PJM_emissions', 'emissions']

        # start 2nd of Sept 2022, end 4 of Sept 2022
        if self.name in market_based_control_names:
            hs_market = hs # assume 3 days of simulation
            ts_market = 3600 # assume hourly timesteps
            LMP_path_string = os.path.join(os.getcwd(), 'inputs','market_prices')
            day_ahead_path_string = os.path.join(LMP_path_string,"Dayahead_LMP_hourly_2022.csv")
            print(f'day_ahead_path_string: {day_ahead_path_string}')
            open(day_ahead_path_string, 'r')
            DA_LMP_hourly_2022 = pd.read_csv(day_ahead_path_string)
            start_lmp_index = DA_LMP_hourly_2022.index[DA_LMP_hourly_2022['Time'] == '2022-09-02 00:00:00-04:00'].tolist()[0]
            end_lmp_index = start_lmp_index + int(hs_market/ts_market) - 1
            if self.name == "TOU_withoutV":            
                #prices_import = [[0.2]]* n_timesteps
                prices_import = DA_LMP_hourly_2022.iloc[start_lmp_index:end_lmp_index]['TOU_cent'].values
            elif self.name == "DA_withoutV":
                #prices_import = [[0.15]] * n_timesteps
                prices_import = DA_LMP_hourly_2022[start_lmp_index:end_lmp_index]['LMP_cent_adjust'].values
            elif self.name == "LMP_dayahead":
                # adjusted values are scaled such that the time averaged LMP is equalt to the time averaged TOU
                prices_import = DA_LMP_hourly_2022[start_lmp_index:end_lmp_index]['LMP_cent_adjust'].values
            elif self.name == "LMP_realtime":
                RT_LMP_hourly_2022 = pd.read_csv(os.path.join(LMP_path_string,"RT_LMP_hourly_2022.csv"))
                start_lmp_index = RT_LMP_hourly_2022[RT_LMP_hourly_2022['Time'] == '2022-09-01 00:00:00-04:00'].tolist()[0]
                end_lmp_index = start_lmp_index + int(hs_market/ts_market) - 1
                # adjusted values are scaled such that the time averaged LMP is equal to the time averaged TOU
                prices_import = RT_LMP_hourly_2022[start_lmp_index:end_lmp_index]['LMP_cent_adjust'].values
            else:
                prices_import = DA_LMP_hourly_2022.iloc[start_lmp_index:end_lmp_index]['TOU_cent'].values #[[0.15]] * n_timesteps
        elif self.name in emissions_based_control_names:
            # for the emissions based control, read from a csv
            # and then scale to be similar to typical (divide by 100)
            demand_charge = 0
            carbon_file_path = os.path.join(os.getcwd(), 'inputs', 'market_prices', "Cambium23_MidCase_hourly_PJM_East_2030_srmer_co2e.csv")
            carbon_df = pd.read_csv(carbon_file_path)

            # Convert timestamp to datetime
            carbon_df['timestamp'] = pd.to_datetime(carbon_df['timestamp'])

            # Extract month, day, and hour from the timestamp
            carbon_df['month'] = carbon_df['timestamp'].dt.month
            carbon_df['day'] = carbon_df['timestamp'].dt.day
            carbon_df['hour'] = carbon_df['timestamp'].dt.hour
            # get data for Sept 2, 2030
            carbon_df = carbon_df[carbon_df['month']==9]
            carbon_df = carbon_df[carbon_df['day']==2]
            prices_import = carbon_df['srmer_co2e'].values
            # get the timestep and horizon for the carbon values
            ts_market = (carbon_df['hour'].values[1] - carbon_df['hour'].values[0])*3600
            hs_market = ts*len(prices_import)
        else: #self.name == 'greedy':
            prices_import = [[0.1]]* n_timesteps #[]
            hs_market = hs
            ts_market = ts
        print(f'import pricing based on {self.name} style control: {prices_import}')

        # pull lists of non-dispatchable and dispatchable assets from opendss model

        # create an evse storage asset for each charge event with
        # timeseries and Emin and Emax which dictate if the vehicle is plugged
        # and the end state of charge
        evse_types = pd.read_csv(evse_types_file)
        evse_assets = self.evse_df#self.stationary_storage_df
        n_evse = len(self.evse_df['SE_id'])
        evse_assets['Emax'] = [0 for i in range(n_evse)] # assume no evs at station until updated with charge events
        evse_assets['Emin'] = [0 for i in range(n_evse)] 
        evse_assets['Pmax'] = [0 for i in range(n_evse)] 
        evse_assets['Pmin'] = [0 for i in range(n_evse)] 
        evse_assets['Pmax_abs'] = max(evse_assets['Pmax'])*n_evse
        evse_assets['E0'] = [np.zeros((n_timesteps,1)) for i in range(n_evse)] #.1*evse_assets['Emax']
        #TODO: add integration of charge event limits
        #final_energy[-1] = 70 #assume wanting full charge for 70kWh battery
        evse_assets['ET'] = [np.zeros((n_timesteps,1)) for i in range(n_evse)] # this becomes a list
        evse_assets['dt_ems'] = [ts/3600]*n_evse
        evse_assets['T_ems'] = [24*12]*n_evse
        evse_assets['Pnet'] = [np.zeros((n_timesteps,1)) for i in range(n_evse)]
        evse_assets['Qnet'] = [np.zeros((n_timesteps,1)) for i in range(n_evse)]
        evse_assets['c_deg_lin'] = [0 for i in range(n_evse)] 
        evse_assets['eff'] = [np.ones(100)]*n_evse
        evse_assets['eff_opt'] = [1]*n_evse
        evse_assets['bus_id'] = evse_assets['node_id']

        #evse_assets = pd.DataFrame(evse_assets)
        # read in the EV specs to get the battery sizes
        ev_df = pd.read_csv(pev_file)
        # first read in the CE (charge event) file and determine the start and stop times 
        # as well as the final energy needed
        ce_df = pd.read_csv(ce_file)
        # iterate through SE (supply equipment) and get the list of CE for it
        se_iter = 0
        for se_id in evse_assets['SE_id']:
            ce_at_se = ce_df[ce_df['SE_id']==se_id]
            se_type = evse_assets[evse_assets['SE_id']==se_id]['SE_type'].item()
            max_evse = evse_types[evse_types['EVSE_type']==se_type]['AC/DC_power_limit_kW'].item()
            evse_assets['Pmax'][se_iter] = max_evse
            # assume there is more than one
            pev_battery_sizes = [0]
            print(f'ces at se {len(ce_at_se)}')
            for _, ce_i in ce_at_se.iterrows():
                ev_battery_size = ev_df[ev_df['EV_type'] == ce_i['pev_type']]['usable_battery_size_kWh'].item()
                pev_battery_sizes.append(ev_battery_size)
                start_charge = ce_i['soc_i']*ev_battery_size
                end_charge = ce_i['soc_f']*ev_battery_size
                start_timestep = int(np.floor((ce_i['start_time'] * 3600 - self.start_time_sec)/ ts)) # convert from hours to timestep index
                end_timestep = int(np.floor((ce_i['end_time_chg'] * 3600 - self.start_time_sec) / ts))
                # figure out the minimum you need to charge
                # if your charge session ends beyond the horizon, calc the actual end charge energy for end of horizon
                if end_timestep <0 or start_timestep > n_timesteps:
                    # if the end time is before the simulation start time, then don't include this in the sim
                    continue 
                if end_timestep >= n_timesteps:
                    end_charge = max(0, end_charge - (end_timestep - n_timesteps) * max_evse)
                    end_timestep = n_timesteps - 1
                evse_assets['E0'][se_iter][start_timestep,0]= start_charge
                evse_assets['ET'][se_iter][end_timestep,0]= end_charge #final_energy_i
                #print(f'start_timestep {start_timestep} start_charge {start_charge} vs {sum(sum(evse_assets["E0"][se_iter]))}')
                #print(f'end_timestep: {end_timestep} end_charge: {end_charge} vs {sum(sum(evse_assets["ET"][se_iter]))}')
            evse_assets.loc[evse_assets['SE_id']==se_id,['Emax']] = max(pev_battery_sizes)
            se_iter += 1

        self.evse_assets = evse_assets[evse_assets['Emax']>0] # only load the ones with vehicles at them to save on variables/constraints
        #evse_assets = self.evse_df
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

        # TODO: find E_T for each charging station according to charging vehicles.
        self.nondispatch_assets = nondispatch_assets
        nda_list = nondispatch_assets.to_dict('records')
        for nda in nda_list:
            if isinstance(nda['bus_id'], str):
                nda['number'] = self.bus_df[self.bus_df['name']==nda['bus_id']]['number'][0]
        bus_id_market = 0 # all buses have the same market, so this identifier can be 0

        self.network = Network_3ph(self.feeder_name)
        # the market prices need to be interpolated to match the timestep of the optimization
        timesteps_market = [i*ts_market for i in range(len(prices_import))]
        timesteps_opt = [i*ts for i in range(n_timesteps)]
        prices_import = np.interp(timesteps_opt, timesteps_market, prices_import)
        #print(f'interpolating from {timesteps_market} to {timesteps_opt} with prices {prices_import}')
        market = Market(bus_id_market, prices_export, prices_import, demand_charge, Pmax_market, Pmin_market, ts, hs)
        #print(f'market.T_market is {market.T_market} market.dt is {market.dt_market}')
        self.market = market
        energy_system = EnergySystem(nda_list, self.network, market, ts, hs, ts, hs, evse_assets=self.evse_assets)
        #print(f'energy system T: {energy_system.T}, energy_system.dt {energy_system.dt}, energy_system.T_ems {energy_system.T_ems}')
        self.energy_system = energy_system


    def solve(self, DSS_state_info_dict, federate_time=0):
        # this function solves the control parameters
        ev_control_setpoints = {}
        hs = (self.horizon_sec-self.start_time_sec)
        ts = self.timestep_sec
        start_time = self.start_time_sec
        i_timestep = int(np.floor((federate_time-start_time)/ts))

        # first get the updated grid status
        #voltages = json.loads(h.helicsInputGetString(self.subscriptions[0]))
        #currents = json.loads(h.helicsInputGetString(self.subscriptions[1]))
    
        #non_disp_loads = json.loads(h.helicsInputGetString(self.subscriptions[2]))
        # update the non dispatchable asset loads
        self.nondispatch_assets.loc['Pnet'] = self.nondispatch_assets['Pnet_pred']
        self.nondispatch_assets.loc['Qnet'] = self.nondispatch_assets['Qnet_pred']
        if OpenDSS_message_types.get_basenetloads in DSS_state_info_dict.keys():
            non_disp_loads = DSS_state_info_dict[OpenDSS_message_types.get_basenetloads]
            for busname_i in non_disp_loads.keys():
                self.nondispatch_assets.loc[self.nondispatch_assets['bus_id']==busname_i,['Pnet_pred']] = non_disp_loads[busname_i][0]
                self.nondispatch_assets.loc[self.nondispatch_assets['bus_id']==busname_i,['Qnet_pred']] = non_disp_loads[busname_i][1]
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
        network.N_phases = 1
        network.N_buses = 1
        self.energy_system.network = network
        i_line_unconst_list = []#list(range(network.N_lines))   
        v_bus_unconst_list = [] # no voltage constraints 
        t0 = int(np.floor((federate_time-start_time)/self.energy_system.dt))

        # run the optimization for each bus
        for bus_i in self.bus_df['name']:
            evse_bus_i = self.evse_assets[self.evse_assets['bus_id']==bus_i] # the evse at bus_i
            if len(evse_bus_i)>0:
                if sum(sum(evse_bus_i['ET']))>0:
                    #print(f'running opt for {evse_bus_i}')
                    self.energy_system = EnergySystem(nda_list, self.network, self.market, ts, hs, ts, hs, evse_assets=evse_bus_i)
                    #print(f'energy_system.evse_assets in market_control_block {self.energy_system.evse_assets}')

                    EMS_output = self.energy_system.simulate_network_3phPF('3ph',\
                                                i_unconstrained_lines=\
                                                i_line_unconst_list,\
                                                v_unconstrained_buses=\
                                                v_bus_unconst_list, t0=t0)  
                    i_es = 0
                    for se_id in evse_bus_i['SE_id']:                                              
                        ev_control_setpoints[se_id] = EMS_output['P_EVSE_ems'][i_timestep, i_es]                       
                        i_es = i_es+1
        #pickle.dump(EMS_output, open(join(EMS_path_string, normpath('Month' + str(Case_Month) + '_Day' + str(day) + '_EMS_output_' + str(x) + '.p')), "wb"))    

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EMS Result post-processing and Analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Step 1:  Extract the original EMS results  
        #PF_network_res = EMS_output['PF_network_res']
        #P_import_ems = EMS_output['P_import_ems']
        #P_export_ems = EMS_output['P_export_ems']
        #P_EVSE_ems = EMS_output['P_EVSE_ems']
        ##print(f'P_EVSE_ems line 209 market_control_block: {P_EVSE_ems}')
        ##P_demand_ems = EMS_output['P_demand_ems']  
        #i_es = 0
        #for SE_id in self.evse_assets['SE_id'].values:
        #    # take the first timestep in the horizon opt
        #    ev_control_setpoints[SE_id] = P_EVSE_ems[i_timestep,i_es]#P_import_ems[load_name] - P_export_ems[load_name]
        #   i_es = i_es+1
        #print(f'updating setpoints line 219 market_control_block {ev_control_setpoints}')
        self.control_setpoints = ev_control_setpoints
        return ev_control_setpoints
