
from math import floor
from tokenize import group

from Caldera_globals import SE_setpoint
from global_aux import Caldera_message_types, OpenDSS_message_types, input_datasets, container_class
from control_templates import typeB_control

import numpy as np
import pandas as pd
from btm_control import BTM_Control
import time
import datetime, calendar

import sys

'''
Notice: This computer software was prepared by Alliance for Sustainable Energy, LLC, 
hereinafter the Contractor, under Contract DE-AC36-08GO28308 with the Department of 
Energy (DOE). All rights in the computer software are reserved by DOE on behalf of 
the United States Government and the Contractor as provided in the Contract. You are
authorized to use this computer software for Governmental purposes but it is not to 
be released or distributed to the public. NEITHER THE GOVERNMENT NOR THE CONTRACTOR 
MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS 
SOFTWARE. This notice including this sentence must appear on any copies of this computer
software.
'''


class btms_control(typeB_control):
    
    def __init__(self, base_dir, simulation_time_constraints, input_se_csv='inputs/SE_.csv', ce_ext_strategy='ext0002'):        
        super().__init__(base_dir, simulation_time_constraints)
        
        self.ext_strategy = ce_ext_strategy
        self.control_timestep_min = 15
        self.request_state_lead_time_min = 10.1
        self.send_control_info_lead_time_min = (simulation_time_constraints.grid_timestep_sec + 0.5)/60
        se_df = pd.read_csv(input_se_csv)
        se_df =se_df[['SE_id','node_id', 'SE_group']]
        self.se_group_by_bus = se_df.set_index('SE_id').to_dict()
        self.se_groups = list(set(se_df['SE_group']))
        self.evse_df = pd.read_csv('EVSE_inputs.csv')
        
    
    def get_input_dataset_enum_list(self):
        return [input_datasets.SE_group_configuration, input_datasets.SE_group_charge_event_data, input_datasets.SEid_to_SE_type, input_datasets.external_strategies]


    def load_input_datasets(self, datasets_dict):
        # datasets_dict is a dictionary with input_datasets as keys.
        self.datasets_dict = datasets_dict
    
    
    def terminate_this_federate(self):
        #print(self.datasets_dict[input_datasets.external_strategies])
        if self.ext_strategy in self.datasets_dict[input_datasets.external_strategies]:
            print(f'running with btms_ld_l2 federate')
            return False
        #elif "ext0001q" in self.datasets_dict[input_datasets.external_strategies]:
        #    return False

        return True
    
    
    def initialize(self):
        #-------------------------------------
        #    Calculate Timing Parameters
        #-------------------------------------        
        X = container_class()
        X.control_timestep_min = self.control_timestep_min
        X.request_state_lead_time_min = self.request_state_lead_time_min
        X.send_control_info_lead_time_min = self.send_control_info_lead_time_min
        #self._calculate_timing_parameters(X, self.__class__.__name__)
    
        deptTime_lookup = {}
        #-------------------------------------
        #      Get Charge Event Data
        #-------------------------------------
        charge_events = self.datasets_dict[input_datasets.SE_group_charge_event_data]
        self.charge_events = charge_events
        for ces in self.charge_events:
            for ce in ces.charge_events:
                ID = ce.charge_event_id
                deptTime = ce.departure_unix_time
                deptTime_lookup[ID]= deptTime
                
        self.deptTime_lookup = deptTime_lookup

        #-------------------------------------
        for CE_group in charge_events:
            SE_group_id = CE_group.SE_group_id
            
            for CE in CE_group.charge_events:
                id = CE.charge_event_id
                SE_id = CE.SE_id
                arrival_time = CE.arrival_unix_time
                arrival_SOC = CE.arrival_SOC
                #...

        #-------------------------------------
        #      Get Supply Equipment Data
        #-------------------------------------
        SE_group_config = self.datasets_dict[input_datasets.SE_group_configuration]  
        self.se_group_config =SE_group_config      
        SEid_to_SE_type = self.datasets_dict[input_datasets.SEid_to_SE_type]
        self.SEid_to_SE_type = SEid_to_SE_type
  
       
    def log_data(self):
        pass

    
    def get_messages_to_request_state_info_from_Caldera(self, next_control_timestep_start_unix_time):
        return_dict = {}
        #print(f'groups by bus: {self.se_groups}')
        return_dict[Caldera_message_types.get_active_charge_events_by_SE_groups] = self.se_groups #Grid teams to update this
        return_dict[Caldera_message_types.get_active_charge_events_by_extCS] = [self.ext_strategy] #['ext0003']


        # The return value (return_dict) must be a dictionary with Caldera_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict
    
    
    #def get_messages_to_request_state_info_from_Caldera(self, next_control_timestep_start_unix_time):
    #    return_dict = {}
    #    
    #    # The return value (return_dict) must be a dictionary with Caldera_message_types as keys.
    #    # If there is nothing to return, return an empty dictionary.
    #    return return_dict
    
    
    def get_messages_to_request_state_info_from_OpenDSS(self, next_control_timestep_start_unix_time):
        return_dict = {}
        return_dict[OpenDSS_message_types.get_all_DER] = None
        return_dict[OpenDSS_message_types.get_basenetloads] = None
        #return_dict[OpenDSS_message_types.get_all_node_voltages] = None      
        
        # The return value (return_dict) must be a dictionary with OpenDSS_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict
    
   
 
    def solve(self, next_control_timestep_start_unix_time, Caldera_control_info_dict, DSS_state_info_dict):
        # this solves a timeseries optimization of energy storage dispatch with the dwell period as the horizon
        # energy storage is used to attempt to flatten total net load over the horizon period
        # 
        # Caldera_control_info_dict is a dictionary with Caldera_message_types as keys.
        # DSS_state_info_dict is a dictionary with OpenDSS_message_types as keys.                 
        
        
        #======================================================
        #      Get data from OpenDSS
        #======================================================

        if OpenDSS_message_types.get_all_DER in DSS_state_info_dict.keys():
            DER_data = DSS_state_info_dict[OpenDSS_message_types.get_all_DER]
            #print(DER_data)
            Storage_SOC = DER_data["storage_SOC"]
            Storage_Capacity = DER_data["storage_cap_kwh"]
            #Net_load = DER_data["Net_load"]
            storage_buses = DER_data["bus_name"]
            #print(f'DER_data stor soc: {DER_data["storage_SOC"]}')
            Net_load = DSS_state_info_dict[OpenDSS_message_types.get_basenetloads]
        else:
            print(f'no der info at timestep {next_control_timestep_start_unix_time}')
            DSS_control_info_dict = {}
            DER_data = {}
            return (Caldera_control_info_dict, DSS_control_info_dict, DER_data)

        #=================================
        #      Define BTM controller
        #=================================
        btm_control = BTM_Control(time_step_mins=5, ess_size=Storage_Capacity, max_power_ess=40, min_power_ess=-40, 
                                  max_power_l2=17.66, min_power_l2=1.5, time_horizon=1)
        # time_horizon is in hours
        #btm_control = BTM_Control(time_step_mins=5, ess_size=40, max_power_ess=40, min_power_ess=-40, 
        #                          max_power_l2=6.6, min_power_l2=1.5, time_horizon=5)
        
        Emin = list(np.zeros(int(btm_control.time_horizon*60/btm_control.time_step_mins)))
        Emax = list(np.zeros(int(btm_control.time_horizon*60/btm_control.time_step_mins)))
        '''
        for (node_id, puV) in node_voltages.items():
            print('node_id:{}  puV:{}'.format(node_id, puV))
        '''

        CE_by_SE_groups = Caldera_control_info_dict[Caldera_message_types.get_active_charge_events_by_SE_groups]
        CE_by_ext = Caldera_control_info_dict[Caldera_message_types.get_active_charge_events_by_extCS]
        #print(CE_by_ext)
        # make a list of active SE_ids:
        active_SEs = []
        for CE in CE_by_ext[self.ext_strategy]:
            active_SEs.append(CE.SE_id)

        Setpoint_updated = {} #Make sure this is empty at the beginning of calling this federate. 
        storage_powers_setpoints={} #Make sure this is empty at the beginning of calling this federate. 

          
        
        i= -1     
        for (group_id, active_CEs) in CE_by_SE_groups.items():
            i = i + 1   #i for interating trhough group node
            i_storage = -1
            storage_at_bus = False
            
            if bool(active_CEs):
                print('BTM: Found active charge event')
                len_active_CEs= 0
                # find storage at buses with actie charge events
                
                storages_involved = []
                bus_load_without_ev_ess = []
                #ces_with_storage = []
                active_buses = []
                i_ce = 0
                for CE in active_CEs:
                    if CE.SE_id in active_SEs:
                        len_active_CEs = len_active_CEs + 1
                        #print(f'btms CE.SE_id: {CE.SE_id}')
                        ce_bus = self.se_group_by_bus['node_id'][CE.SE_id]
                        if ce_bus in storage_buses:
                            i_storage = storage_buses.index(ce_bus)
                            storages_involved.append(i_storage)
                            bus_load_without_ev_ess.append(Net_load[ce_bus])
                            #ces_with_storage.append(i_ce)
                            active_buses.append(ce_bus)
                            storage_at_bus = True
                        else:
                            bus_load_without_ev_ess = list(np.zeros(int(btm_control.time_horizon*60/btm_control.time_step_mins)))
                        i_ce = i_ce+1
                # remove duplicates
                storages_involved = list(set(storages_involved))
                
                active_buses = list(set(active_buses))
                    #print('group_id at line 240 = ', group_id)
                #print('length of active_CE = ', len_active_CEs)

                SE_id_CE = np.zeros((len(CE_by_SE_groups), len_active_CEs),dtype=np.int64)
                #charge_event_id = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.int64)
                #event_nodes = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.int64)
                now_unix_time = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                now_soc = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64) # this is vehicle SOC
                now_charge_energy_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                energy_of_complete_charge_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                remaining_charge_energy_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                arrival_time = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                departure_time = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
 
                j = 0  #i for interating through active_CE
                for CE in active_CEs:
                    if CE.SE_id in active_SEs:
                        SE_id_CE[i,j] = CE.SE_id
                        #print(SE_id_CE)
                        #charge_event_id[i,j] = CE.charge_event_id
                        #print(charge_event_id)
                        now_unix_time[i,j] = CE.now_unix_time
                        #now_soc[i,j] = CE.now_soc
                        now_charge_energy_ackWh[i,j] = CE.now_charge_energy_ackWh
                        energy_of_complete_charge_ackWh[i,j] = CE.energy_of_complete_charge_ackWh
                        #print(f'now_charge: {CE.now_charge_energy_ackWh}, energy_complete: {CE.energy_of_complete_charge_ackWh}')
                        remaining_charge_energy_ackWh[i,j] = energy_of_complete_charge_ackWh[i,j] - now_charge_energy_ackWh[i,j]
                        arrival_time[i,j] = now_unix_time[i,j]
                        #departure_time[i,j] = df_lookup._get_value((charge_event_id[i,j]-1),'end_time_prk')
                        departure_time[i,j] = self.deptTime_lookup[CE.charge_event_id]
                        j = j + 1
                    
                btm_control.num_of_active_evse = j
            else:
                print('BTM: No active charge event')
            
            
            #Step 0 : All the input data required for BTM control
        
            if len_active_CEs>0:
                group_id_Parse = group_id
                #charge_event_id_Parse = charge_event_id[i,:]
                SE_id_Parse = SE_id_CE[i,:]
                se_type = self.SEid_to_SE_type[SE_id_Parse[0]]
                btm_control.max_power_l2 = self.evse_df[self.evse_df['EVSE_type']==se_type]['AC/DC_power_limit_kW'].values[-1]
 
                now_unix_time_Parse = now_unix_time[i,:]
                arrival_times_Parse = arrival_time[i,:]
                #print('arrival times', arrival_times_Parse)
                departure_times_Parse = departure_time[i,:]
                #print('departure_times', departure_times_Parse)
                charging_times_Parse = departure_times_Parse - arrival_times_Parse
                #print('charging times', charging_times_Parse)
                energy_used_Parse = remaining_charge_energy_ackWh[i,:]
                #SE_group_NODE_Parse = SE_group_NODE[group_id]
                
                Storage_SOC_Parse = Storage_SOC
                Storage_Capacity_Parse = Storage_Capacity
        
        
                #SetPoint_Dict = {} #for EV power setpoint
                Setpoint_storage = {} #for storage power setpoint
                #Step 1 : BTM control
                
                #=====================================================
                # Calculate Emin/Emax for the active EVSE's
                #=====================================================
                for x1, x2, x3 in zip(now_unix_time_Parse, departure_times_Parse, energy_used_Parse):
                    [Emin_, Emax_] = btm_control.energy_constraints_evse(x1, x2, x3)
                    Emin = [sum(x) for x in zip(Emin, Emin_)]         
                    Emax = [sum(x) for x in zip(Emax, Emax_)]
                
                #=====================================================
                # Calculate Emin/Emax for ESS
                #=====================================================
                tm = time.localtime(now_unix_time[i,0])
                ess_target_time = (datetime.datetime(tm.tm_year, tm.tm_mon, tm.tm_mday, 0, 0, 0, 0) \
                                  + datetime.timedelta(hours=24)).timestamp()  ## midnight
                if (ess_target_time - now_unix_time[i,0])/60/60 < btm_control.time_horizon:
                    ess_target_time = (datetime.datetime(tm.tm_year, tm.tm_mon, tm.tm_mday, 0, 0, 0, 0) \
                                      + datetime.timedelta(hours=48)).timestamp()

                available_energy_ess = 0
                i_storage_soc = 0
                for i_storage in storages_involved:   
                    remaining_energy_ess = (1 - Storage_SOC_Parse[i_storage])*Storage_Capacity_Parse[i_storage]
                    available_energy_ess = (Storage_SOC_Parse[i_storage] - btm_control.ess_soc_min)*Storage_Capacity_Parse[i_storage]
                    Emin_, Emax_ = btm_control.energy_constraints_ess(now_unix_time[i,0], ess_target_time,  
                                                                    remaining_energy_ess, available_energy_ess)
                    i_storage_soc = Storage_SOC_Parse[i_storage]
                    Emin = [sum(x) for x in zip(Emin, Emin_)] 
                    Emax = [sum(x) for x in zip(Emax, Emax_)]
                
                #=====================================================
                # Solve optimization
                #=====================================================
                print(f'there are {len(active_CEs)} charge events and {len(storages_involved)} storages \n Emin is {Emin} \n Emax is {Emax}')
                #p = non_pev_loads_forecast - pv_powers_forecast
                p = bus_load_without_ev_ess
                print(f'power without ev and ess is {p}')
                x0 = np.zeros(len(p)) #np.ones(len(p))
                results = btm_control.solve_optimization(x0, p, Emin, Emax) 

                #print(f'opt results: {results.x}')
                
                #=====================================================
                # Allocate setpoint
                #=====================================================
                Setpoint_storage, SetPoint_list = btm_control.allocate_setpoint(results.x[0], i_storage_soc, available_energy_ess, departure_times_Parse, energy_used_Parse)
                
                #Update SetPoint_Dict within the control with key being SE_id and value being Power
                #Update Setpoint_storage within the control with key being node_ID and value being Power
                
                #Step 2: Store results to pass to Caldera
                j = 0  #i for interating through active_CE
                for SEs in SE_id_Parse:
                    Setpoint_updated[SEs] = SetPoint_list[j]
                    j = j + 1
        
                #Step 3: Store results to pass to OpenDSS and update dataframe tracking soc
                i_setp = 0
                for i_storage in storages_involved:
                    ce_bus = active_buses[i_setp]
                    if not ce_bus in storage_powers_setpoints.keys():
                        storage_powers_setpoints[ce_bus] = []
                    storage_powers_setpoints[ce_bus].append(Setpoint_storage[i_setp])

                    #Step 3.5: if the storage isn't passed to OpenDSS save it in the dataframe
                    ess_energy_used = Setpoint_storage[i_setp]*self.control_timestep_min/60
                    DER_data['storage_SOC'][i_storage] = DER_data['storage_SOC'][i_storage] - ess_energy_used/DER_data['storage_cap'][i_storage]
                    i_setp = i_setp + 1


      
        #=================================
        #      Control PEV Charging
        #=================================        
        next_control_timestep_start_hrs = (next_control_timestep_start_unix_time/3600)        
        clock_mins = 60*(next_control_timestep_start_hrs - floor(next_control_timestep_start_hrs))
        
        
        #-----------------------------
        PQ_setpoints = []        
        
        if bool(Setpoint_updated):
            for (group_id, active_CEs) in CE_by_SE_groups.items():
                for CE in active_CEs:
                    if CE.SE_id in active_SEs:
                        remaining_charge_energy_ackWh = CE.energy_of_complete_charge_ackWh - CE.now_charge_energy_ackWh
                        if 0 < remaining_charge_energy_ackWh:
                            X = SE_setpoint()
                            X.SE_id = CE.SE_id
                            for (SE_id, setpoint) in Setpoint_updated.items():
                                if  SE_id ==  X.SE_id:
                                    X.PkW = setpoint
                                    X.QkVAR = 0
                                    PQ_setpoints.append(X)

            
        #-----------------------------
        
        Caldera_control_info_dict = {}
        #print(f'pq pev setpoints {PQ_setpoints}')
        if len(PQ_setpoints) > 0:
            Caldera_control_info_dict[Caldera_message_types.set_pev_charging_PQ] = PQ_setpoints
            print(f'first in list {PQ_setpoints[0].PkW}')
        
 
        DSS_control_info_dict = {}
        
        if len(storage_powers_setpoints) > 0:
            DSS_control_info_dict[OpenDSS_message_types.set_storage_power] = storage_powers_setpoints
        # Caldera_control_info_dict must be a dictionary with Caldera_message_types as keys.
        # DSS_control_info_dict must be a dictionary with OpenDSS_message_types as keys.
        # If either value has nothing to return, return an empty dictionary.
        sys.stdout.flush()
        return (Caldera_control_info_dict, DSS_control_info_dict, DER_data)
        
               
       

