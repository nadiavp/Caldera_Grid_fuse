
from math import floor
from tokenize import group

from Caldera_globals import SE_setpoint
from global_aux import Caldera_message_types, OpenDSS_message_types, input_datasets, container_class
from control_templates import typeB_control

import numpy as np
from btm_control import BTM_Control
import time
import datetime, calendar

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
    
    def __init__(self, base_dir, simulation_time_constraints):        
        super().__init__(base_dir, simulation_time_constraints)
        
        self.control_timestep_min = 15
        self.request_state_lead_time_min = 10.1
        self.send_control_info_lead_time_min = (simulation_time_constraints.grid_timestep_sec + 0.5)/60
        
    
    def get_input_dataset_enum_list(self):
        return [input_datasets.SE_group_configuration, input_datasets.SE_group_charge_event_data, input_datasets.SEid_to_SE_type, input_datasets.external_strategies]


    def load_input_datasets(self, datasets_dict):
        # datasets_dict is a dictionary with input_datasets as keys.
        self.datasets_dict = datasets_dict
    
    
    def terminate_this_federate(self):
        print(self.datasets_dict[input_datasets.external_strategies])
        if "ext_btms_ld_l2" in self.datasets_dict[input_datasets.external_strategies]:
            print(f'running with btms_ld_l2 federate')
            return False
        elif "ext0002" in self.datasets_dict[input_datasets.external_strategies]:
            return False

        return True
    
    
    def initialize(self):
        #-------------------------------------
        #    Calculate Timing Parameters
        #-------------------------------------        
        X = container_class()
        X.control_timestep_min = self.control_timestep_min
        X.request_state_lead_time_min = self.request_state_lead_time_min
        X.send_control_info_lead_time_min = self.send_control_info_lead_time_min
        self._calculate_timing_parameters(X, self.__class__.__name__)
    
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
        SEid_to_SE_type = self.datasets_dict[input_datasets.SEid_to_SE_type]
  
       
    def log_data(self):
        pass

    
    def get_messages_to_request_state_info_from_Caldera(self, next_control_timestep_start_unix_time):
        return_dict = {}
        #return_dict[Caldera_message_types.get_active_charge_events_by_SE_groups] = [10] #Grid teams to update this
        return_dict[Caldera_message_types.get_active_charge_events_by_extCS] = ['ext0002', 'ext_btsm_ld_l2']

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
        return_dict[OpenDSS_message_types.get_all_node_voltages] = None      
        
        # The return value (return_dict) must be a dictionary with OpenDSS_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict
    
   
 
    def solve(self, next_control_timestep_start_unix_time, Caldera_control_info_dict, DSS_state_info_dict):
        # Caldera_control_info_dict is a dictionary with Caldera_message_types as keys.
        # DSS_state_info_dict is a dictionary with OpenDSS_message_types as keys.                 
        
        #=================================
        #      Define BTM controller
        #=================================
        btm_control = BTM_Control(time_step_mins=5, ess_size=40, max_power_ess=40, min_power_ess=-40, 
                                  max_power_l2=6.6, min_power_l2=1.5, time_horizon=5)
        
        Emin = list(np.zeros(int(btm_control.time_horizon*60/btm_control.time_step_mins)))
        Emax = list(np.zeros(int(btm_control.time_horizon*60/btm_control.time_step_mins)))
        
        
        #======================================================
        #      Get data from OpenDSS
        #======================================================

        if OpenDSS_message_types.get_all_DER in DSS_state_info_dict.keys():
            DER_data = DSS_state_info_dict[OpenDSS_message_types.get_all_DER]
            Net_load = DER_data["Net_load"]
            Storage_SOC = DER_data["storage_SOC"]
            Storage_Capacity = DER_data["storage_cap"]
        else:
            print(f'no der info at timestep {next_control_timestep_start_unix_time}')
            print(f'DSS_state_info keys: {DSS_state_info_dict.keys()}')
            DSS_control_info_dict = {}
            return (Caldera_control_info_dict, DSS_control_info_dict)

        
        '''
        for (node_id, puV) in node_voltages.items():
            print('node_id:{}  puV:{}'.format(node_id, puV))
        '''

        CE_by_SE_groups = Caldera_control_info_dict[Caldera_message_types.get_active_charge_events_by_SE_groups]
        Setpoint_updated = {} #Make sure this is empty at the beginning of calling this federate. 
        storage_powers_setpoints={} #Make sure this is empty at the beginning of calling this federate. 

          
        
        i= -1     
        for (group_id, active_CEs) in CE_by_SE_groups.items():
            i = i + 1   #i for interating trhough group node
            if bool(active_CEs):
                print('BTM: Found active charge event')
                len_active_CEs= 0
                for CE in active_CEs:
                    len_active_CEs = len_active_CEs + 1
                    #print('group_id at line 240 = ', group_id)
                #print('length of active_CE = ', len_active_CEs)

                SE_id_CE = np.zeros((len(CE_by_SE_groups), len_active_CEs),dtype=np.int64)
                charge_event_id = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.int64)
                event_nodes = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.int64)
                now_unix_time = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                now_soc = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                now_charge_energy_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                energy_of_complete_charge_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                remaining_charge_energy_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                arrival_time = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                departure_time = np.zeros((len(CE_by_SE_groups),len_active_CEs),dtype=np.float64)
                
 
                j = 0  #i for interating through active_CE
                for CE in active_CEs:
                    
                    SE_id_CE[i,j] = CE.SE_id
                    #print(SE_id_CE)
                    charge_event_id[i,j] = CE.charge_event_id
                    #print(charge_event_id)
                    now_unix_time[i,j] = CE.now_unix_time
                    now_soc[i,j] = CE.now_soc
                    now_charge_energy_ackWh[i,j] = CE.now_charge_energy_ackWh
                    energy_of_complete_charge_ackWh[i,j] = CE.energy_of_complete_charge_ackWh
                    remaining_charge_energy_ackWh[i,j] = energy_of_complete_charge_ackWh[i,j] - now_charge_energy_ackWh[i,j]
                    arrival_time[i,j] = now_unix_time[i,j]
                    #departure_time[i,j] = df_lookup._get_value((charge_event_id[i,j]-1),'end_time_prk')
                    departure_time[i,j] = self.deptTime_lookup[CE.charge_event_id]
                    j = j + 1
                    
                btm_control.num_of_active_evse = j
            else:
                print('BTM: No active charge event')
            
            
            #Step 0 : All the input data required for BTM control
        
            if bool(active_CEs):
                group_id_Parse = group_id
                charge_event_id_Parse = charge_event_id[i,:]
                SE_id_Parse = SE_id_CE[i,:]
 
                now_unix_time_Parse = now_unix_time[i,:]
                arrival_times_Parse = arrival_time[i,:]
                #print('arrival times', self.arrival_times)
                departure_times_Parse = departure_time[i,:]
                #print('departure_times', self.departure_times)
                charging_times_Parse = departure_times_Parse - arrival_times_Parse
                energy_used_Parse = remaining_charge_energy_ackWh[i,:]
                #SE_group_NODE_Parse = SE_group_NODE[group_id]
                
                Storage_SOC_Parse = Storage_SOC
                Storage_Capacity_Parse = Storage_Capacity
        
        
                SetPoint_Dict = {} #for EV power setpoint
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
                    
                remaining_energy_ess = (1 - Storage_SOC_Parse[i])*Storage_Capacity_Parse[i]
                available_energy_ess = (Storage_SOC_Parse[i] - btm_control.ess_soc_min)*Storage_Capacity_Parse[i]
                Emin_, Emax_ = btm_control.energy_constraints_ess(now_unix_time[i,0], ess_target_time,  
                                                                    remaining_energy_ess, available_energy_ess)
                Emin = [sum(x) for x in zip(Emin, Emin_)]         
                Emax = [sum(x) for x in zip(Emax, Emax_)]
                
                #=====================================================
                # Solve optimization
                #=====================================================
                #p = non_pev_loads_forecast - pv_powers_forecast
                p = Net_load[i]
                x0 = np.ones(len(p))
                results = btm_control.solve_optimization(x0, p, Emin, Emax) 
                
                #=====================================================
                # Allocate setpoint
                #=====================================================
                Setpoint_storage, SetPoint_list = btm_control.allocate_setpoint(results.x[0], Storage_SOC_Parse[i], available_energy_ess, departure_times_Parse, energy_used_Parse)
                
                #Update SetPoint_Dict within the control with key being SE_id and value being Power
                #Update Setpoint_storage within the control with key being node_ID and value being Power
                
                #Step 2: Store results to pass to Caldera
                j = 0  #i for interating through active_CE
                for SEs in SE_id_Parse:
                    Setpoint_updated[SEs] = SetPoint_list[j]
                    j = j + 1
        
                #Step 3: Store results to pass to OpenDSS
                if group_id_Parse not in storage_powers_setpoints.keys():
                    storage_powers_setpoints[group_id_Parse] = []
                storage_powers_setpoints[group_id_Parse].append(Setpoint_storage)


      
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
        if len(PQ_setpoints) > 0:
            Caldera_control_info_dict[Caldera_message_types.set_pev_charging_PQ] = PQ_setpoints
        
 
        DSS_control_info_dict = {}
        
        if len(storage_powers_setpoints) > 0:
            DSS_control_info_dict[OpenDSS_message_types.set_storage_power] = storage_powers_setpoints
        # Caldera_control_info_dict must be a dictionary with Caldera_message_types as keys.
        # DSS_control_info_dict must be a dictionary with OpenDSS_message_types as keys.
        # If either value has nothing to return, return an empty dictionary.
        return (Caldera_control_info_dict, DSS_control_info_dict)
        
               
       

