
from math import floor
import numpy as np
import pandas as pd

from scipy.optimize import minimize
import sys
import copy

from Caldera_globals import SE_setpoint
from global_aux import Caldera_message_types, OpenDSS_message_types, input_datasets, container_class
from control_templates import typeA_control
from Consensus_inputs import DataInputs
from Consensus_Opt import Hierarchical

from __main__ import *

#from Consensus_Controller import Consensus_Controller.DataInputs
#from Consensus_Controller.Consensus_Controller import Hierarchical


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

class voltwatt_control(typeA_control):

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
        if "ext0001q" in self.datasets_dict[input_datasets.external_strategies]:
            return False
        elif "voltwatt_ld_l2" in self.datasets_dict[input_datasets.external_strategies]:
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
    
        #-------------------------------------
        #      Get Charge Event Data
        #-------------------------------------
        charge_events = self.datasets_dict[input_datasets.SE_group_charge_event_data]
        self.charge_events = charge_events

        #-------------------------------------
        #      Get Supply Equipment Data
        #-------------------------------------
        SE_group_config = self.datasets_dict[input_datasets.SE_group_configuration]  
        #SEid_to_SE_type = self.datasets_dict[input_datasets.SEid_to_SE_type]
        
        SE_group_NODE  = {}
        # Example
        for SE_group in SE_group_config:
            #SE_group_id = SE_group.SE_group_id
            #print('SE_group_id type = ', (SE_group_id))
            
            for SE in SE_group.SEs:
                SE_group_id = int(SE_group.SE_group_id)
                SE_node_id = SE.grid_node_id
                #print('SE_group_id  = ', (SE_group_id))
                #print('SE_node_id  = ', (grid_node_id))
                SE_group_NODE.update({SE_group_id:SE_node_id})
        
        #print('SE_group_NODE = ', SE_group_NODE)
        self.se_group_node = SE_group_NODE  
         
        #SE_group_NODE = {10: '856.2', 20: '826.2'} #Manually mapping SE_group_id with node_id for Sandbox case. 

       
    def log_data(self):
        pass
        
    
    def get_messages_to_request_state_info_from_Caldera(self, next_control_timestep_start_unix_time):
        return_dict = {}
        return_dict[Caldera_message_types.get_active_charge_events_by_SE_groups] = [2,10,20,30]
        #return_dict[Caldera_message_types.get_active_charge_events_by_extCS] = ['ext0001q', 'voltwatt_ld_l2']
        
        # The return value (return_dict) must be a dictionary with Caldera_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict
    
    
    def get_messages_to_request_state_info_from_OpenDSS(self, next_control_timestep_start_unix_time):
        return_dict = {}
        return_dict[OpenDSS_message_types.get_all_node_voltages] = None
        
        # The return value (return_dict) must be a dictionary with OpenDSS_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict
    
    
    def solve(self, next_control_timestep_start_unix_time, Caldera_state_info_dict, DSS_state_info_dict):
        # Caldera_state_info_dict is a dictionary with Caldera_message_types as keys.
        # DSS_state_info_dict is a dictionary with OpenDSS_message_types as keys. 
        
        #current_time = next_control_timestep_start_unix_time
        

        #=================================
        #      Display Node Voltages
        #=================================
        node_voltages = DSS_state_info_dict[OpenDSS_message_types.get_all_node_voltages]
        
        #for (node_id, puV) in node_voltages.items():
            #print('node_id:{}  puV:{}'.format(node_id, puV))
         

        # Step 0: user-defined inputs
        CE_by_SE_groups = Caldera_state_info_dict[Caldera_message_types.get_active_charge_events_by_SE_groups]
        Setpoint_updated = {} #Make sure this is empty at the beginning of calling this federate. 
        
        '''
        Quick fix to retrieve arrival ad departure time of active_CEs from CE_*.csv file
        until messages are received from Caldera federate directly.
        '''
        #Grid teams to modify filename with the address of CE_200/CE-example csv file.
        #filename = 'Inputs\CE_200 pevs.csv'
        #If there is any error in reading csv file, provide absolute path like below:
        #filename = r'C:\Users\KCHAUDHA\Documents\RECHARGE\Ext_Control\Windows-10\Inputs\CE_200 pevs.csv'
        
        
        #df_lookup = pd.read_csv(filename)

        Node_pu = 0.98  # Backup declaration to avoid error
        Vpu_Lower = 0.95
        Vpu_Higher = 1.05
        Adjust_Percent = 0.15  #Customizable charging load adjustment based on grid analysis results
        Voltage_multiplier = 1-Adjust_Percent
       #comm_delay = 0.0050         # 5 ms per communication iteration / communication delay
        comm_delay = 0              #put 0 for RECHARGE
        max_charge = 6.6            # max charging rate kW
        num_charging_stations = 1  # In RECHARGE Context it means Node (DO NOT MODIFY)
        
        
        i= -1
        for (group_id, active_CEs) in CE_by_SE_groups.items():
            i = i + 1
            if bool(active_CEs):
                print('Concensus: Found active charge event')
                len_active_CEs= len(active_CEs)
                #for CE in active_CEs:
                #    len_active_CEs = len_active_CEs + 1
                #    #print('group_id at line 240 = ', group_id)
                #print('length of active_CE = ', len_active_CEs)

                SE_id_CE = np.zeros((len(CE_by_SE_groups), len_active_CEs))
                charge_event_id = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                #event_nodes = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                now_unix_time = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                now_soc = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                now_charge_energy_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                energy_of_complete_charge_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                remaining_charge_energy_ackWh = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                arrival_time = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                departure_time = np.zeros((len(CE_by_SE_groups),len_active_CEs))
                
 
                j = 0
                for CE in active_CEs:
                    
                    SE_id_CE[i,j] = CE.SE_id
                    #print(SE_id_CE)
                    charge_event_id[i,j] = CE.charge_event_id
                    #print(charge_event_id)
                    now_unix_time[i,j] = CE.now_unix_time/3600
                    now_soc[i,j] = CE.now_soc
                    now_charge_energy_ackWh[i,j] = CE.now_charge_energy_ackWh
                    energy_of_complete_charge_ackWh[i,j] = CE.energy_of_complete_charge_ackWh
                    remaining_charge_energy_ackWh[i,j] = energy_of_complete_charge_ackWh[i,j] - now_charge_energy_ackWh[i,j]
                    arrival_time[i,j] = now_unix_time[i,j]
                    #departure_time[i,j] = df_lookup._get_value((charge_event_id[i,j]-1),'end_time_prk')
                    for m in range(len(self.charge_events)):
                        for k in range(len(self.charge_events[m].charge_events)):
                            if self.charge_events[m].charge_events[k].charge_event_id == CE.charge_event_id:
                                departure_time[i,j] = self.charge_events[m].charge_events[k].departure_unix_time/3600   
                   
                    j = j + 1
                    

            else:
                print('Concensus: No active charge event')
                
            PQ_setpoints = []
            if bool(active_CEs):
                group_id_Parse = group_id
                charge_event_id_Parse = charge_event_id[i,:]
                SE_id_Parse = SE_id_CE[i,:]
    
                arrival_times_Parse = arrival_time[i,:]
                #print('arrival times', arrival_time)
                departure_times_Parse = departure_time[i,:]
                #print('departure_times', departure_time)
                charging_times_Parse = departure_times_Parse - arrival_times_Parse
                energy_used_Parse = remaining_charge_energy_ackWh[i,:]
                SE_group_NODE_Parse = self.se_group_node[group_id]
                for (node_id, puV) in node_voltages.items():
                    if  node_id == SE_group_NODE_Parse:
                        Node_pu = puV
                SetPoint_Dict = {} 
                # Step 1: load data
                # Active CEs for each node passed.  all times are in hours
                data = DataInputs(comm_delay,num_charging_stations,Node_pu,Vpu_Lower,Vpu_Higher,Voltage_multiplier,CE_by_SE_groups,charge_event_id_Parse,SE_id_Parse,arrival_times_Parse,departure_times_Parse,charging_times_Parse,energy_used_Parse,group_id_Parse,SetPoint_Dict)
                # Step 2: Run opt code
                opt_hierarchical = Hierarchical(data,max_charge,next_control_timestep_start_unix_time) 
                SetPoint_Dict = opt_hierarchical.solve_real_time()
                #Step 3: Store results
                for key, value in SetPoint_Dict.items():
                    Setpoint_updated[key] = value
                    print(f'powerlimit id: {key} set to {value}')
                    setpoint = SE_setpoint()
                    setpoint.SE_id = int(key)
                    setpoint.PkW = value
                    setpoint.QkVAR = 0
                    PQ_setpoints.append(setpoint)
        
        #======================================================
        #      Open-DSS Feedback Control algorithm ends here
        #======================================================         
            
       
        #=================================
        #      Control PEV Charging
        #=================================        
        #next_control_timestep_start_hrs = (next_control_timestep_start_unix_time/3600)        
        #clock_mins = 60*(next_control_timestep_start_hrs - floor(next_control_timestep_start_hrs))
        
            
        #-----------------------------
        
        Caldera_control_info_dict = {}
        if len(PQ_setpoints) > 0:
            Caldera_control_info_dict[Caldera_message_types.set_pev_charging_PQ] = PQ_setpoints
        
        DSS_control_info_dict = {}
        
        # Caldera_control_info_dict must be a dictionary with Caldera_message_types as keys.
        # DSS_control_info_dict must be a dictionary with OpenDSS_message_types as keys.
        # If either value has nothing to return, return an empty dictionary.
        sys.stdout.flush()
        return (Caldera_control_info_dict, DSS_control_info_dict)

