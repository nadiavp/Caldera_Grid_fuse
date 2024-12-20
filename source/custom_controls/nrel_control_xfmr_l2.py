import os 
from os.path import normpath, join
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import pickle
import warnings
warnings.filterwarnings('ignore')

from Caldera_globals import SE_setpoint
from global_aux import Caldera_message_types, OpenDSS_message_types, input_datasets, container_class
import opendssdirect as dss


class transformer_control():
    def __init__(self, base_dir, simulation_time_constraints, input_se_csv='inputs/SE_.csv', 
        name='equal_sharing', helics_config_path='', feeder_name='ieee34', input_ce_csv='inputs/CE_.csv', ce_ext_strategy='ext0004', se_group=[10]):
        # name options: ['first_come_first_served', 'fcfs_with_minimum', 'equal_sharing']
        self.name=name
        self.se_file = input_se_csv
        self.ce_file = input_ce_csv
        self.se_group = se_group # this is a code which matches the code in the SE input file. It determines which SE locations get this control
        self.ce_ext_strategy = ce_ext_strategy# this is a code which matches the code in the CE input file and determines which vehicles get this control
        self.timestep_sec = simulation_time_constraints.grid_timestep_sec
        self.horizon_sec = simulation_time_constraints.end_simulation_unix_time
        self.start_time_sec = simulation_time_constraints.start_simulation_unix_time
        self.simulation_time_constraints = simulation_time_constraints
        self.time = -1
        self.helics_config_path = helics_config_path
        self.fed = None
        self.publications = []
        self.subscriptions = []
        # transformer ratings sorted by node name are needed
        self.xfmr_ratings_by_node_name = []



    def get_input_dataset_enum_list(self):
        return [input_datasets.SE_group_configuration, input_datasets.SE_group_charge_event_data, input_datasets.SEid_to_SE_type, input_datasets.external_strategies]


    def load_input_datasets(self, datasets_dict):
        # datasets_dict is a dictionary with input_datasets as keys.
        self.datasets_dict = datasets_dict
    
    
    def terminate_this_federate(self):
        if self.ce_ext_strategy in self.datasets_dict[input_datasets.external_strategies]:
            print(f'running with transformer overload mitigation federate')
            return False
        return True

    def log_data(self):
        pass

    
    def get_messages_to_request_state_info_from_Caldera(self, next_control_timestep_start_unix_time):
        return_dict = {}
        return_dict[Caldera_message_types.get_active_charge_events_by_extCS] = [self.ce_ext_strategy]
        # The return value (return_dict) must be a dictionary with Caldera_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict
    
    
    def get_messages_to_request_state_info_from_OpenDSS(self, next_control_timestep_start_unix_time):
        return_dict = {}
        return_dict[OpenDSS_message_types.get_all_DER] = None
        return_dict[OpenDSS_message_types.get_basenetloads] = None
        #return_dict[OpenDSS_message_types.get_all_node_voltages] = None      
        
        # The return value (return_dict) must be a dictionary with OpenDSS_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict

    def initialize(self, ):
        self.management_system = ChargingManagementSystem(capacity_series_1min['Capacity'], allocation_method=method)
        # get a list of nodes with SE
        SE_df = pd.read_csv(self.se_file)
        SE_df['trns_node'] = SE_df['node_id'].split('_')[0]
        #get transformer ratings
        # first compile the feeder
        dss.Command("clear")
        dss.Command(f'Redirect "{self.opendss_dir}"')
        dss.Command("solve")
        # next get the transformer names and names of buses they are connected to
        trns_names = dss.Transformers.AllNames()
        trns_kva = {}
        for trns in trns_names:
            dss.Transformers.Name(trns)
            node_names = dss.CktElement.BusNames()
            trns_kva[trns] = dss.Transformers.kVA()
            trns_nodes[trns] = node_names
        # map the nodes to SE_id and create a transformer by SE_id list
        # there could be multiple SE_ids that go to one transformer
            for node_name in node_names:
                SE_id = SE_df[SE_df['node_id']==node_name]['SE_id']
                if len(SE_id)>0:
                    trns_by_seid[SE_id] = trns
        self.trns_kva = trns_kva
        self.trns_by_seid = trns_by_seid
        self.nodes_by_trns = trns_nodes

    def solve(self, federate_time, Caldera_state_info_dict, DSS_state_info_dict):#df_unique_vehicles_per_transformer, CEs_feeder_day, tf_capacity_available_KW):
        ev_control_setpoints = {}
        PQ_setpoints = []
        
        # Creating a capacity time series
        Sim_start_time = self.start_time
        Sim_end_time = self.horizon
        time_index = pd.date_range(start=self.start_time_sec, end=self.horizon_sec, freq='T')
        all_ev_power_profiles = pd.DataFrame(index=time_index)
        all_ev_energy_profiles = pd.DataFrame(index=time_index)
        all_charging_events_evaluation = []
        # get CEs from caldera into charging stations
        CE_by_extCS = Caldera_state_info_dict[Caldera_message_types.get_active_charge_events_by_extCS]
        if self.ce_ext_strategy in CE_by_extCS.keys():
            active_CEs = CE_by_extCS[self.ce_ext_strategy]
        else: 
            ev_control_setpoints = {}
            Caldera_control_info_dict, DSS_state_info_dict, ev_control_setpoints 

        # get baseload by node and transformer info from opendss
        baseloads = DSS_state_info_dict[OpenDSS_message_types.get_basenetloads]
        # TODO: make method that get get transformer info
        unique_vehicles_per_transformer = {}
        for CE in active_CEs:
            SE_id = CE.SE_id
            trns_name = self.trns_by_seid[SE_id]
            if not trns_name in unique_vehicles_per_transformer:
                unique_vehicles_per_transformer[trns_name] = [SE_id]
            else:
                unique_vehicles_per_transformer[trns_name].append(SE_id)
            

        for transformer_id in unique_vehicles_per_transformer.keys():
            # one management system for each group of EVSE
            self.management_system = ChargingManagementSystem(capacity_series_1min['Capacity'], allocation_method=method)
            transformer_id_str = str(transformer_id)
            se_ids = unique_vehicles_per_transformer[transformer_id]#row['Unique_Vehicles']
            
            print(f"Processing transformer ID: {transformer_id} with {len(se_ids)} unique active chargers")
            
            # Check if the transformer_id exists in tf_capacity_available_KW
            if transformer_id_str not in self.trns_kva.keys():#tf_capacity_available_KW.columns:
                print(f"Warning: Transformer ID {transformer_id} not found in tf_capacity_available_KW. Skipping this transformer.")
                continue

            # Extract the available capacity time series for the corresponding transformer from the tf_capacity_available_KW
            capacity = self.trns_kva[transformer_id_str]#tf_capacity_available_KW[transformer_id_str].to_frame(name='Capacity')
            for node in self.nodes_by_trns[transformer_id]:
                if node in baseloads.keys():
                    capacity = capacity + np.array(baseloads[node])
            
            #capacity_series_monthly.index = pd.to_datetime(capacity_series_monthly.index.tz_localize(None))
            #capacity_series_monthly.index = pd.to_datetime(capacity_series_monthly.index)
            # subtract the baseload for all subsequent nodes 
            # capacity_series_1min = capacity_series_monthly.loc[Sim_start_time:Sim_end_time].resample('1T').ffill()
            # Ensure the capacity series covers the full simulation period
            full_index = pd.date_range(start=self.start_time_sec, end=self.horizon_sec, freq='1T')
            capacity_series_1min = capacity_series_monthly.reindex(full_index, method='ffill')

            
            event_index = 0
            for CE in active_CEs:
                if self.trns_by_seid[CE.SE_id] == transformer_id:
                    park_start_time = pd.to_datetime(event['park_start_timestamp']).floor('T').time()
                    park_start_timestamp = datetime.combine(Sim_start_time.date(), park_start_time)
                    duration = round(event['park_time_seconds'] / 60)  # Total connection time in minutes

                    energy_need = event['energy_kwh']
                    start_soc = event['start_soc']
                    max_charge_rate = event['Max AC Power kW']
                    charging_time_uncontrol = event['charging_time_uncontrol']

                    # Create EV object with event_row_index
                    ev = EV(
                        transformer_id=transformer_id,
                        ev_id=CE.SE_id,
                        premise_id=CE.SE_id,
                        plug_in_time=park_start_timestamp,
                        duration=duration,
                        start_SOC=start_soc,
                        energy_need=energy_need,
                        max_charge_rate=max_charge_rate,
                        event_index=event_index,
                        charging_time_uncontrol=charging_time_uncontrol)  # Add this line
                    print(f"EV {ev.ev_id} added to the CMS list, Charge Event Index: {ev.event_index}, Transformer ID: {ev.transformer_id}, Plug-in Time: {ev.plug_in_time}, Duration: {ev.duration}, Start SOC: {ev.start_SOC}, Energy Need: {ev.energy_need}, Max Charge Rate: {ev.max_charge_rate} ")
            
                    management_system.add_ev(ev)
                    event_index = event_index+1

                management_system.simulate(Sim_start_time, Sim_end_time, time_step=timedelta(seconds=self.timestep_sec))  # simulation resolution 1 min
            
                ev_power_profiles, ev_energy_profiles = management_system.get_ev_data()
                charging_events_evaluation = management_system.get_charging_events_evaluation()

                # get into Caldera accepted formate of SE_setpoint
                for ev in management_system.station.connected_evs:
                    X = SE_setpoint()
                    X.SE_id = int(ev.ev_id)
                    ev_event_id = f"{ev.ev_id}_{ev.event_index}"
                    X.PkW = management_system.ev_power_series.at[management_system.current_time, ev_event_id]
                    X.QkVAR = 0
                    PQ_setpoints.append(X)
            

        # send to caldera
        Caldera_control_info_dict = {}
        if len(PQ_setpoints)>0:
            Caldera_control_info_dict[Caldera_message_types.set_pev_charging_PQ] = PQ_setpoints

        return Caldera_control_info_dict, DSS_state_info_dict, ev_control_setpoints 
        #return all_ev_power_profiles, all_ev_energy_profiles, all_charging_events_evaluation




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Heuristic-based EV SCM for Transformer Overloading Mitigation (SCM Core Function)
# EV energy tracking version ( EV energy requirment known, stop charging when energy requirement is met)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class EV:
    def __init__(self, transformer_id, ev_id, premise_id, plug_in_time, duration, start_SOC, energy_need,max_charge_rate=9.6,charging_time_uncontrol=None, event_index=None):
        # Initialize EV attributes
        self.transformer_id = transformer_id
        self.ev_id = ev_id  # Unique identifier for the electric vehicle
        self.event_index = event_index
        self.premise_id = premise_id
        self.plug_in_time = plug_in_time  # Time when the vehicle is plugged in for charging
        self.duration = duration # Total connection time in minutes
        self.allocated_power = 0  # Currently allocated charging power
        self.start_SOC = start_SOC  # State of Charge of the vehicle's battery
        self.energy_need = energy_need  # Total energy demanded by the vehicle in kWh
        self.energy_charged = 0  # Total energy charged at the current timestamp in kWh
        self.max_charge_rate = max_charge_rate  # in kW, set based on EVSE capability
        #self.actual_charging_time = 0  # New attribute to track actual charging time
        self.full_charge_flag = 0
        self.charging_time_uncontrol = charging_time_uncontrol
        
        self.managed_charging_time = None
    
    def is_connected(self, current_time):  
        # Check if the EV is connected at the current timestamp
        # Connect time (dwell period) can be longer than charging time
        return self.plug_in_time <= current_time <= (self.plug_in_time + timedelta(minutes=self.duration))

class ChargingStation:
    def __init__(self):
        self.connected_evs = []

    def add_ev(self, ev):
        self.connected_evs.append(ev)

    def remove_ev(self, ev):
        self.connected_evs.remove(ev)

class ChargingManagementSystem:
    def __init__(self, capacity_series, allocation_method='uncontrol'):
        self.station = ChargingStation()
        self.current_time = None
        self.capacity_series = capacity_series
        self.available_capacity_series = pd.Series(dtype=float, index=capacity_series.index)
        self.ev_power_series = pd.DataFrame(index=capacity_series.index)
        self.ev_energy_series = pd.DataFrame(index=capacity_series.index)
        self.allocation_method = allocation_method
        self.charging_events_evaluation = []  # New attribute to store evaluation results

    def update_time(self, new_time):
        self.current_time = new_time

    def add_ev(self, ev):
        self.station.add_ev(ev)
        
        # Use the row_index as the event_id
        event_id = f"{ev.ev_id}_{ev.event_index}"
        
        # Initialize EV power series with zeros
        if event_id not in self.ev_power_series.columns:
            self.ev_power_series[event_id] = 0
        
        # Initialize EV energy need series
        energy_need_id = f"{event_id} Energy Need"
        if energy_need_id not in self.ev_energy_series.columns:
            self.ev_energy_series[energy_need_id] = 0
        
        mask = (self.ev_energy_series.index >= ev.plug_in_time) & (self.ev_energy_series.index < ev.plug_in_time + timedelta(minutes=ev.duration))
        self.ev_energy_series.loc[mask, energy_need_id] = ev.energy_need
        
        # Initialize EV energy charged series with zeros
        energy_charged_id = f"{event_id} Energy Charged"
        if energy_charged_id not in self.ev_energy_series.columns:
            self.ev_energy_series[energy_charged_id] = 0
        
        #print(f"Added EV event: {event_id}, Plug-in Time: {ev.plug_in_time}, Duration: {ev.duration}, Energy Need: {ev.energy_need}")
    
    
    def get_ev_data(self):
        return self.ev_power_series, self.ev_energy_series
        

    def allocate_power_uncontrol(self):
        connected_evs = [ev for ev in self.station.connected_evs if ev.is_connected(self.current_time)]
        
        for ev in connected_evs:
            ev.allocated_power = ev.max_charge_rate


    def allocate_power_first_come_first_served(self):
        remaining_capacity = self.capacity_series[self.current_time]
        sorted_evs = sorted(self.station.connected_evs, key=lambda ev: ev.plug_in_time)
        
        for ev in sorted_evs:
            if ev.is_connected(self.current_time):
                if remaining_capacity >= 1.44:
                    allocatable_power = max(min(remaining_capacity, ev.max_charge_rate), 1.44)
                    ev.allocated_power = allocatable_power
                    remaining_capacity -= ev.allocated_power
                else:
                    ev.allocated_power = 0  # Not enough capacity to meet minimum requirement

    def allocate_power_fcfs_with_minimum(self):
        connected_evs = [ev for ev in self.station.connected_evs if ev.is_connected(self.current_time)]
        num_connected_evs = len(connected_evs)
        
        if num_connected_evs == 0:
            return
        
        total_capacity = self.capacity_series[self.current_time]
        min_power = 1.44  # Minimum power allocation (kW)
        average_power = 0.5 * total_capacity / num_connected_evs
        
        # Sort EVs by plug-in time (FCFS)
        sorted_evs = sorted(connected_evs, key=lambda ev: ev.plug_in_time)
        
        # First pass: Allocate average power if it's at least the minimum power
        remaining_capacity = total_capacity
        if average_power >= min_power:
            for ev in sorted_evs:
                ev.allocated_power = min(average_power, ev.max_charge_rate)
                remaining_capacity -= ev.allocated_power
        else:
            # If average power is less than minimum, allocate no power
            for ev in sorted_evs:
                ev.allocated_power = 0
        
        # Second pass: Distribute remaining capacity
        if remaining_capacity > 0:
            for ev in sorted_evs:
                additional_power = min(remaining_capacity, ev.max_charge_rate - ev.allocated_power)
                ev.allocated_power += additional_power
                remaining_capacity -= additional_power
                
                if remaining_capacity <= 0:
                    break
        
        # Final check: If any EV got less than minimum power, set it to 0
        for ev in connected_evs:
            if 0 < ev.allocated_power < min_power:
                ev.allocated_power = 0

    # Calculates the equal power share for all connected EVs.
    # If the equal share is above the minimum power, it allocates this share to all EVs (limited by their max charge rate).
    # If the equal share is below the minimum power, it allocates the minimum power to as many EVs as possible.
    # It then redistributes any remaining capacity among the EVs that received power.
    # This approach maintains the principle of equal sharing while being more computationally efficient. It ensures that:
    # When capacity is sufficient, all EVs receive an equal share.
    # When capacity is insufficient for all EVs to receive the minimum power, it allocates power to as many as possible.
    # Any remaining capacity is distributed equally among charging EVs.

    def allocate_power_equal_sharing(self):
        connected_evs = [ev for ev in self.station.connected_evs if ev.is_connected(self.current_time)]
        if not connected_evs:
            return

        total_capacity = self.capacity_series[self.current_time]
        num_evs = len(connected_evs)
        equal_power = total_capacity / num_evs
        min_power = 1.44  # Minimum power allocation (kW)

        if equal_power >= min_power:
            # Allocate equal power to all EVs
            for ev in connected_evs:
                ev.allocated_power = min(equal_power, ev.max_charge_rate)
        else:
            # Allocate minimum power to as many EVs as possible
            evs_to_charge = int(total_capacity / min_power)
            for i, ev in enumerate(connected_evs):
                if i < evs_to_charge:
                    ev.allocated_power = min(min_power, ev.max_charge_rate)
                else:
                    ev.allocated_power = 0

        # Redistribute any remaining capacity
        remaining_capacity = total_capacity - sum(ev.allocated_power for ev in connected_evs)
        if remaining_capacity > 0:
            charging_evs = [ev for ev in connected_evs if ev.allocated_power > 0]
            additional_power = remaining_capacity / len(charging_evs)
            for ev in charging_evs:
                ev.allocated_power = min(ev.allocated_power + additional_power, ev.max_charge_rate)

            
    def allocate_power_based_on_soc(self):
        connected_evs = [ev for ev in self.station.connected_evs if ev.is_connected(self.current_time)]
        if connected_evs:
            remaining_capacity = self.capacity_series[self.current_time]
            total_inverse_soc = sum(1 - ev.current_soc for ev in connected_evs)
        
            # First round allocation based on SOC
            for ev in connected_evs:
                proportion = (1 - ev.current_soc) / total_inverse_soc
                allocatable_power = min(proportion * remaining_capacity, ev.max_charge_rate)
                ev.allocated_power = max(allocatable_power, 1.44) if allocatable_power >= 1.44 else 0
                remaining_capacity -= ev.allocated_power
  
            # Second round allocation
            evs_under_max_power = [ev for ev in connected_evs if ev.allocated_power < ev.max_charge_rate]
            if evs_under_max_power and remaining_capacity > 0:
                total_inverse_soc_under_max = sum(1 - ev.current_soc for ev in evs_under_max_power)
                for ev in evs_under_max_power:
                    proportion = (1 - ev.current_soc) / total_inverse_soc_under_max
                    additional_power = min(proportion * remaining_capacity, ev.max_charge_rate - ev.allocated_power)
                    ev.allocated_power += additional_power
                    remaining_capacity -= additional_power


    def allocate_power_priority_factors(self):
        connected_evs = [ev for ev in self.station.connected_evs if ev.is_connected(self.current_time)]
        total_priority = sum(ev.priority_factor for ev in connected_evs)
        remaining_capacity = self.capacity_series[self.current_time]
        for ev in connected_evs:
            allocatable_power = min((ev.priority_factor / total_priority) * remaining_capacity, 9.6)
            ev.allocated_power = allocatable_power

    def evaluate_charging_event(self, ev):
        energy_satisfaction_ratio = ev.energy_charged / ev.energy_need if ev.energy_need > 0 else 1.0
        
        if self.allocation_method != 'uncontrol' and ev.managed_charging_time is not None:
            charging_time_ratio = ev.managed_charging_time / ev.charging_time_uncontrol
        else:
            charging_time_ratio = 1.0  # For uncontrolled charging or if managed_charging_time is not set
        
        return {
            'ev_id': ev.ev_id,
            'event_index': ev.event_index,
            'transformer_id': ev.transformer_id,
            'full_charge_flag': ev.full_charge_flag,
            'energy_satisfaction_ratio': energy_satisfaction_ratio,
            'charging_time_ratio': charging_time_ratio,
            'allocation_method': self.allocation_method
        } 
    
    def simulate(self, start_time, end_time, time_step=timedelta(minutes=1)):
        self.current_time = start_time
        
        while self.current_time <= end_time:
            self.update_time(self.current_time)

            if self.allocation_method == 'uncontrol':
                self.allocate_power_uncontrol()
            if self.allocation_method == 'first_come_first_served':
                self.allocate_power_first_come_first_served()
            if self.allocation_method == 'fcfs_with_minimum':
                self.allocate_power_fcfs_with_minimum()    
            elif self.allocation_method == 'equal_sharing':
                self.allocate_power_equal_sharing()
            elif self.allocation_method == 'soc_priority':
                self.allocate_power_based_on_soc()
            elif self.allocation_method == 'priority_factors':
                self.allocate_power_priority_factors()
            
        
            for ev in list(self.station.connected_evs):
                # Check if the EV is still connected or if this is its last time step
                if ev.is_connected(self.current_time) or self.current_time + time_step > ev.plug_in_time + timedelta(minutes=ev.duration):
                    ev.energy_charged += ev.allocated_power * (time_step.total_seconds() / 3600)  # energy = power * time (hours)
                    
                    event_id = f"{ev.ev_id}_{ev.event_index}"
                    energy_charged_id = f"{event_id} Energy Charged"

                    self.ev_power_series.at[self.current_time, event_id] = ev.allocated_power if ev.is_connected(self.current_time) else 0
                    self.ev_energy_series.at[self.current_time, energy_charged_id] = ev.energy_charged
                    
                    # Check if the EV is fully charged or if this is its last time step
                    if ev.energy_charged >= ev.energy_need or self.current_time + time_step > ev.plug_in_time + timedelta(minutes=ev.duration):
                        ev.full_charge_flag = 1 if ev.energy_charged >= ev.energy_need else 0
                        ev.managed_charging_time = (self.current_time - ev.plug_in_time).total_seconds() / 60
                        evaluation = self.evaluate_charging_event(ev)
                        self.charging_events_evaluation.append(evaluation)
                        
                        if ev.full_charge_flag:
                            print(f"EV {ev.ev_id}, Charge Event Index: {ev.event_index} fully charged. Time taken: {ev.managed_charging_time:.2f} minutes, allocation method: {self.allocation_method}, Evaluated: {evaluation}")
                        else:
                            print(f"EV {ev.ev_id}, Charge Event Index: {ev.event_index} did not get fully charged. Time taken: {ev.managed_charging_time:.2f} minutes, allocation method: {self.allocation_method}, Evaluated: {evaluation}")
                        
                        self.station.remove_ev(ev)

            self.current_time += time_step

    def get_charging_events_evaluation(self):
        return pd.DataFrame(self.charging_events_evaluation)
 
