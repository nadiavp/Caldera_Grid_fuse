
from math import floor
import numpy as np
import pandas as pd

from scipy.optimize import minimize
import sys
import copy
from __main__ import *
#======================================================
#      Open-DSS Feedback Control algorithm starts here
#======================================================     

class DataInputs():

    #def __init__(self,filename,N,comm_delay,num_charging_stations,max_ev,Node_pu,Vpu_Lower,Vpu_Higher,Voltage_multiplier,charge_event_id,arrival_time,departure_time,remaining_charge_energy_ackWh):
    def __init__(self,comm_delay,num_charging_stations,Node_pu,Vpu_Lower,Vpu_Higher,Voltage_multiplier,CE_by_SE_groups,charge_event_id_Parse,SE_id_Parse,arrival_times_Parse,departure_times_Parse,charging_times_Parse,energy_used_Parse,group_id_Parse,SetPoint_Dict):
        
        #self.filename = filename                           # filename with all of the events
        #df_lookup = pd.read_csv(filename)
        self.CE_by_SE_groups = CE_by_SE_groups
        #self.N = N                                          # number of charging events
        self.comm_delay = comm_delay                        # communication delay
        self.t = np.arange(0, 5 + 1, 1)                     # time (minutes)
        #self.load_data() 
        #print(self.load_data())                                   # custom data loader
        #self.compute_power()                                # compute total power consumed without control
        #self.compute_peak()                                 # compute peak power with no control
    
        # for hierarchical control
        self.num_charging_stations = 1  # number of charging stations
        #self.max_ev = max_ev                                # maximum number of vehicles at one charging station
        #self.arrival_times = arrival_time
        #self.end_times = departure_time
        # which charging events go to which charging stations
        self.Node_pu = Node_pu
        self.Vpu_Lower = Vpu_Lower
        self.Vpu_Higher = Vpu_Higher
        self.Voltage_multiplier = Voltage_multiplier
        self.SetPoint_Dict = SetPoint_Dict
        self.SE_id = SE_id_Parse
        
    #def load_data(self):
               
        self.group_id = group_id_Parse
        self.charge_event_id = charge_event_id_Parse
        N = len(self.charge_event_id)
        self.N = N                                          # number of charging events
        self.max_ev = N                                # maximum number of vehicles at one charging station
        self.max_ev = len(self.charge_event_id)                                # maximum number of vehicles at one charging station
        self.SE_id = SE_id_Parse

        self.arrival_times = arrival_times_Parse * 60 #converting time to minutes
        #print('arrival times', self.arrival_times)
        self.departure_times = departure_times_Parse * 60 #converting time to minutes
        #print('departure_times', self.departure_times)

        self.charging_times = charging_times_Parse

              
        self.energy_used = energy_used_Parse          
        
        
    #def compute_power(self):
        #print('self.N = ', self.N)
        power = np.zeros((self.N, 5))
        for i in range(self.N):
            #idx_start = np.min(np.where(self.t >= self.arrival_times[i]))
            #idx_end = np.min(np.where(self.t >= self.end_times[i]))
            idx_start = self.arrival_times[i]
            idx_end = self.departure_times[i]
            hours = (self.departure_times[i] - self.arrival_times[i]) / 60 #Retained division by 60
            #print('hour = ', hours)
            #print('self.energy_used[i]=' , self.energy_used[i])
            Nt = idx_end - idx_start
            #power[i, idx_start:idx_end] = int((self.energy_used[i] / hours) * np.ones(Nt))
            power[i, 0:5] = (self.energy_used[i] / hours)
            
        self.power = power
        #print('power =', self.power)
        #print('size of power =', len(self.power))
    
    #def compute_peak(self):
    
        self.peak = np.max(np.sum(self.power,axis=0))
        #print('peak =', self.peak)
 
 