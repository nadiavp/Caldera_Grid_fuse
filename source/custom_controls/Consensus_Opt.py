

from math import floor
import numpy as np
import pandas as pd
import time

from scipy.optimize import minimize
import sys
import copy
from __main__ import *



class Hierarchical():
        
        def __init__(self,data,max_charge,next_control_timestep_start_unix_time):
        
            # problem setup
            self.data = data
            self.max_charge = max_charge
        
            # construct network of charging stations
            self.generate_network()
        
            # algorithmic constants
            self.maxiter = 5
            self.rho = 10
            self.lbda = 0.1
            self.res = 10000.0
        
            # initialize optimal start and end times
            self.opt_start = self.data.arrival_times
            self.opt_stop = self.data.departure_times
            self.comms = 0
            current_time = next_control_timestep_start_unix_time/3600
            self.next_control_timestep_start_unix_time = next_control_timestep_start_unix_time
        
        def generate_network(self):       
            # TODO: generate network based on location or based on TNC (currently all charging stations can talk to each other)
            self.A = np.zeros((self.data.num_charging_stations,self.data.num_charging_stations))
        
            for i in range(self.data.num_charging_stations):
                for j in range(self.data.num_charging_stations):
                    if i != j:
                        self.A[i,j] = 1.0
        
        def peak_charging_total(self, x, z, u, idx, power_chargers):
        
            # number of design variables
            N_station = int(len(x)/2)
        
            # extract power used
            power = np.zeros((N_station, len(self.data.t)))
        
            for i in range(N_station):
                #idx_start = np.min(np.where(self.data.t >= x[i]))
                idx_start = np.min(x[i])
                # impose bounds
                #if x[i + N_station] > np.max(self.data.t):
                    #x[i + N_station] = np.max(self.data.t)
        
                idx_end = np.min(x[i + N_station])
                hours = (x[i + N_station] - x[i]) /60 #Retained division by 60
                Nt = idx_end - idx_start
                if Nt < 0:
                    power[i, idx_start:idx_start + 1] = 1000000000
                else:
                    power[i, idx_start:idx_end] = (self.data.energy_used[idx[i]] / hours) * np.ones(Nt)
        
            # make sure the peak power output does not exceed x_peak
        
            # objective function
            output = np.max(np.sum(power, axis=0) + power_chargers)
        
            # compare the peak with the peak from its neighbors
            x_peak = np.max(np.sum(power, axis=0) + power_chargers)
        
            # consensus piece
            for j in range(len(z)):
                output = output + (self.rho / 2) * (x_peak - z[j] + u[j]) ** 2
        
            return output
        
        def initial_conditions(self,idx):
        
            # number of events at each charging station
            N_station = len(idx)
            x0 = np.zeros(2 * N_station)
            for i, j in enumerate(idx):
                x0[i] = self.opt_start[j]
            for i, j in enumerate(idx):
                x0[i + N_station] = self.opt_stop[j]
        
            return x0
        
        def generate_constraints(self,idx):
        
            # number of events at each charging station
            N_station = len(idx)
        
            # bounds on arrival and end times
            bnds = []
            for j in idx:
                bnds.append((self.data.arrival_times[j], self.data.departure_times[j]))
            
            # arrival time always has to be smaller than departure time
            # must be possible to charge fully in the time between arrival and departure times
            cons = []
            for i in range(N_station):
                cons.append({'type': 'ineq', 'fun': lambda x: x[i + N_station] - x[i] - 1 })
        
            return bnds, cons

        def x_update(self, z, u, opt_power):
        
            # solve the optimization for each charging station
            for k in range(self.data.num_charging_stations):
        
                # find events that occur at charging station k
                # TODO: decide which charging station each event happens at in real-time
                idx = np.copy(self.data.charge_event_id)
                N_stations = len(idx)  # number of charging events at each station
        
                # determine optimal start and stop times for each event
                x0 = self.initial_conditions(idx)
        
                # setup bounds and constraints for the problem
                bnds, cons = self.generate_constraints(idx)
        
                # pass each charging station the sum of all the other charging stations
                power_chargers = np.sum(opt_power, axis=0) - np.sum(opt_power[idx, :], axis=0)
                residual = minimize(self.peak_charging_total, x0,
                                    args=(z[k, :], u[k, :], idx, power_chargers), method='SLSQP', bounds=bnds, constraints=cons)
        
                # parse optimized powers
                x = residual.x
        
                # extract power used
                for i in range(N_stations):
                    #idx_start = np.min(np.where(self.data.t >= x[i]))
                    #idx_end = np.min(np.where(self.data.t >= x[i + N_stations]))
                    idx_start = np.min(x[i])
                    idx_end = np.min(x[i + N_stations])
                    hours = (x[i + N_stations] - x[i]) / 60 #Retained division by 60
                    Nt = idx_end - idx_start
                    opt_power[idx[i], idx_start:idx_end] = (self.data.energy_used[idx[i]] / hours) * np.ones(Nt)
                    self.opt_start[idx[i]] = x[i]
                    self.opt_stop[idx[i]] = x[i + N_stations]
        
           
            # extract the peak of each charging station
            x_out = np.zeros(self.data.num_charging_stations)
            for k in range(self.data.num_charging_stations):
                #idx = np.where(self.data.event_nodes == k)[0]
                idx = np.copy(self.data.charge_event_id)
                x_out[k] = np.max(np.sum(opt_power[idx, :], axis=0) + power_chargers[k,:])
        
            return x_out, opt_power
        
        def z_update(self,x_peak,z,u):
        
            # update z which holds a copy of the communications
            # TODO: embed the network structure into the updates of z
            for i in range(self.data.num_charging_stations):
                for j in range(self.data.num_charging_stations):
                    theta_init = 1 - (
                                self.lbda / (self.rho * np.linalg.norm((x_peak[i] + u[i, j]) - (x_peak[j] + u[j, i]))))
                    theta = np.max([theta_init, 0.5])
                    z[i, j] = theta * (x_peak[i] + u[i, j]) + (1 - theta) * (x_peak[j] + u[j, i])
                    z[j, i] = (1 - theta) * (x_peak[i] + u[i, j]) + theta * (x_peak[j] + u[j, i])
                    self.comms = self.comms + 1
        
            return z
        
        def u_update(self,x_peak,z,u):
        
            # update u
            for i in range(self.data.num_charging_stations):
                for j in range(self.data.num_charging_stations):
                    u[i, j] = u[i, j] + (x_peak[i] - z[i, j])
        
            return u
         
        
        def peak_charging_real_time(self,x,z,u,idx,current_time,energy_left):
        
            # print(x)
            # minimize the peak charging at each charging station
            x = np.nan_to_num(x)
            #expectedDeparture = (energy_left / x) + self.data.arrival_times[idx]
            expectedDeparture = (energy_left / x) + self.data.arrival_times

            output = np.sum(x)  # adding in a time component
            # add in time component
            for j in range(len(energy_left)-1):
                
                if (energy_left[j] > 0):
                    output = output + 0.0001 * (expectedDeparture[j] - self.data.arrival_times[[j]])
            x_peak = np.sum(x)
        
            # consensus piece
            for j in range(len(z)):
                output = output + (self.rho / 2) * (x_peak - z[j] + u[j]) ** 2
        
            return output
        
        def x_update_real_time(self, x_peak, z, u, current_time, energy_left):
        
            # x_peak is the current measurement of the peak load - agree on the peak load over charging footprint
        
            # this is for one time step - get the optimal charging values for each location
            opt_power = np.zeros(self.data.N)  # optimal charging power for each vehicle
        
            for k in range(self.data.num_charging_stations):
                # find which events correspond to which charging stations
                #idx = np.where(self.data.event_nodes == k)[0]
                idx = np.copy(self.data.charge_event_id)
                
                #print('idx =', idx)
                #print('size of idx', len(idx) )
                bnds = []
                x0 = np.zeros(len(idx))
                #print('x0 = ', x0)
                # 1. determine the bounds of each of the charging stations
                for i in range(len(idx)):
                    #print('i =', i)
                    #print('energy_left =', energy_left[i] )
                    if (energy_left[i] > 0):
                        time_left = (self.data.departure_times[i] - current_time)
                        energy = energy_left[i]
                        minCharge = np.min([np.max([(energy / time_left), 0]), self.max_charge])
                        bnds.append((minCharge, self.max_charge))
                        x0[i] = minCharge
                    else:
                        bnds.append((0, 0))
                        x0[i] = 0.0
        
                # 2. optimize power at each of the charging stations using
                #residual = minimize(self.peak_charging_real_time, x0,args=(z[k, :], u[k, :], idx, current_time, energy_left[idx]), method='SLSQP', bounds=bnds, options={'ftol': 0.01})
                
                residual = minimize(self.peak_charging_real_time, x0,args=(z[k, :], u[k, :], idx, current_time, energy_left), method='SLSQP', bounds=bnds, options={'ftol': 0.01})            
                
                # record the peak charge at the charging station
                optCharge = residual.x
                #print('OptCharge = ', optCharge)
                #print('idx = ', idx)

                x_peak[k] = np.sum(optCharge)
                            
                # energy_left[idx] = energy_left[idx] - optCharge
                opt_power = optCharge

                #opt_power[idx] = optCharge
                #print('opt_power', opt_power[idx-1])
            return x_peak, opt_power
        
        def solve_real_time(self):
        
            # consensus on the peak load - solve in real-time
        
            print('SCM Consensus: Solving Consensus algorithm... ')
        
            # consensus initialization
            x_peak = np.zeros(self.data.num_charging_stations)
            z = np.zeros((self.data.num_charging_stations, self.data.num_charging_stations))
            u = np.zeros((self.data.num_charging_stations, self.data.num_charging_stations))
        
            # initialize x with power from each charging station with no control
            for i in range(self.data.num_charging_stations):
        
                # find which events correspond to which charging stations
                #idx = np.copy(self.data.charge_event_id)
                x_peak[i] = np.max(np.sum(self.data.power, axis=0))
                z[i, i] = np.max(np.sum(self.data.power, axis=0))
        
        
            # estimate of the peak of the charging stations
            peak_est = 0.0
        
            # energy needed by each agent
            #energy_left = copy.deepcopy(self.data.energy_used) * 60
            energy_left = self.data.energy_used 
            
            opt_power = np.zeros((self.data.N, len(self.data.t)))
        
            # check for connections in adjacency matrix
            num_connections = np.sum(self.A)
            comms = 0
            num_comms = 0
        
#            t1 = time.time()
            for n in range(len(self.data.t)):
        
        
                # reinitialize x_peak
                x_peak = np.zeros(self.data.num_charging_stations)
        
                for k in range(0, self.maxiter):
        
                    comms = comms + num_connections * self.data.comm_delay
                    num_comms = num_comms + num_connections
        
                    # Step 1: x-update / cell update
                    #current_time = self.data.t[n]
                    current_time = self.next_control_timestep_start_unix_time/3600
        
                    # update the estimate of the peak
                    if np.sum(x_peak) > peak_est:
                        peak_est = np.sum(x_peak)
        
                    x_peak, opt_charge = self.x_update_real_time(x_peak, z, u, current_time, energy_left)
        
                    # Step 2: z-update
                    for i in range(self.data.num_charging_stations):
                        for j in range(self.data.num_charging_stations):
                            theta_init = 1 - (self.lbda / (self.rho * np.linalg.norm((x_peak[i] + u[i, j]) - (x_peak[j] + u[j, i]))))
                            theta = np.max([theta_init, 0.5])
        
                            z[i, j] = theta * (x_peak[i] + u[i, j]) + (1 - theta) * (x_peak[j] + u[j, i])
                            z[j, i] = (1 - theta) * (x_peak[i] + u[i, j]) + theta * (x_peak[j] + u[j, i])
        
        
                    # Step 3: u-update
                    for i in range(self.data.num_charging_stations):
                        for j in range(self.data.num_charging_stations):
                            u[i, j] = u[i, j] + (x_peak[i] - z[i, j])
        
        
                # record the optimal charging rates and the energy left to fully charge the vehicle
                opt_power[:, n] = opt_charge
                if  self.data.Node_pu < self.data.Vpu_Lower:
                    opt_power[:, n] = self.data.Voltage_multiplier*opt_power[:, n]
                    opt_charge = self.data.Voltage_multiplier*opt_charge
                
                if  self.data.Node_pu > self.data.Vpu_Higher:
                    opt_power[:, n] = 1.25*opt_power[:, n]
                    opt_charge = 1.25*opt_charge
                    
                    
                energy_left = energy_left - opt_charge
                #print('opt_power = ', opt_power)
                #print('OptCharge = ', opt_charge)
                
#            t2 = time.time()
        
#            self.time_to_execute_TS = (t2 - t1) / self.data.num_charging_stations
            self.peak_TS = np.max(np.sum(opt_power, axis=0))
            self.opt_power_TS = opt_power
#            self.comms_TS = copy.deepcopy(num_comms)
            Set_power = opt_power.mean(axis=1)

            print('SCM Consensus: Sending Consensus Control Setpoints... ')

            SetPoint_Dict = dict(zip(self.data.SE_id,Set_power))
            
            
        ''' CORE Consensus control code ends here *** DO NOT MODIFY *** '''
