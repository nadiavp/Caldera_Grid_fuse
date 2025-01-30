import numpy as np
import math
import scipy.optimize as opt
import time
import datetime, calendar
#======================================================
#      BTM control algorithm starts here
#======================================================     
class BTM_Control():
    def __init__(self, time_step_mins, ess_size, max_power_ess, min_power_ess, 
                 max_power_l2, min_power_l2, time_horizon):
        self.time_step_mins = time_step_mins   # Emin/Emax time step: 5 minutes
        self.ess_size       = ess_size        # in kWh
        self.max_power_ess  = max_power_ess   # 1C
        self.min_power_ess  = min_power_ess
        self.max_power_l2   = max_power_l2
        self.min_power_l2   = min_power_l2
        self.time_horizon   = time_horizon      # in hour
        
        self.ess_soc_min = 0.2
        self.num_of_active_evse = 0
                 
    def objective(self, x, p):
        x = x + p
        m = max(x)/np.mean(x)
        sum = 0
        for i in x:
            sum += i
    
        return m/sum + m
    
    def constraints(self, x, Emin, Emax):
        e = np.zeros(len(Emin))
        #for i in range(0,len(Emin)-1):
        for i in range(0,len(Emin)-1):
            e[i+1] = e[i] + x[i]*self.time_step_mins/60
    
        p_const1 = []
        p_const2 = []
        e_const1 = []
        e_const2 = []
        
        ### EVSE constraints >= 0
        for i in range(len(Emin)):
            p_const1 += [ x[i] - self.min_power_l2 - self.min_power_ess - 0*1e-09 ]
            p_const2 += [ self.max_power_l2*self.num_of_active_evse + self.max_power_ess - 0*1e-09 - x[i] ]
            e_const1 += [e[i] - Emin[i]]
            e_const2 += [Emax[i] - e[i]]
    
        return np.array(p_const1 + p_const2 + e_const1 + e_const2)
    
    
    def solve_optimization(self, x0, p, Emin, Emax):
        print('Solving BTM Control algorithm... ')
        myoptions={'disp':False, "maxiter": 500}
        cons = ({'type': 'ineq', 'fun': self.constraints})
        results = opt.minimize(self.objective, x0, args=(p,),
                       constraints={'type': 'ineq', 'fun': lambda x, Emin, Emax: self.constraints(x, Emin, Emax), 'args': (Emin, Emax)},
                       method='COBYLA', options = myoptions) 
                    
        return results
    
    def energy_constraints_ess(self, current_time, depart_time, energy_remaining, available_energy):
        max_power = self.max_power_ess
        min_power = self.min_power_ess
        print(f'energy_remaining: {energy_remaining} \n available_energy:{available_energy} \n max_power:{max_power}')
        hr_to_complete_charge    = (energy_remaining+available_energy)/max_power
        hr_to_complete_discharge = available_energy/max_power
        print(f'depart_time:{depart_time} hr_to_complete_charge:{hr_to_complete_charge}, current_time:{current_time}')
        min_to_start = int((math.floor(depart_time - hr_to_complete_charge*60*60) - current_time)/60)
        
    
    
        offset = list(range(self.time_step_mins, self.time_horizon*60+int(self.time_step_mins*2/3), self.time_step_mins))  # 6 hr window with 15 minute resolution
        #offset = list(range(0,time_horizon*60+int(time_step_mins*2/3),time_step_mins))
        Emax = []
        Emin = []
    
        for i in offset:
            Emax.append(min(max_power*i/60, energy_remaining))
            if i <= min_to_start:
                ##Emin.append(0.0)
                Emin.append(max(min_power*i/60, -available_energy))
            else:
                Emin.append(min(-available_energy+max_power*(i-min_to_start)/60, energy_remaining))
                
        if hr_to_complete_charge > (depart_time-current_time)/60/60:
            return Emax, Emax
        else:
            return Emin, Emax
        
    def energy_constraints_evse(self, current_time, depart_time, energy_remaining):
        max_power = self.max_power_l2
        min_power = self.min_power_l2
        
        hr_to_complete_charge = energy_remaining/max_power
        
        min_to_start = math.floor((depart_time - hr_to_complete_charge*60*60 - current_time)/60)
        if min_to_start < 0 and energy_remaining>0:
            #print(f'WARNING: not enough time to complete charge at {current_time}, energy_kWh:{energy_remaining}, evse_maxkW:{max_power}, min_to_start:{min_to_start}')
            min_to_start = 0

    
        offset = list(range(self.time_step_mins, self.time_horizon*60+1, self.time_step_mins))
        Emax = []
        Emin = []
        energy_so_far = 0.0
        last_idx = offset[0]
    
        for i in offset:
            Emax_i = min(max_power*i/60, energy_remaining)
            if i < min_to_start-self.time_step_mins:
                ##Emin.append(0.0)
                Emin.append(min(min_power*i/60, energy_remaining))
                energy_so_far = min(min_power*i/60, energy_remaining)
                last_idx = i
            else:
                #Emin.append(min(max_power*(i-min_to_start)/60, energy_remaining))
                Emin.append(min(energy_so_far + max_power*(i-last_idx)/60, 
                                energy_remaining))
            Emax.append(max(Emax_i, Emin[-1]))
        
        if hr_to_complete_charge > (depart_time-current_time)/60/60:
            return(Emax, Emax)
        else:
            return(Emin, Emax)
        
        
    def allocate_setpoint(self, setpoint, soc_ess, available_energy_ess, depart_time_evse, energy_remaining_evse):
        setpoint_ess = 0.0
        setpoint_evse = [0.0]*len(energy_remaining_evse)
        
        if -self.min_power_ess*self.time_step_mins/60 < available_energy_ess:
            setpoint_ess = self.min_power_ess
        else:
            setpoint_ess = -available_energy_ess*60/self.time_step_mins
            
        remaining_pwr = setpoint - setpoint_ess
        tmp = [y/x for x, y in zip(depart_time_evse, energy_remaining_evse)]
        priority = [x/sum(tmp) for x in tmp]
        
        setpoint_evse = [min(remaining_pwr*x, self.max_power_l2) for x in priority]
        
        index = np.flip(np.argsort(setpoint_evse))
        for k in range(len(index)-1):
            if setpoint_evse[index[k]] > self.max_power_l2:
                setpoint_evse[index[k+1]] = setpoint_evse[index[k+1]] + setpoint_evse[index[k]] - self.max_power_l2
                setpoint_evse[index[k]] = self.max_power_l2
        
        print('Sending BTM Control Setpoints... ')

        return setpoint_ess, setpoint_evse
#======================================================
#      BTM control algorithm ends here
#======================================================        
