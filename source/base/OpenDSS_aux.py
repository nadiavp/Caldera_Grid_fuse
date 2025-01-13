
import opendssdirect as dss
import os, math
import numpy as np
import pandas as pd

from global_aux import OpenDSS_message_types, input_datasets, non_pev_feeder_load


class open_dss:

    # NOTE: The 'open_dss' class does two things: Use OpenDSS, and do logging.
    # if the boolean 'use_opendss' is false, then this class will just do the logging.
    # TODO: Separate the logging into a separate federate.
    def __init__(self, io_dir, use_opendss):

        #if use_opendss == True:
        print(f'use_opendss is {use_opendss}')
        self.helper = open_dss_helper(io_dir)
        #else:
        #    self.helper = logger_helper(io_dir)


    def get_input_dataset_enum_list(self):
        return self.helper.get_request_list()


    def load_input_datasets(self, datasets_dict):
        self.helper.load_input_datasets(datasets_dict)


    def initialize(self):       
        return self.helper.initialize()
    

    def process_control_messages(self, simulation_unix_time, message_dict, der_busnames=[]):        
        return self.helper.process_control_messages(simulation_unix_time, message_dict, der_busnames)

    
    def set_caldera_pev_charging_loads(self, node_pevPQ):
        self.helper.set_caldera_pev_charging_loads(node_pevPQ)


    def set_der_charge_controlb(self, controlb_setpoint):
        self.helper.set_der_charge_controlb(controlb_setpoint)    
    
    def get_pu_node_voltages_for_caldera(self):
        return self.helper.get_pu_node_voltages_for_caldera()

    def get_node_load_profile_for_controlb(self, t_now, t_horizon, t_step, der_busnames=[]):
        return self.helper.get_node_load_profile_for_controlb(t_now, t_horizon, t_step, der_busnames)

    #def get_netload_node_for_controlb(self):
    #    return self.helper.get_netload_node_for_controlb()

    def get_der_soc_for_controlb(self):
        return self.helper.get_der_soc_for_controlb()    
    
    def solve(self, simulation_unix_time):
        self.helper.solve(simulation_unix_time)        
    

    def log_data(self, simulation_unix_time):
        self.helper.log_data(simulation_unix_time)


    def post_simulation(self):
        self.helper.post_simulation()


class open_dss_helper:

    def __init__(self, io_dir):
        self.io_dir = io_dir
        self.dss_file_name = 'Master.dss' #'ieee34.dss'
        self.feeder_name = self.io_dir.feeder_name 
        self.scenario = self.io_dir.scenario_name 

    def get_request_list(self):
        return [input_datasets.baseLD_data_obj, input_datasets.all_caldera_node_names, input_datasets.HPSE_caldera_node_names]

    def load_input_datasets(self, datasets_dict):
        # datasets_dict is a dictionary with input_datasets as keys.
        self.datasets_dict = datasets_dict

    def initialize(self):

        baseLD_data_obj = self.datasets_dict[input_datasets.baseLD_data_obj]
        all_caldera_node_names = self.datasets_dict[input_datasets.all_caldera_node_names]
        HPSE_caldera_node_names = self.datasets_dict[input_datasets.HPSE_caldera_node_names]
        
        self.dss_core = open_dss_core(self.io_dir, self.dss_file_name, self.feeder_name, baseLD_data_obj)
        self.dss_Caldera = open_dss_Caldera(self.io_dir, all_caldera_node_names, HPSE_caldera_node_names)
        self.dss_external_control = open_dss_external_control()
        
        #-------------------------------------------
        #         Load and Check .dss file
        #-------------------------------------------
        is_successful = self.dss_core.load_dss_file()
    
        if is_successful:
            is_successful = self.dss_Caldera.check_caldera_node_names()
        
        #----------------------------------------------------
        #  Create logger object:   
        #----------------------------------------------------
        self.dss_logger = None
        if(is_successful):
            self.dss_logger = open_dss_logger_A(self.io_dir, self.scenario, self.feeder_name, all_caldera_node_names, HPSE_caldera_node_names)

        return is_successful

    def process_control_messages(self, simulation_unix_time, message_dict, der_busnames):        
        return self.dss_external_control.process_control_messages(simulation_unix_time, message_dict, der_busnames)
    
    def set_caldera_pev_charging_loads(self, node_pevPQ):
        self.node_pevPQ = node_pevPQ
        self.dss_Caldera.set_caldera_pev_charging_loads(node_pevPQ)

    def set_der_charge_controlb(self, controlb_setpoint):
        self.dss_external_control.set_der_charge_controlb(controlb_setpoint)

    def get_pu_node_voltages_for_caldera(self):    
        return self.dss_Caldera.get_pu_node_voltages_for_caldera()

    def get_der_soc_for_controlb(self):
        return self.dss_external_control.get_der_soc_for_controlb()

    def get_node_load_profile_for_controlb(self, t_now, t_horizon, t_step, der_busnames):
        return self.dss_external_control.get_node_load_profile_for_controlb(t_now, t_horizon, t_step, der_busnames)

    #def get_netload_node_for_controlb(self):
    #    return self.dss_external_control.get_netload_node_for_controlb()

    def solve(self, simulation_unix_time):
        self.dss_core.solve(simulation_unix_time)

    def log_data(self, simulation_unix_time):
        self.dss_logger.log_data(simulation_unix_time, self.node_pevPQ)
    
    def post_simulation(self):
        pass

class logger_helper:

    def __init__(self, io_dir):
        self.io_dir = io_dir

    def get_request_list(self):
        return [input_datasets.baseLD_data_obj, input_datasets.all_caldera_node_names]

    def load_input_datasets(self, datasets_dict):
        # datasets_dict is a dictionary with input_datasets as keys.
        self.datasets_dict = datasets_dict

    def initialize(self):

        self.baseLD_data_obj = self.datasets_dict[input_datasets.baseLD_data_obj]
        self.all_caldera_node_names = self.datasets_dict[input_datasets.all_caldera_node_names]
                
        self.logger_obj = logger(self.io_dir, self.baseLD_data_obj, self.all_caldera_node_names)

        is_successful = True
        return is_successful
    
    def process_control_messages(self, simulation_unix_time, message_dict): 

        return_dict = {}
        
        for (msg_enum, parameters) in message_dict.items():
            if msg_enum == OpenDSS_message_types.get_all_node_voltages:
                return_dict[msg_enum] = self.get_pu_node_voltages_for_caldera()
            else:
                raise ValueError('Invalid message in caldera_ICM_aux::process_message.')
        
        # The return value (return_dict) must be a dictionary with OpenDSS_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict

    
    def set_caldera_pev_charging_loads(self, node_pevPQ):
        self.node_pevPQ = node_pevPQ

    def get_pu_node_voltages_for_caldera(self):
        return_dict = {}
        for node_name in self.all_caldera_node_names:
            return_dict[node_name] = 1.0
        
        return return_dict

    def solve(self, simulation_unix_time):
        self.logger_obj.compute_total_load_profiles(self.node_pevPQ, simulation_unix_time)
        
    def log_data(self, simulation_unix_time):
        return None

    def post_simulation(self):
        self.logger_obj.write_data_to_disk()


class open_dss_external_control:

    def __get_all_node_voltages(self):
        return_dict = {}
        
        all_V = dss.Circuit.AllBusMagPu()
        all_node_names = dss.Circuit.AllNodeNames()
        
        for i in range(len(all_V)):            
            return_dict[all_node_names[i]] = all_V[i]
        
        return return_dict

    def get_der_soc_for_controlb(self):
        storage_soc = {'Net_load':[], 'storage_SOC':[], 'storage_cap_kwh':[], 'storage_power_kw':[], 'names': [], 'bus_name': []}
        storage_soc['names'] = dss.Storages.AllNames()
        #print(f'storage names: {dss.Storages.AllNames()}')
        for storage_name in dss.Storages.AllNames():
            # make the stoage the active element
            dss.Storages.Name(storage_name)
            storage_soc['storage_power_kw'].append(dss.CktElement.Powers()[0]) # this is the output real power
            storage_soc['storage_SOC'].append(dss.Storages.puSOC())
            storage_soc['storage_cap_kwh'].append(dss.Storages.kVARated())
            # make the bus the active element
            storage_soc['bus_name'].append(dss.CktElement.BusNames()[0])
            dss.Buses.Name(dss.CktElement.BusNames()[0])
            storage_soc['Net_load'].append(dss.CktElment.TotalPowers()[0]) # this is real power netload
        return storage_soc

    def set_der_charge_controlb(self, controlb_setpoint=[]):
        for storage in dss.Storages.AllNames():
            dss.Storages.Name(storage)
            dss.Storages.kW(controlb_setpoint[storage])

    def get_node_load_profile_for_controlb(self, t_now, t_horizon, t_step, der_busnames):
        # this gets the daily load profile only at nodes with storage
        # it returns the sum of only the fixed loads and fixed pv
        # the storage and evse are not added here
        # this does not include reactive loads
        # t_now is the time according to the co-simulation in minutes
        # t_horizon is the horizon of the optimization in minutes
        # t_step is the timestep of the co-simulation (and therefore also optimization) in minutes
        # der_bunames is the list of bus names for any storage that was added but not 
        # put in the opendss model. This can be taken from the der_data['bus_name'] list
        # used by the nrel_control_btms_ld_l2 controller
        n_steps = int(np.ceil((t_horizon-t_now)/t_step))
        #t_now = int(np.floor(t_now))
        t_horizon = int(np.floor(t_horizon))
        #t_step = int(round(t_step))
        time_steps_desired = [ts*t_step + t_now for ts in range(n_steps)]

        #print(f'getting node load profiles for controlb:')
        #print(f'storages: {dss.Storages.AllNames()} and {der_busnames} \n ')
        netload = {}
        # first add the storage that is in the dss model
        for storage in dss.Storages.AllNames():
            # set active element to this storage
            dss.Storages.Name(storage)
            # get bus name
            bus_name = dss.CktElement.BusNames()[0]
            netload[bus_name] = np.zeros(n_steps)
        # then add the storage that's not in the model, but was added by the siting script
        for busname in der_busnames:
            netload[busname] = np.zeros(n_steps)
        # now get the load profiles for those buses
        i_load = dss.Loads.First()
        while i_load>0:
            # determine if there are loads attached to that
            loadbusname = dss.CktElement.BusNames()[0]
            if loadbusname in netload.keys():
                # if there is a loadshape then use it, check for daily or yearly profile
                # then if no loadshape, just use the base kva for the full profile
                loadshape_name = dss.Loads.Daily()
                if len(loadshape_name)<1:
                    loadshape_name = dss.Loads.Yearly()
                base_kva_0 = dss.Loads.kVABase()
                if len(loadshape_name)<1:
                    netload[loadbusname] = netload[loadbusname] + base_kva_0
                    print(f'loadshape not available, using basekva {base_kva_0}')
                else:
                    dss.LoadShape.Name(loadshape_name)
                    load_profile = base_kva_0*dss.LoadShape.PMult()
                    # sample only the ones for this time
                    loadshape_tstep = int(dss.LoadShape.MinInterval())
                    load_profile = np.interp(time_steps_desired, range(0, loadshape_tstep*len(load_profile),loadshape_tstep), load_profile)
                    # store the sample
                    netload[loadbusname] = netload[loadbusname] + load_profile
                    print(f'loadshape at {loadbusname} with profile {load_profile}')
            i_load = dss.Loads.Next()
        # do the same for the pv systems
        i_pv = dss.PVsystems.First()
        while i_pv>0:
            pvbusname = dss.CktElement.BusNames()[0]
            if pvbusname in netload.keys():
                loadshape_name = dss.PVsystems.daily()
                kva_rated = dss.PVsystems.kVARated()
                dss.LoadShape.Name(loadshape_name)
                pv_profile = kva_rated*dss.LoadShape.PMult()
                # sample only the ones for this time
                loadshape_tstep = int(dss.LoadShape.HrInterval())
                pv_profile = np.interp(time_steps_desired, range(0, loadshape_tstep*len(pv_profile),loadshape_tstep), pv_profile)
                # storage the sample
                netload[bus_name] = netload[bus_name] + load_profile
            i_pv = dss.PVsystems.Next()
        # if there arent any pv or storage, do it for all load buses
        #print(f'netloads: {netload} \n getting all load bus loadshapes')
        if len(netload.keys())== 0:
            i_load = dss.Loads.First()
            while i_load>0:
                loadshape_name = dss.Loads.Daily()
                base_kva_0 = dss.Loads.kVABase()
                loadbusname = dss.CktElement.BusNames()[0]
                if len(loadshape_name)<1:
                    loadshape_name = dss.Loads.Yearly()
                if len(loadshape_name)>1:
                    dss.LoadShape.Name(loadshape_name)
                    load_profile = base_kva_0*dss.LoadShape.PMult()
                    # sample only the ones for this time
                    loadshape_tstep = int(dss.LoadShape.MinInterval())
                    #print(f'load: {dss.Loads.Name()}, loadshape_tstep:{loadshape_tstep}, load_profile:{load_profile}')
                    load_profile = np.interp(time_steps_desired, range(0, loadshape_tstep*len(load_profile),loadshape_tstep), load_profile)
                else: # if there is no loadshape, assume a constant load
                    load_profile = np.array([base_kva_0]*n_steps)
                # store the sample
                if loadbusname in netload.keys():
                    netload[loadbusname] = netload[loadbusname] + load_profile
                else:
                    netload[loadbusname] = load_profile
                #print(f'loadshape at {loadbusname} with profile {load_profile}')
                i_load = dss.Loads.Next()
        # convert numpy arrays to lists
        for busname in netload.keys():
            netload[busname] = netload[busname].tolist()
        #print(f'returning netload: {netload}')
        return netload


    #def get_netload_node_controlb(self):
    #    # this gets the netloads only at buses with storage
    #    # this only returns the real load at the current time
    #    netload = {}
    #    buses_w_storage = []
    #    for storage in dss.Storages.AllNames():
    #        # set active element to this storage
    #        dss.Storages.Name(storage)
    #        # get connected bus
    #        bus_name = dss.CktElement.BusNames()[0]
    #        # set the ative element to the bus
    #        dss.Buses.Name(bus_name)
    #        # get the real power at he node
    #        netload[bus_name] = dss.CktElement.TotalPowers()[0]
    #
    #    return netload

    
    def process_control_messages(self, simulation_unix_time, msg_dict, der_busnames):
        # msg_dict is a dictionary with OpenDSS_message_types as keys.
        
        return_dict = {}
        
        for (msg_enum, parameters) in msg_dict.items():
            if msg_enum == OpenDSS_message_types.get_all_node_voltages:
                return_dict[msg_enum] = self.__get_all_node_voltages()
            elif msg_enum == OpenDSS_message_types.get_all_DER:
                return_dict[msg_enum] = self.get_der_soc_for_controlb()
            elif msg_enum == OpenDSS_message_types.get_basenetloads:
                return_dict[msg_enum] = self.get_node_load_profile_for_controlb(t_now=simulation_unix_time/3600,t_horizon=simulation_unix_time/3600+12,t_step=0.25, der_busnames=der_busnames)
            else:
                raise ValueError('Invalid message in caldera_ICM_aux::process_message.')
        
        # The return value (return_dict) must be a dictionary with OpenDSS_message_types as keys.
        # If there is nothing to return, return an empty dictionary.
        return return_dict



class open_dss_core:
    

    def __init__(self, io_dir, dss_file_name, feedername, baseLD_data_obj):
        self.io_dir = io_dir
        self.output_path = self.io_dir.outputs_dir
        self.dss_file_name = dss_file_name
        self.feeder_load = non_pev_feeder_load(baseLD_data_obj)
        self.feeder_name = feedername
        self.ref_feeder_kW = 1
    
    
    def __configuring_dss(self):
        #dss.run_command('set controlMode=TIME') 
        dss.Solution.ControlMode(2) # 0->STATIC, 1->EVENT, 2->TIME
        dss.Solution.SolveSnap()
        (ref_feeder_kW, ref_feeder_kVAR) = dss.Circuit.TotalPower()
        self.ref_feeder_kW = abs(ref_feeder_kW)
    
    
    def load_dss_file(self):
        is_successful = True
        
        dss_filepath = os.path.join( self.io_dir.base_dir, 'opendss', self.feeder_name, self.dss_file_name )
        opendss_input_file_exists = os.path.isfile(dss_filepath)
        
        #-----------------------
        
        if not opendss_input_file_exists or not dss.Basic.Start(0):
            is_successful = False
        
            if not opendss_input_file_exists:
                print('OpenDSS input file does not exist.  Path to file: {}'.format(dss_filepath))
            else:
                print('OpenDSS not started!')
        
        #-----------------------
        
        if is_successful:        
            dss.Basic.ClearAll()
            dss.Basic.DataPath(self.output_path)
            opendss_load_status = dss.run_command('Compile ['+ dss_filepath + ']')
            
            if opendss_load_status != '':
                is_successful = False
                print('Unable to Compile OpenDSS input file. Error message: {}'.format(opendss_load_status))
        
        #-----------------------
        
        if is_successful:
            self.__configuring_dss()
        
        return is_successful
    
    
    def solve(self, simulation_unix_time):
        # Scaling the non-pev feeder load
        feeder_load_akW = self.feeder_load.get_non_pev_feeder_load_akW(simulation_unix_time)
        load_multiplier =  feeder_load_akW / self.ref_feeder_kW
        dss.Solution.LoadMult(load_multiplier)
        
        # Solving Powerflow
        hours = math.floor(simulation_unix_time/3600)
        seconds = simulation_unix_time - hours*3600
        dss.Solution.Hour(hours)
        dss.Solution.Seconds(seconds)
        dss.Solution.SolveSnap()

        converged = dss.Solution.Converged()
        if not converged:
            print('OpenDSS simulation did NOT converge at simulation time: {} hours.'.format(simulation_unix_time/3600))


class open_dss_Caldera:
    
    def __init__(self, io_dir, all_caldera_node_names, HPSE_caldera_node_names):
        self.io_dir = io_dir
        self.all_caldera_node_names = all_caldera_node_names
        self.HPSE_caldera_node_names = HPSE_caldera_node_names
        
    
    def check_caldera_node_names(self):
        errors = []
        
        all_node_names = dss.Circuit.AllNodeNames()
        all_bus_names = dss.Circuit.AllBusNames()
        all_load_names = dss.Loads.AllNames()
        
        for node_name in self.all_caldera_node_names:
            if node_name in self.HPSE_caldera_node_names:
                if node_name not in all_bus_names:
                    errors.append('Error: The Caldera node_id: ({}) for fast charging does not correspond to a bus in Open-DSS.'.format(node_name))
                
                load_name = 'pev3p_' + node_name
                if load_name not in all_load_names:
                    errors.append('Error: The fast charging loads on Caldera node_id: ({}) does not have a corresponing Open-DSS load object.  The Open-DSS load object should be named ({}).'.format(node_name, load_name))
            else:
                if node_name not in all_node_names:
                    errors.append('Error: The Caldera node_id: ({}) for L1 and L2 charging does not correspond to a node in Open-DSS.'.format(node_name))
                
                load_name = 'pev1p_' + node_name
                if load_name not in all_load_names:
                    errors.append('Error: The L1 & L2 loads on Caldera node_id: ({}) does not have a corresponing Open-DSS load object.  The Open-DSS load object should be named ({}).'.format(node_name, load_name))
        
        #-----------------------
        
        is_successful = True
        f_invalid_nodes = open( os.path.join( self.io_dir.inputs_dir, 'error_invalid_node_in_SE_file.txt' ), 'w')
        
        if len(errors) > 0:
            is_successful = False
            print('Invalid OpenDSS nodes in Caldera SE file.')
            
            for msg in errors:
                f_invalid_nodes.write(msg + '\n')
        
        f_invalid_nodes.close()
        
        #-----------------------
        
        return is_successful
    
    
    def get_pu_node_voltages_for_caldera(self):
        return_dict = {}
        
        #---------------
        # Single Phase
        #---------------
        all_V = dss.Circuit.AllBusMagPu()
        all_node_names = dss.Circuit.AllNodeNames()
        
        for i in range(len(all_V)):
            if all_node_names[i] in self.all_caldera_node_names:
                return_dict[all_node_names[i]] = all_V[i]
        
        #---------------
        # Avg 3 Phase
        #---------------
        for bus_name in self.HPSE_caldera_node_names:    
            dss.Circuit.SetActiveBus(bus_name)
            V_complex_pu = dss.Bus.PuVoltage()
            num_nodes = int(round(len(V_complex_pu)/2))
            pu_V = 0
            
            for i in range(num_nodes):
                pu_V += (V_complex_pu[2*i] **2 + V_complex_pu[2*i+1] **2) ** (1/2) # complex absolute value
                #pu_V += dss.CmathLib.cabs(V_complex_pu[2*i], V_complex_pu[2*i+1])
            
            return_dict[bus_name] = pu_V / num_nodes
        
        return return_dict


    def set_caldera_pev_charging_loads(self, node_pevPQ): 
        for (node_id, (P_kW, Q_kVAR)) in node_pevPQ.items(): 
            if node_id in self.HPSE_caldera_node_names:
                dss.Loads.Name('pev3p_' + node_id)
                dss.Loads.kW(P_kW)
                dss.Loads.kvar(Q_kVAR)
            else:
                dss.Loads.Name('pev1p_' + node_id)
                dss.Loads.kW(P_kW)
                dss.Loads.kvar(Q_kVAR)



class open_dss_logger_A:

    def __init__(self, io_dir, scenario, feeder_name, all_caldera_node_names, HPSE_caldera_node_names):
        
        node_voltages_to_log = [] #'810.2', '822.1', '826.2', '856.2', '864.1', '848.1', '848.2', '848.3', '840.1', '840.2', '840.3', '838.2', '890.1', '890.2', '890.3']
        #node_pev_charging_to_log = ['810.2', '826.2', '856.2', '838.2']
        node_pev_charging_to_log = list(all_caldera_node_names)
        
        #------------------------------
        
        self.io_dir = io_dir
        self.all_caldera_node_names = all_caldera_node_names
        self.HPSE_caldera_node_names = HPSE_caldera_node_names
        
        openDSS_node_names = set()
        X = dss.Circuit.AllNodeNames()
        for x in X:
            openDSS_node_names.add(x)
        
        #output_path = self.io_dir.outputs_dir
        output_folder_path = self.io_dir.outputs_dir
        sc_path = os.path.join(output_folder_path, scenario)
        if not os.path.exists(sc_path):
            os.mkdir(sc_path)
        output_path = os.path.join(sc_path, feeder_name)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        
        #--------------------------------------
        #           feeder_PQ.csv
        #--------------------------------------
        self.f_feeder = open(output_path + '/feeder_PQ.csv', 'w')
        self.f_feeder.write('simulation_time_hrs, feeder_kW, pev_kW, feeder_kVAR, pev_kVAR' + '\n')
        
        #--------------------------------------
        #     Selected_Node_Voltages.csv
        #--------------------------------------
        self.node_voltages_to_log = []
        
        header = 'simulation_time_hrs'
        for node_id in node_voltages_to_log:
            if node_id in openDSS_node_names:
                self.node_voltages_to_log.append(node_id)
                header += ', _' + node_id 
        
        self.f_V = open(output_path + '/Selected_Node_Voltages.csv', 'w')
        self.f_V.write(header + '\n')
        
        #--------------------------------------
        #       node_pev_charging_to_log
        #--------------------------------------
        self.f_node_pev_charging = {}
        for x in node_pev_charging_to_log:
            if x in self.all_caldera_node_names:
                self.f_node_pev_charging[x] = open(output_path + '/' + x + '.csv', 'w')
                self.f_node_pev_charging[x].write('simulation_time_hrs, pu_Vrms, pev_kVAR, pev_kW' + '\n')
    
    
    def __del__(self):
        self.f_feeder.close()
        self.f_V.close()
        
        for (node_id, f_node) in self.f_node_pev_charging.items():
            f_node.close()
    
    
    def log_data(self, simulation_unix_time, node_pevPQ):
        simulation_time_hrs = simulation_unix_time/3600
    
        (feeder_kW, feeder_kVAR) = dss.Circuit.TotalPower()
        feeder_kW = -feeder_kW
        feeder_kVAR = -feeder_kVAR
        
        pev_kW = 0
        pev_kVAR = 0        
        for (node_id, (P_kW, Q_kVAR)) in node_pevPQ.items():
            pev_kW += P_kW
            pev_kVAR += Q_kVAR
        
        node_puV = {}
        all_V = dss.Circuit.AllBusMagPu()
        all_node_names = dss.Circuit.AllNodeNames()
        
        for i in range(len(all_V)):
            node_puV[all_node_names[i]] = all_V[i]
        
        #--------------------------------------
        #           feeder_PQ.csv
        #--------------------------------------
        tmp_str = '{}, {}, {}, {}, {}'.format(simulation_time_hrs, feeder_kW, pev_kW, feeder_kVAR, pev_kVAR)
        self.f_feeder.write(tmp_str + '\n')
        
        #--------------------------------------
        #     Selected_Node_Voltages.csv
        #--------------------------------------
        node_voltage_str = '{}'.format(simulation_time_hrs)
        for node_id in self.node_voltages_to_log:
            node_voltage_str += ', {}'.format(node_puV[node_id])
        
        self.f_V.write(node_voltage_str + '\n')
    
        #--------------------------------------
        #       node_pev_charging_to_log
        #--------------------------------------
        for (node_id, f_node) in self.f_node_pev_charging.items():   
            (pevP_kW, pevQ_kVAR) = node_pevPQ[node_id]
            
            if node_id in node_puV:
                node_V = node_puV[node_id]
            else:
                node_V = 0.0
                dss.Circuit.SetActiveBus(node_id)
                X = dss.Bus.Nodes()
                for x in X:
                    node_name = node_id + "." + str(x)
                    node_V += node_puV[node_name]

                node_V = node_V/len(X)
            
            tmp_str = '{}, {}, {}, {}'.format(simulation_time_hrs, node_V, pevQ_kVAR, pevP_kW)
            f_node.write(tmp_str + '\n')

class logger:

    def __init__(self, io_dir, baseLD_data_obj, all_caldera_node_names):

        self.io_dir = io_dir
        self.all_caldera_node_names = all_caldera_node_names
        self.baseLD_data_obj = baseLD_data_obj
        #print("all_caldera_node_names : {}".format(all_caldera_node_names))
        
        self.real_power_profiles = {}
        self.reactive_power_profiles = {}
        self.real_power_profiles["simulation_time_hrs"] = []
        self.real_power_profiles["base_load_kW"] = []
        self.real_power_profiles["total_demand_kW"] = []
        
        self.reactive_power_profiles["simulation_time_hrs"] = []
        self.reactive_power_profiles["base_load_kW"] = []
        self.reactive_power_profiles["total_demand_kW"] = []

        for node_name in all_caldera_node_names:
            self.real_power_profiles[node_name] = []
            self.reactive_power_profiles[node_name] = []


    def compute_total_load_profiles(self, node_pevPQ, simulation_unix_time):
  
        simulation_time_hrs = simulation_unix_time/3600.0

        index = math.floor((simulation_unix_time - self.baseLD_data_obj.data_start_unix_time) / self.baseLD_data_obj.data_timestep_sec)

        if (index < 0) or (index >= len(self.baseLD_data_obj.actual_load_akW)):
            print("Error : base_LD index computed not in data range")
            exit()

        base_LD_kW = self.baseLD_data_obj.actual_load_akW[index]
        self.real_power_profiles["simulation_time_hrs"].append(simulation_time_hrs)
        self.real_power_profiles["base_load_kW"].append(base_LD_kW)

        self.reactive_power_profiles["simulation_time_hrs"].append(simulation_time_hrs)
        self.reactive_power_profiles["base_load_kW"].append(base_LD_kW)

        total_P_kW = 0.0
        total_Q_kVAR = 0.0
        for (node_name, (P_kW, Q_kVAR)) in node_pevPQ.items():
            self.real_power_profiles[node_name].append(P_kW)
            self.reactive_power_profiles[node_name].append(Q_kVAR)
            total_P_kW += P_kW
            total_Q_kVAR += Q_kVAR
    
        self.real_power_profiles["total_demand_kW"].append(total_P_kW)
        self.reactive_power_profiles["total_demand_kW"].append(total_Q_kVAR)

    def write_data_to_disk(self):
        output_dir = self.io_dir.outputs_dir

        df = pd.DataFrame(self.real_power_profiles)
        df.to_csv( os.path.join( output_dir, "real_power_profiles.csv" ), index=False)

        df = pd.DataFrame(self.reactive_power_profiles)
        df.to_csv( os.path.join( output_dir, "reactive_power_profiles.csv" ), index=False)

    def get_pu_node_voltages_for_caldera(self):

        return_dict = {}
        for node_name in self.all_caldera_node_names:
            return_dict[node_name] = 1.0
        
        return return_dict