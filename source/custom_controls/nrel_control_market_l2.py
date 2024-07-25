import pandas as pd
import helics as h
import json
import numpy as np
import logging
logger = logging.getLogger(__name__)

from Caldera_globals import SE_setpoint
from global_aux import Caldera_message_types, OpenDSS_message_types, input_datasets, container_class
from control_templates import typeB_control

from market_control_block import LPMarketController
"""
This class is to hold the controller 
It gets inputs on grid status from the grid sim 
It is initialized with info on electricity pricing and forecasting
It takes inputs on plug-in time, departure time, energy needs from the mobility analysis module
It then does some SCM control and sends those EV control setpoints to the EV Sim
"""

class market_control(typeB_control):
    def __init__(self, base_dir, simulation_time_constraints, input_se_csv='inputs/SE_.csv',
        name='market_control', helics_config_path='', timestep_sec=60*5, feeder_name='ieee_34', horizon_sec=24*60*60):
        super().__init__(base_dir, simulation_time_constraints)
        # add important params here
        self.name = name
        self.feeder_name = feeder_name
        #logging.basicConfig(filename=f'{name}.log', encoding='utf-8', level=logging.DEBUG)
        self.pricing = []
        self.baseload_forecast = []
        self.charge_events = []
        self.control_setpoints = {'evse0':0, 'evse1':0, 'evse2':0}
        self.voltages = []
        # these are for if you want time-step based sim
        self.timestep_sec = timestep_sec
        self.horizon_sec = horizon_sec
        self.time = -1
        self.helics_config_path = helics_config_path
        self.fed = None
        self.publications = []
        self.subscriptions = []

    def get_input_dataset_enum_list(self):
        return [input_datasets.SE_group_configuration, input_datasets.SE_group_charge_event_data, input_datasets.SEid_to_SE_type, input_datasets.external_strategies]


    def load_input_datasets(self, datasets_dict):
        # datasets_dict is a dictionary with input_datasets as keys.
        self.datasets_dict = datasets_dict
    
    
    def terminate_this_federate(self):
        #print(self.datasets_dict[input_datasets.external_strategies])
        if "ext_market_l2" in self.datasets_dict[input_datasets.external_strategies]:
            print(f'running with market_l2 federate')
            return False
        elif "ext0001q" in self.datasets_dict[input_datasets.external_strategies]:
            print(f'running with market_l2 federate')
            return False

        return True

    def log_data(self):
        pass

    
    def get_messages_to_request_state_info_from_Caldera(self, next_control_timestep_start_unix_time):
        return_dict = {}
        return_dict[Caldera_message_types.get_active_charge_events_by_SE_groups] = [2] #Grid teams to update this
        #return_dict[Caldera_message_types.get_active_charge_events_by_extCS] = ['ext0003', 'ext_market_l2']

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

    def initialize(self):
        # this function loads the pricing, forecast, 
        # loads plug-in time, departure time, and energy needs
        dss_file_name = 'opendss/'+self.feeder_name+'/Master.dss'     
        self.market = LPMarketController(dss_file_name=dss_file_name, feeder_name=self.feeder_name, helics_config_path=self.helics_config_path, horizon_sec=self.horizon_sec, timestep_sec=self.timestep_sec )
        prices_export = [[0.8]]*int(24*3600/self.timestep_sec)
        self.market.setup_market_controller(prices_export=prices_export, demand_charge=0)

        #if self.helics_config_path == '':
        #    #fedinfo = h.helicsCreateFederateInfo()
        #    #h.helicsFederateInfoSetCoreName(fedinfo, self.name)
        #    #h.helicsFederateInfoSetCoreInitString(fedinfo, "--federates=1")
        #    #self.fed = h.helicsCreateValueFederate(self.name, fedinfo)
        #    #self.publications.append(h.helicsFederateRegisterPublication(self.fed, 'control_setpoints', h.HelicsDataType.STRING)) # this is published as a json string of a dictionary
        #    #self.subscriptions.append(h.helicsFederateRegisterSubscription(self.fed, f'{self.feeder_name}/voltages', ""))
        #    #self.subscriptions.append(h.helicsFederateRegisterSubscription(self.fed, f'{self.feeder_name}/currents', ""))
        #    #self.subscriptions.append(h.helicsFederateRegisterSubscription(self.fed, f'{self.feeder_name}/netloads', ""))
        #else:
        #    print(f'creating market controller federate from {self.helics_config_path}')
        #    self.fed = h.helicsCreateValueFederateFromConfig(self.helics_config_path)
        #h.helicsFederateSetTimeProperty(self.fed, h.helics_property_time_delta, self.timestep_sec)

        return

    def solve(self, federate_time, Caldera_state_info_dict, DSS_state_info_dict):
        # this function solves the control parameters
        ev_control_setpoints = {}

        # first get the updated grid status
        #voltages = json.loads(h.helicsInputGetString(self.subscriptions[0]))
        #currents = json.loads(h.helicsInputGetString(self.subscriptions[1]))
        #DSS_state_info_dict = get_messages_to_request_state_info_from_OpenDSS(self, federate_time)

        ev_control_setpoints = self.market.solve(DSS_state_info_dict)

        return ev_control_setpoints



    def advance_time(self, updated_time):
        while self.time < updated_time:
            self.time = h.helicsFederateRequestTime(self.fed, updated_time)
            

    def output_control_setpoints(self):
        # this function either records setpoints or 
        # for co-simulation sends them as a helics publication
        #pd.DataFrame(self.control_setpoints).to_csv(f'{self.name}_setpoints.csv')
        h.helicsFederateEnterExecutingMode(self.fed)
        h.helicsPublicationPublishString(self.publications[0], json.dumps(self.control_setpoints))
        return 

