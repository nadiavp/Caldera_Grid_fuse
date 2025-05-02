### this script scans the opendss model
### and develops a dictionary of distribution buses with btm storage
### that dictionary includes storage size, 
### max charging and discharging power, and min charging and discharing power
### all values here are in kW or kva
import pandas as pd
import numpy as np
from opendssdirect import dss
import sys
import os

def get_load_points(dss_loads=[], load_threshold=0, pv_threshold=5):
    # only add the load point if the threshold is above 350kw
    # or pv capacity is above 5 kw
    bus_df = pd.DataFrame(columns=['bus_name','load_peak', 'pv_cap'])
    for load_name in dss_loads:
        # if the load is negative, assume pv generation
        # if positive assume normal load
        dss.Loads.Name(load_name)
        kva_base = dss.Loads.kVABase()
        loadshape_name = dss.Loads.Yearly() # try .Daily() if only 24 hour sim
        if not loadshape_name=='':
            loadshape = dss.LoadShape.Name(loadshape_name)
            peak = max(dss.LoadShape.PMult())*kva_base
        else:
            peak = dss.Loads.kW()
        bus_name = dss.CktElement.BusNames()[0]
        if -peak > pv_threshold:
            # first check to see if it's been added yet
            if not (bus_name in bus_df['bus_name']):
                bus_df = pd.concat([bus_df, pd.DataFrame.from_records({'bus_name':bus_name,'load_peak':0, 'pv_cap':0}, index=[0])])
            bus_df.loc[bus_df['bus_name']==bus_name,'pv_cap'] = bus_df.loc[bus_df['bus_name']==bus_name,'pv_cap'] + peak
        elif peak > load_threshold:
            # first check to see if it's been added yet
            if not (bus_name in bus_df['bus_name']):
                bus_df = pd.concat([bus_df, pd.DataFrame.from_records({'bus_name':bus_name,'load_peak':0, 'pv_cap':0}, index=[0])])
            bus_df.loc[bus_df['bus_name']==bus_name,'load_peak'] = bus_df.loc[bus_df['bus_name']==bus_name,'load_peak'] + peak
    return bus_df


def assign_battery_sizes(bus_df):
    # small_residential
    smr_bes = {'storage_cap_kwh':5, 'storage_power_kw':5}
    # large residential
    lgr_bes = {'storage_cap_kwh':30, 'storage_power_kw':5}
    # small commercial
    smc_bes = {'storage_cap_kwh':1100, 'storage_power_kw':350}
    # large commercial
    lgc_bes = {'storage_cap_kwh':1100, 'storage_power_kw':350}
    # dataframe of storage by bus
    bess_df = {'bus_name':[], 'storage_cap_kwh':[], 'storage_power_kw':[], "storage_SOC":[], 'bes_eff':[], 'Net_load':[]}
    for _, bus in bus_df.iterrows():
        bess_df['bus_name'].append(bus['bus_name'])
        if bus['pv_cap'] > 25 or bus['load_peak'] > 1000: # in kW
            cap = lgc_bes['storage_cap_kwh']
            pwr = lgc_bes['storage_power_kw']
        elif bus['pv_cap'] > 12 or bus['load_peak'] > 350:
            cap = smc_bes['storage_cap_kwh']
            pwr = smc_bes['storage_power_kw']
        elif bus['pv_cap'] > 7 or bus['load_peak'] > 15:
            cap = lgr_bes['storage_cap_kwh']
            pwr = smc_bes['storage_power_kw']
        else:
            cap = smr_bes['storage_cap_kwh']
            pwr = smr_bes['storage_power_kw']
        bess_df['storage_cap_kwh'].append(cap)
        bess_df['storage_power_kw'].append(pwr)
        bess_df["storage_SOC"].append(.5)
        bess_df['bes_eff'].append(0.98)
        bess_df['Net_load'].append(0)
    ## convert to dataframe
    #bess_df = pd.DataFrame.from_dict(bess_df)
    return bess_df

def get_btms_siting(opendss_file):
    print(f'stiting storage on dss model {opendss_file}')
    dss_model = dss.run_command(f'Redirect {opendss_file}')
    dss.Solution.Solve()
    dss_loads = dss.Loads.AllNames()
    bus_df = get_load_points(dss_loads)
    bess_df = assign_battery_sizes(bus_df)
    return bess_df

def add_btms_to_opendss_model(opendss_main_file, bess_df):
    print(f'creating dss file of storage from {opendss_file}')
    bess_file_name = 'BTM_BESS_and_PV.dss'
    bess_dss_file = opendss_main_file.replace('Master.dss',bess_file_name)
    # load the opendss model so that you can get the bus phases and voltages later
    if not os.path.exists(opendss_main_file):
        print(f'opendss file: {opendss_main_file} does not exist, continuing to next feeder')
        return
    
    dss.Command(f'Redirect {opendss_main_file}')
    # initialize list of new storage and solar
    storage_pv_str_list = []
    bus_phases = 1
    # add all the storage and PV as a lines in a new file
    bess_keys = bess_df.columns.values
    for index, bess_i in bess_df.iterrows():
        # get the sizing info
        # if the info hasn't been calculated, calculate it
        if not 'node_name' in bess_keys:
            #print(f'index: {index} bess_i: {bess_i}')
            bus_name = bess_i[0]
            peak_kw = bess_i['peak_kW']
            pv_size = np.round(peak_kw*0.25) # a quarter of the peak power
            power_size = pv_size/2 # half the kw of the pv array
            energy_size = power_size*4 # 4 hours of max power
        else:
            bus_name = bess_i['node_name']
            power_size = bess_i['batt_kW']
            energy_size = bess_i['batt_kWh']
            pv_size = bess_i['pv_kW']
        #bus_phases = len(bus_name.split('.'))-1
        # get the voltage from the bus name
        dss.Circuit.SetActiveBus(bus_name)
        bus_voltage = dss.Bus.kVBase()
        bus_phases = len(dss.Bus.Voltages())
        # if you have multiple phases on the bus, pick the first one to attach the pv and storage
        if bus_phases>1:
            bus_phases = 1
            bus_name = bus_name.split('.')[0] + '.' + bus_name.split('.')[1]
        # add storages
        new_storage_str = f"New Storage.{bus_name} phases={bus_phases} Bus1={bus_name} kV={bus_voltage}  kW={power_size}  kWrated={power_size}  kWhrated={energy_size} dispmode=external"
        storage_pv_str_list.append(new_storage_str)
        # add PV
        if pv_size>0:
            new_pv_str = f"New PVSystem.{bus_name} bus1={bus_name} phases={bus_phases} kV={bus_voltage} kVA={pv_size} Pmpp={pv_size}" # setting Pmpp same as kVA means no temperature degradation
            storage_pv_str_list.append(new_pv_str)
    # save the storage details to a file
    with open(bess_dss_file, 'w') as dss_file:
        for row in storage_pv_str_list:
            dss_file.write(row+"\n")    
    # make the main redirect to the BTM_BESS_and_PV.dss
    new_main_lines = []
    with open(opendss_main_file) as main_file:
        main_lines = main_file.readlines()
        for line in main_lines:
            if line.startswith('Set Voltagebases'):
                new_main_lines.append(f'Redirect {bess_file_name} \n')
            new_main_lines.append(line)
    with open(opendss_main_file, 'w') as main_file:
        for line in new_main_lines:
            main_file.write(line)

if __name__ == "__main__":
    # first load the loadshape file
    if len(sys.argv)>1:
        opendss_file = sys.argv[1]
        dss_model = dss.run_command(f'Redirect {opendss_file}')
        dss.Solution.Solve()
        dss_loads = dss.Loads.AllNames()
    else:
        dss_loads = []
    #bus_df = get_load_points(dss_loads)
    #bess_df = assign_battery_sizes(bus_df)
    #print(bess_df)


    # read all sheets in xcel file and create the btms sizing dss file
    # for pandas version >= 0.21.0
    file_name = 'mhdv_depot_pv_and_storage_sizing.xlsx'
    sheet_to_df_map = pd.ExcelFile(file_name)#, sheet_name=None)

    ## for pandas version < 0.21.0
    #sheet_to_df_map = pd.read_excel(file_name, sheetname=None)
    for feeder in sheet_to_df_map.sheet_names:
        bess_df = sheet_to_df_map.parse(feeder, index=False)
        opendss_file = f'opendss/{feeder}/Master.dss'
        add_btms_to_opendss_model(opendss_file, bess_df)
