### this script scans the opendss model
### and develops a dictionary of distribution buses with btm storage
### that dictionary includes storage size, 
### max charging and discharging power, and min charging and discharing power
### all values here are in kW or kva
import pandas as pd
from opendssdirect import dss
import sys

def get_load_points(dss_loads=[], load_threshold=25, pv_threshold=5):
    # only add the load point if the threshold is above 350kw
    # or pv capacity is above 5 kw
    bus_df = pd.DataFrame({'bus_name':[],'load_peak':[], 'pv_cap':[]})
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
                bus_df._append({'bus_name':bus_name,'load_peak':0, 'pv_cap':0})
            bus_df.loc[bus_df['bus_name']==bus_name,'pv_cap'] = bus_df.loc[bus_df['bus_name']==bus_name,'pv_cap'] + peak
        elif peak > load_threshold:
            # first check to see if it's been added yet
            if not (bus_name in bus_df['bus_name']):
                bus_df = bus_df._append({'bus_name':bus_name,'load_peak':0, 'pv_cap':0}, ignore_index=True)
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
        bess_df['bus_name'].append(bus)
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
    # convert to dataframe
    bess_df = pd.DataFrame.from_dict(bess_df)
    return bess_df

def get_btms_siting(opendss_file):
    print(f'stiting storage on dss model {opendss_file}')
    dss_model = dss.run_command(f'Redirect {opendss_file}')
    dss.Solution.Solve()
    dss_loads = dss.Loads.AllNames()
    bus_df = get_load_points(dss_loads)
    bess_df = assign_battery_sizes(bus_df)
    return bess_df

if __name__ == "__main__":
    # first load the loadshape file
    if len(sys.argv)>1:
        opendss_file = sys.argv[1]
        dss_model = dss.run_command(f'Redirect {opendss_file}')
        dss.Solution.Solve()
        dss_loads = dss.Loads.AllNames()
    else:
        dss_loads = []
    bus_df = get_load_points(dss_loads)
    bess_df = assign_battery_sizes(bus_df)
    print(bess_df)