### this script scans the opendss model
### and develops a dictionary of distribution buses with btm storage
### that dictionary includes storage size, 
### max charging and discharging power, and min charging and discharing power
import pandas as pd
from opendssdirect import dss
import sys

def get_load_points(dss_loads=[], load_threshold=350, pv_threshold=5):
    # only add the load point if the threshold is above 350kw
    # or pv capacity is above 5 kw
    bus_df = pd.DataFrame({'bus_name':[],'load_peak':[], 'pv_cap':[]})
    for bus_load in dss_loads:
        # if the load is negative, assume pv generation
        # if positive assume normal load
        kva_base = bus_load.kVABase()
        loadshape_name = bus_load.Yearly() # try .Daily() if only 24 hour sim
        loadshape = dss.Loadshape.Name(loadshape_name)
        peak = max(loadshape.PMult())*kva_base
        bus_name = bus_load.bus()[0]
        if -peak > pv_threshold*1000:
            # first check to see if it's been added yet
            if not (bus_name in bus_df['bus_name']):
                bus_df._append({'bus_name':bus_name,'load_peak':0, 'pv_cap':0})
            bus_df[bus_df['bus_name']==bus_name]['pv_cap'] += peak
        elif peak > load_threshold*1000:
            # first check to see if it's been added yet
            if not (bus_name in bus_df['bus_name']):
                bus_df._append({'bus_name':bus_name,'load_peak':0, 'pv_cap':0})
            bus_df[bus_df['bus_name']==bus_name]['load_peak'] += peak
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
    bess_df = {'bus_name':[], 'storage_cap_kwh':[], 'storage_power_kw':[]}
    for _, bus in bus_df.iterrows():
        breakpoint()
        bess_df['bus_name'].append(bus)
        if bus['pv_cap'] > 25000 or bus['load_peak'] > 1e6:
            cap = lgc_bes['storage_cap_kwh']
            pwr = lgc_bes['storage_power_kw']
        elif bus['pv_cap'] > 12000 or bus['load_peak'] > 350000:
            cap = smc_bes['storage_cap_kwh']
            pwr = smc_bes['storage_power_kw']
        elif bus['pv_cap'] > 7000 or bus['load_peak'] > 15000:
            cap = lgr_bes['storage_cap_kwh']
            pwr = smc_bes['storage_power_kw']
        else:
            cap = smr_bes['storage_cap_kwh']
            pwr = smr_bes['storage_power_kw']
        bess_df['storage_cap_kwh'].append(cap)
        bess_df['storage_power_kw'].append(pwr)
    # convert to dataframe
    bess_df = pd.DataFrame.from_dict(bess_df)
    return bess_df

if __name__ == "__main__":
    # first load the loadshape file
    if len(sys.argv)>1:
        opendss_file = sys.argv[1]
        dss_model = dss.Command(f'Redirect {opendss_file}')
        dss.Solve()
        dss_loads = dss.Loads
    else:
        dss_loads = []
    bus_df = get_load_points(dss_loads)
    bess_df = assign_battery_sizes(bus_df)
    print(bess_df)