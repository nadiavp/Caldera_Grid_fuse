import pandas as pd 


# update the next line with the CE file you need to update
ce_file_name = 'CE_longdwell_Sep_Hanover_01359.csv'#'CE_Sep_Shellbank_22700.csv'

# names must be of length 7 and start with ext
# if they include reactive power control they must be of length 8 and end with q
# names must be in format ext0001q where ext_string[3:7] can be converted to a positive integer
market_control = 'ext0001'
emissions_control = 'ext0002'
depot_control = 'ext0003'
transformer_control = 'ext0004'
voltwatt_control = 'ext0005q'

# update the next 4 lines to match the control types you want
light_duty_control = market_control
transit_control = emissions_control
school_control = depot_control
local_control = depot_control

transit_ev_types = ['BEB_1000_kWh', 'BEB_1500_kWh', 'BEB_2000_kWh', 'Proterra ZX5 + 40-feet', 'Proterra ZX5 Max 40-feet']
school_ev_types = ['Type A_1', 'Type C_1', 'Type D_1', 'Type A_2', 'Type C_2', 'Type D_2', 'Type A_3', 'Type C_3','Type C_4', 'Type C_5']
local_ev_types = ['StraightTruck_3','StraightTruck_4','StraightTruck_5','StraightTruck_6', 'Cutaway_2', 'Cutaway_3', 'CargoVan_2', 'CargoVan_3',
    'StepVan_3', 'StepVan_4', 'StepVan_5', 'StepVan_6']
# assume everything not in a list above is a ldv

ce_df = pd.read_csv(ce_file_name)

for i, charge_event in ce_df.iterrows():
    control_code = 'ext0001q'
    if charge_event['pev_type'] in transit_ev_types:
        control_code = transit_control
    elif charge_event['pev_type'] in school_ev_types:
        control_code = school_control
    elif charge_event['pev_type'] in local_ev_types:
        control_code = local_control
    else:
        control_code = light_duty_control
    ce_df.at[i,'Ext_strategy'] = control_code
    ce_df.at[i,'VS_strategy'] = 'NA'
    ce_df.at[i,'ES_strategy'] = 'NA'
    # for versions before 0.21.0 use set_value
    # ce_df.set_value(i, 'Ext_strategy', control_code)

# save the updated csv
ce_df.to_csv(ce_file_name, index=False)

