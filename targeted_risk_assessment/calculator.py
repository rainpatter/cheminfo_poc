
import pandas as pd

df = pd.read_csv('../data/ECETOC-TRAworker-version3.2-final.XLSX - TRAlookup.csv')

# # define parameters

# substance_name = 'ethanol'
# cas_number = '64-17-5'
# #g/mol
# mol_weight = 46.069

# #DNEL OR OEL
# # mg/m3
# long_term_inhalation = 950
# # mg/kg bw/day
# long_term_dermal = 206

# #mg/m3
# short_term_inhalation = 1900

# # ug/cm2
# local_dermal = 100000

# vap_pressure_at_operating_temp = 7832.4225

example_user_inputs = {
    'substance_name': 'ethanol',
    'cas_number': '64-17-5',
    'mol_weight': 46.069,
    'long_term_inhalation': 950,
    'long_term_dermal': 206,
    'short_term_inhalation': 1900,
    'local_dermal': 10000,
    'vap_pressure_at_operating_temp': 7832.4225,
    'proc': 'PROC7',
    'ind_prof': 'ind',
    'phys_state': 'liquid',
    'fugacity': 'low',
    'ventilation': 'indoors - good ventilation',
    'duration': '15min-1hr',
    'concentration': '>25%',
    'lev': 'no',
    'rpe_mask': 'no RPE',
    'ppe_gloves': 'PPE95%',
    'lev_dermal': 'yes'
}


def calc_fugacity_band(dict):
    vap = dict['vap_pressure_at_operating_temp']
    if vap < 0.01:
        fugacity_band = 'very low'
    elif 500 > vap >= 0.01:
        fugacity_band = 'low'
    elif 10000 >= vap >= 500:
        fugacity_band = 'medium'
    else:
        fugacity_band = 'high'
    return fugacity_band


def calculate_ventilation_reduction_factor(dict):
    vent = dict['ventilation']
    if vent == 'outdoors' or 'indoors - good ventilation':
        vrf = 0.7
    elif vent == 'indoors - enhanced ventilation':
        vrf = 0.3
    else:
        vrf = 1
    return vrf


def calculate_duration_reduction_factor_inhalation(dict):
    dur = dict['duration']
    if dur == '<15min':
        drfi = 0.1
    elif dur == '15min-1hr':
        drfi = 0.2
    elif dur == '1-4hr':
        drfi = 0.6
    else:
        drfi = 1
    return drfi

# =IF(AND(D12="solid",OR(W12="medium",W12="high")),1,IF(AND(D12="liquid",OR(W12="very low",W12="low")),1,Y12))


def calculate_duration_reduction_factor_dermal(dict, calc_vrf):
    phys = dict['phys_state']
    fug = dict['fugacity']
    if (phys == 'solid') and (fug == 'medium' or 'high'):
        drfd = 1
    elif (phys == 'liquid') and (fug == 'very low' or 'low'):
        drfd = 1
    else:
        drfd = calc_vrf
    return drfd


def calculate_concentration_reduction_factor(dict):
    conc = dict['concentration']
    if conc == '<1%':
        crf = 0.1
    elif conc == '1-5%':
        crf = 0.2
    elif conc == '5-25%':
        crf = 0.2
    else:
        crf = 1
    return crf


def calculate_rpe_reduction_factor(dict):
    rpe = dict['rpe_mask']
    if rpe == 'RPE90%':
        rrf = 0.1
    elif rpe == 'RPE95%':
        rrf = 0.05
    else:
        rrf = 1
    return rrf

# =IF(K12="PPE80%",0.2, IF(K12="PPE90%",0.1, IF(AND(C12="prof",K12="PPE95%"),0.1, IF(K12="PPE95%",0.05,1))))


def calculate_ppe_reduction_factor(dict):
    ip = dict['ind_prof']
    ppe = dict['ppe_gloves']
    if ppe == 'PPE80%':
        prf = 0.2
    elif ppe == 'PPE90%':
        prf = 0.1
    elif (ip == 'prof') and (ppe == 'PPE95%'):
        prf = 0.1
    elif ppe == 'PPE95%':
        prf = 0.05
    else:
        prf = 1
    return prf

# IF(AND(AND(D12="liquid",W12="very low"),NOT(OR(B12="PROC7",B12="PROC11",B12="PROC17",B12="PROC18",AND(B12="PROC10",I12="no"),AND(B12="PROC19",I12="no")))),1,4)


def calcule_multiplier_short_term(dict):
    phys = dict['phys_state']
    fug = dict['fugacity']
    proc = dict['proc']
    lev = dict['lev']
    if (phys == 'liquid' and fug == 'very low') and not (proc in ['PROC7', 'PROC11', 'PROC17', 'PROC18'] and (proc == 'PROC10' and lev == 'no') and (proc == 'PROC19' and lev == 'no')):
        mst = 1
    else:
        mst = 4
    return mst

def generate_lookup_descriptor(dict):
    proc = dict['proc']
    phys = dict['phys_state']
    lev = dict['lev']
    fug = dict['fugacity']
    ip = dict['ind_prof']
    string = proc+phys+lev+fug+ip
    return string

print(calc_fugacity_band(example_user_inputs))
print(calculate_ventilation_reduction_factor(example_user_inputs))
vrf = calculate_duration_reduction_factor_inhalation(example_user_inputs)
print(vrf)
print(calculate_duration_reduction_factor_dermal(example_user_inputs, vrf))
print(calculate_concentration_reduction_factor(example_user_inputs))
print(calculate_rpe_reduction_factor(example_user_inputs))
print(calculate_ppe_reduction_factor(example_user_inputs))
print(calcule_multiplier_short_term(example_user_inputs))
print(generate_lookup_descriptor(example_user_inputs))

print(df)
print(df.columns.values)