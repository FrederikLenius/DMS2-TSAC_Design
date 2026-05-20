# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 21:40:42 2024

@author: joj
"""
# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse
import pandas as pd
import os

plt.rcParams.update({'font.size':22, 'font.family': 'serif', 'pgf.texsystem': 'pdflatex'})



# =============================================================================
# Functions
# =============================================================================

def ReadPatchResultFile(filename):
    print('Msg: Data file being loaded: %s' %filename)
    Data = np.genfromtxt(filename, dtype=None,
                         delimiter=",", skip_header=1, usecols = (0,1,2,3,4,5,6,7,8),
                         names=['Fx', 'Fy', 'Fz', 'BF_Max' , 'FoS_sleeve_bdl', 'N_fr_max', 'FoS_fr', 'FoS_BearingStress_w_patch', 'sol_err_var'])    
    return Data


def ReadAnsysForcesFile(filename, sheet_name):
    print('Msg: Data file being loaded: %s' % filename)
    fx_values = []
    fy_values = []
    fz_values = []
    hole_values = []
    direction_values = []

    if not os.path.exists(filename):
        print('Msg: File not found: %s' % filename)
        return (
            np.array(fx_values, dtype=float),
            np.array(fy_values, dtype=float),
            np.array(fz_values, dtype=float),
            hole_values,
            direction_values,
        )

    try:
        df = pd.read_excel(filename, sheet_name=sheet_name, header=None)

        current_hole = None

        for _, row in df.iterrows():
            label = row.iloc[0]

            if pd.isna(label):
                continue

            label = str(label).strip()

            if 'hole' in label.lower():
                current_hole = label
                continue

            if current_hole is not None:
                numeric_values = pd.to_numeric(row.iloc[1:4], errors='coerce')
                if numeric_values.isna().any():
                    continue

                fx = abs(float(numeric_values.iloc[0]))
                fy = abs(float(numeric_values.iloc[1]))
                fz = abs(float(numeric_values.iloc[2]))

                hole_values.append(current_hole)
                direction_values.append(label)
                fx_values.append(fx)
                fy_values.append(fy)
                fz_values.append(fz)

        print('Msg: Loaded Ansys force cases: %d rows' % len(fx_values))
    except Exception as e:
        print('Msg: Failed to read %s: %s' % (filename, e))

    return (
        np.array(fx_values, dtype=float),
        np.array(fy_values, dtype=float),
        np.array(fz_values, dtype=float),
        hole_values,
        direction_values,
    )


ResultData = ReadPatchResultFile('Result_PatchFailureData.csv')

# =============================================================================
# Identify Failuire Envelope
# =============================================================================
FoS_low_bound = 0.99
FoS_up_bound = 1.01

ind_FoS_fr = np.where((ResultData['FoS_fr']>=FoS_low_bound) & (ResultData['FoS_fr']<=FoS_up_bound))
print('Found FoS_fr: ', ResultData['FoS_fr'][ind_FoS_fr])

ind_FoS_sl = np.where((ResultData['FoS_sleeve_bdl']>=FoS_low_bound) & (ResultData['FoS_sleeve_bdl']<=FoS_up_bound))
print('Found FoS_sl: ', ResultData['FoS_sleeve_bdl'][ind_FoS_sl])

ind_FoS_bear = np.where((ResultData['FoS_BearingStress_w_patch']>=FoS_low_bound) & (ResultData['FoS_BearingStress_w_patch']<=FoS_up_bound))
print('Found FoS_bear: ', ResultData['FoS_BearingStress_w_patch'][ind_FoS_bear])





print('')
# --- Read Ansys Forces.xlsx ---
ansys_fx, ansys_fy, ansys_fz, ansys_holes, ansys_directions = ReadAnsysForcesFile(
    'Ansys Forces.xlsx',
    'Ark1'
)
# =============================================================================
# --- Find smallest safety factor only among the points that are plotted ---
# =============================================================================
ansys_fx_valid = ansys_fx
ansys_fz_valid = ansys_fz
ansys_holes_valid = np.array(ansys_holes, dtype=object)
ansys_directions_valid = np.array(ansys_directions, dtype=object)

plotted_indices = np.unique(np.concatenate((ind_FoS_sl[0], ind_FoS_fr[0], ind_FoS_bear[0])))
PlottedResultData = ResultData[plotted_indices]

result_forces = np.column_stack((PlottedResultData['Fx'], PlottedResultData['Fz']))

min_sf = np.inf
min_sf_ansys = None
min_sf_result = None
min_sf_case = None

for fx_i, fz_i, hole_i, direction_i in zip(
    ansys_fx_valid,
    ansys_fz_valid,
    ansys_holes_valid,
    ansys_directions_valid,
):
    d2 = (result_forces[:, 0] - fx_i) ** 2 + (result_forces[:, 1] - fz_i) ** 2
    nearest_idx = np.argmin(d2)
    result_fx, result_fz = result_forces[nearest_idx]

    ans_mag = np.hypot(fx_i, fz_i)
    res_mag = np.hypot(result_fx, result_fz)
    sf = res_mag / ans_mag if ans_mag > 0 else np.inf

    if sf < min_sf:
        min_sf = sf
        min_sf_ansys = (fx_i, fz_i)
        min_sf_result = (result_fx, result_fz)
        min_sf_case = (hole_i, direction_i)

print('Msg: Minimum SF = %.3f' % min_sf)
print('Msg: Case = hole="%s", direction="%s"' % min_sf_case)
print('Msg: Ansys force Fx=%.2f N, Fz=%.2f N' % min_sf_ansys)
print('Msg: Result force Fx=%.2f N, Fz=%.2f N' % min_sf_result)

# =============================================================================
# Plot Patch Failure Envelope
# =============================================================================
S_bearing = 450                                   # [MPa], bearing strength of bulk material
d1 = 15.0 
t = 3.0
F_cr_bearing_wo_patch = S_bearing*d1*t              # [N] Max bearing force of material {d1,t} wo patch

f, ax = plt.subplots(figsize=(10,8))

ax.plot(np.array([1,1])*F_cr_bearing_wo_patch*1.0e-3, [-0.3,0.3], color='m', lw=5, alpha=0.5)


ax.plot(ResultData['Fx'][ind_FoS_sl]*1e-3,ResultData['Fz'][ind_FoS_sl]*1e-3, 'ok', label=r'$FoS_{sl}=[%2.2f; %2.2f]$' %(FoS_low_bound,FoS_up_bound))
ax.plot(ResultData['Fx'][ind_FoS_fr]*1e-3,ResultData['Fz'][ind_FoS_fr]*1e-3, '*r', label=r'$FoS_{fr}=[%2.2f; %2.2f]$' %(FoS_low_bound,FoS_up_bound))
ax.plot(ResultData['Fx'][ind_FoS_bear]*1e-3,ResultData['Fz'][ind_FoS_bear]*1e-3, '>b', label=r'$FoS_{bear}=[%2.2f; %2.2f]$' %(FoS_low_bound,FoS_up_bound))

ax.grid()
ax.set_xlim([-0.5,12.0])
ax.set_ylim([-0.5,5.0])
ax.set_xlabel(r'$F_x, \; [kN]$')
ax.set_ylabel(r'$F_z,\; [kN]$')

# --- Plot forces from OptimizationResults.csv (final_fx, final_fz) ---
# --- Plot forces from Ansys Forces.xlsx ---
if len(ansys_fx) > 0 and len(ansys_fz) > 0:
    ax.plot(ansys_fx*1e-3, ansys_fz*1e-3, 'sg', markersize=6, alpha=0.9, label='Forces from Ansys simulation')

ax.legend(loc='upper right', fontsize=12)

# Try to save as PGF for LaTeX, with fallback to PNG
plt.savefig('Plot_PatchFailureEnvelope.png', bbox_inches='tight')
print('Saved as PNG file')

plt.show()

