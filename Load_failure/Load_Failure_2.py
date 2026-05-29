# -*- coding: utf-8 -*-
"""
Created on thu Api 28 09:13:50 2026

@author: Rasmus Tanghus Thorup
"""
# =============================================================================
# Import Python Packages
# =============================================================================
import os
import time
import numpy as np
from scipy.optimize import least_squares
from FE_input import Class_FE as FE
from Cal_Module import Class_CalModule as cm
import pandas as pd
# Library used to read input values from the command line
import matplotlib.pyplot as plt

from optimize_displacements_predictor import run_insert_model
plt.rcParams.update({'font.size':22})

# =============================================================================
# Defining FE inputs
# =============================================================================
print('Msg: Setting up FE input files')


# Kevlar/Epoxy
MatLam = np.array([112400,112400,112400,0.36,0.36,0.36,64890,64890,64890])      # [E1, E2, E3, v12, v13, v23, G12, G13, G23]

FE(MatLam,[],[]).WriteMatFile()

t = 3.2 #sidewall thickness is 3.2 mm buttomwall thickness is 4 mm

# Kevlar/epoxy
SecProp = np.array([t,'OR_0deg','OR_p45deg','OR_0deg','OR_p45deg','OR_0deg','OR_p45deg','OR_p45deg','OR_0deg','OR_p45deg','OR_0deg','OR_p45deg','OR_0deg']) # 12 layer 3.2mm


FE([],SecProp,[]).WriteSectionFile()

# Specifying Boundary Conditions
BCleft = np.array([1,0,3])      # Left edge, [1,0,3] displacements are locked in the x- and z-directions, y-direction free to move.
BCright = np.array([0,0,0])     # Rigth edge is set to free this is to better simulate the edge of a plate
BCtop = np.array([0,0,3])
BCbot = np.array([0,0,3])




# =============================================================================
# Defining Patch Parameters
# =============================================================================
d0 = 6.0
d1 = 12.0
PatchGeometry = np.array([d0, d1, t])               # 1d array [d0,d1,t]  

n_wind = 4
S_bdl = 700.0                                     # [N], CF 6K bundle, tensile strength of a single bundle (estimated)
S_fr = n_wind*S_bdl                                 # [N], tensile strength of the fibre ring consisting of xx bundles
FibreRingParameters = np.array([n_wind, S_fr])      # 1d array [n_wind,S_fr]    S_fr: strength of fibre ring [N]

NoBundlesInBraid = 16
NoSleeves = 1
theta = np.pi/4.0
SleeveParameters = np.array([NoBundlesInBraid,NoSleeves,theta,S_bdl])      # 1d array [NoBundlesInBraid,NoSleeves,theta,S_bdl]         S_bdl: strength of sleeve bundle [N]

S_bearing = 450                                   # [MPa], Woven Kevlar/Epoxy, bearing strength of bulk material, {Lam-3} 
F_cr_bearing_wo_patch = S_bearing*d1*t              # [N] Max bearing force of material {d1,t} wo patch
LaminateStrengthParameters = np.array([S_bearing])              # 1d array [S_bearing]     S_bearing: bearing strength of laminate [MPa] 


# =============================================================================
# SweepLoad_FailEnvelope Function
# =============================================================================
def read_latest_patch_result_row():
    if not os.path.exists('Result_PatchFailureData.csv'):
        return None

    try:
        patch_results = pd.read_csv('Result_PatchFailureData.csv')
    except Exception:
        return None

    if patch_results.empty:
        return None

    return patch_results.iloc[-1]


def write_sweep_result_csv(order_index, analysis_index, fos_mode, scale, fo_s_value, disp_final, fo_s_error, elapsed_sec, guessed_u=0.0, guessed_w=0.0):
    latest_patch_row = read_latest_patch_result_row()

    if latest_patch_row is not None:
        force_x = latest_patch_row['Fx']
        force_y = latest_patch_row['Fy']
        force_z = latest_patch_row['Fz']
        disp_u = latest_patch_row['u']
        disp_v = latest_patch_row['v']
        disp_w = latest_patch_row['w']
    else:
        force_x = np.nan
        force_y = np.nan
        force_z = np.nan
        disp_u = disp_final[0]
        disp_v = 0.0
        disp_w = disp_final[2]

    # Compute elapsed time since last saved result (if available). Falls back to provided elapsed_sec.
    elapsed_to_write = elapsed_sec
    try:
        if os.path.exists('FoS_Results.csv') and os.path.getsize('FoS_Results.csv') > 0:
            # Only read the timestamp column for speed
            prev = pd.read_csv('FoS_Results.csv', usecols=['timestamp'])
            if not prev.empty:
                prev_ts = pd.to_datetime(prev['timestamp'].iloc[-1])
                now_ts = pd.Timestamp.now()
                elapsed_to_write = (now_ts - prev_ts).total_seconds()
    except Exception:
        # If anything goes wrong reading previous file, keep the provided elapsed_sec
        elapsed_to_write = elapsed_sec

    result_row = pd.DataFrame([{
        'order_index': order_index,
        'analysis_index': analysis_index,
        'fos_mode': fos_mode,
        'scale': scale,
        'guessed_u_mm': guessed_u,
        'guessed_v_mm': 0.0,
        'guessed_w_mm': guessed_w,
        'u_mm': disp_u,
        'v_mm': disp_v,
        'w_mm': disp_w,
        'Fx': force_x,
        'Fy': force_y,
        'Fz': force_z,
        'fos_value': fo_s_value,
        'fos_error': fo_s_error,
        'elapsed_sec': elapsed_to_write,
        'timestamp': pd.Timestamp.now().isoformat()
    }])

    write_header = not os.path.exists('FoS_Results.csv')
    result_row.to_csv('FoS_Results.csv', mode='a', index=False, header=write_header)

# =============================================================================
# SweepLoad_FailEnvelope Function
# =============================================================================

def SweepLoad_FailEnvelope(cal_instance, FoS_map, maxiter, u, w, order_index, analysis_index, sweep_start_time):
    """
    Use least-squares optimisation on a single scale factor.
    u and w are coupled: they scale together while v remains fixed at 0.
    
    Re-edited on Thu Apr 28 09:13:50 2026
    @author: Rasmus Tanghus Thorup
    """
    min_disp = 0.001 # 1 micron
    max_disp = 1000.0# 1 meter  Max and min displacement limits for the optimization to prevent unphysical solutions and excessive runtime. Adjust as needed.

    base_disp = np.array([u, 0, w], dtype=float)

    # Define a function to evaluate the FoS for a given scale factor
    def evaluate_scale(scale):
        disp = base_disp * scale
        Temp = cal_instance.RunAnalysis(disp[0], disp[1], disp[2])# Run analysis with scaled displacements
        FoS_sleeve_bdl, FoS_fr, FoS_BearingStress_w_patch = Temp[0], Temp[1], Temp[2]
        
        if FoS_map == 0: # Map 0: Sleeve BDL, Map 1: Fibre Ring, Map 2: Bearing Stress with patch
            FoS_val = FoS_sleeve_bdl
        elif FoS_map == 1:
            FoS_val = FoS_fr
        elif FoS_map == 2:
            FoS_val = FoS_BearingStress_w_patch
        return FoS_val, disp

    #check if magnitude of base_disp is within reasonable limits.
    scale_lower = 0.0
    scale_upper = np.inf
    magnitude = np.linalg.norm(base_disp)
    if magnitude > 1e-12:
        scale_lower = max(scale_lower, min_disp / magnitude)
        scale_upper = min(scale_upper, max_disp / magnitude)

    eval_counter = 0 # Counter to track number of evaluations in the optimization process
    tolerance_reached = False # Flag to indicate if the FoS is within the desired tolerance, allowing for early exit from optimization
    best_scale = None # Variable to store the best scale factor found when tolerance is reached
    best_FoS = None # Variable to store the best FoS value found when tolerance is reached
    best_disp = None # Variable to store the best displacement vector found when tolerance is reached
    best_error = None # Variable to store the best error value found when tolerance is reached


    # Define the residual function for least_squares optimization
    def residual(scale_array):
        nonlocal eval_counter, tolerance_reached, best_scale, best_FoS, best_disp, best_error

        if tolerance_reached:
            return np.array([0.0])

        eval_counter += 1
        scale = float(scale_array[0])
        FoS_val, disp = evaluate_scale(scale) # Evaluate the FoS for the current scale factor
        error = FoS_val - 1.0 # calulate the error of the cost function (FoS - 1.0) to be minimised

        print('---------------------------------------')
        print('SweepLoad_FailEnvelope evaluation')
        print('---------------------------------------')
        print('Eval: ', eval_counter)
        print('Scale: %2.6f' % scale)
        print('FoS:   %2.6f' % FoS_val)
        print('Error: %2.6f' % error)
        print('u-displacement [mm]: %2.3f' % disp[0])
        print('w-displacement [mm]: %2.3f' % disp[2])

        # Exit early if FoS is within 0.1% of 1.0 (i.e., between 0.999 and 1.001)
        if 0.999 <= FoS_val <= 1.001:
            print('FoS within 0.1% tolerance - exiting early')
            tolerance_reached = True
            best_scale = scale
            best_FoS = FoS_val
            best_disp = disp
            best_error = error
            return np.array([0.0])

        return np.array([error])

    # Run least_squares optimization to find the scale factor that brings the FoS closest to 1.0
    result = least_squares(
        residual,
        x0=np.array([1.0]),
        bounds=(np.array([scale_lower]), np.array([scale_upper])),
        diff_step=0.01,
        xtol=1e-2,
        ftol=1e-2,
        gtol=1e-4,
        max_nfev=maxiter
    )

    if tolerance_reached: # If we exited early due to reaching the FoS tolerance, use the best values found. Otherwise, use the optimization result.
        final_scale = float(best_scale)
        FoS_final = float(best_FoS)
        disp_final = best_disp
        FoS_error = float(best_error**2.0)
        success = True
    else:
        final_scale = float(result.x[0])
        FoS_final, disp_final = evaluate_scale(final_scale)
        FoS_error = (FoS_final - 1.0)**2.0
        success = result.success

    # Print final results and write to CSV if successful and within error tolerance
    if success and FoS_error <= 0.001:
        print('SolvePatchFailEnvelope - exit successful')
        elapsed_sec = time.time() - sweep_start_time
        write_sweep_result_csv(
            order_index,
            analysis_index,
            FoS_map,
            final_scale,
            FoS_final,
            disp_final,
            FoS_error,
            elapsed_sec,
            guessed_u=u,
            guessed_w=w
        )
    else:
        print('SolvePatchFailEnvelope - maxiter reached! - check solution')

    print('FoS_error: ', FoS_error)
    print('FoS:  %2.3f' % FoS_final)
    print('u-displacement [mm]: %2.3f' % disp_final[0])
    print('w-displacement [mm]: %2.3f' % disp_final[2])

    return disp_final



# =============================================================================
# Run Multiple Analyses for determining the patch failure envelope
#     Adaptive Bisection Sampling for all 3 FoS values
# =============================================================================
st = time.time()    # get the start time

if os.path.exists('Result_PatchFailureData.csv'): #Delete previous results file Remove this if you want to keep previous results
    os.remove('Result_PatchFailureData.csv')

if os.path.exists('FoS_Results.csv'):
    os.remove('FoS_Results.csv')

n = 20  # number of analyses >= 2
Disp = 1  # standard displacement magnitude (mm) 


# Select which FoS modes to analyse: 0=FoS_sleeve_bdl, 1=FoS_fr, 2=FoS_bearing
# Modify this list to select which FoS values to evaluate
fos_modes = [0, 1]  # Example: [0, 1, 2] analyses all three, [0] only sleeve, [1, 2] for FR and bearing

# Create an instance of the calculation module with the specified inputs
cal_instance = cm(1,BCleft, BCright, BCtop, BCbot, PatchGeometry, FibreRingParameters, SleeveParameters, LaminateStrengthParameters)

analysis_results = {}

# ------------------------------------------------------------
# Construct linear stiffness matrix for initial guesses
def build_linear_predictor(cal_instance, Disp, n):
    """
    Estimate a 2x2 linearized force-response matrix for the insert model.

    The approximation is:

        [F_inplane]   [K11 K12] [u]   [F_inplane0]
        [Fz       ] = [K21 K22] [w] + [Fz0       ]

    Three CalculiX/Python runs are used:
        1. Baseline: u = 0, w = 0
        2. U perturbation: u = CAL_U, w = 0
        3. W perturbation: u = 0, w = CAL_W
    
    This is not the full FEM stiffness matrix.
    It is a 2x2 linearized response matrix relating [u, w] to [F_inplane, Fz].
    """

    print("\n================================")
    print("Building linear predictor")
    print("================================")

    # Baseline run
    #Temp = cal_instance.RunAnalysis(0.0, 0.0, 0.0)
    #fx0, fy0, fz0 = Temp[3], Temp[4], Temp[5]
    #f_inplane0 = np.sqrt(fx0**2 + fy0**2)
    fx0, fy0, fz0 = 0, 0, 0
    f_inplane0 = 0

    # Perturb u only
    Temp = cal_instance.RunAnalysis(Disp, 0.0, 0.0)
    fx_u, fy_u, fz_u = Temp[3], Temp[4], Temp[5]
    f_inplane_u = np.sqrt(fx_u**2 + fy_u**2)
    
    # Perturb w only
    Temp = cal_instance.RunAnalysis(0.0, 0.0, Disp)
    fx_w, fy_w, fz_w = Temp[3], Temp[4], Temp[5]
    f_inplane_w = np.sqrt(fx_w**2 + fy_w**2)

    # Calculate the stiffness coefficients.
    k11 = (f_inplane_u - f_inplane0) / Disp
    k21 = (fz_u - fz0) / Disp

    k12 = (f_inplane_w - f_inplane0) / Disp
    k22 = (fz_w - fz0) / Disp

    K = np.array([
        [k11, k12],
        [k21, k22]
    ])
    F0 = np.array([f_inplane0,fz0])

    # Calculate the force vectors for the U and W perturbations relative to the baseline.
    f_U = np.array([fz_u- fz0, f_inplane_u - f_inplane0])
    f_W = np.array([fz_w- fz0, f_inplane_w - f_inplane0])

    # Calculate the angles of the force vectors for the U and W perturbations, and the angle between them.
    phi_U = np.arctan2(f_U[0], f_U[1])
    phi_W = np.arctan2(f_W[0], f_W[1])
    phi_W_U = phi_W - phi_U


    # Create a list of n forces that are evenly spaced in angle between f_U and f_W, and linearly spaced in magnitude between the magnitudes of f_U and f_W.
    range_phi = np.linspace(0, phi_W_U, n)+phi_U
    range_F = np.linspace(np.linalg.norm(f_U), np.linalg.norm(f_W), n)
    F_List = np.array([range_F * np.cos(range_phi), range_F * np.sin(range_phi)]).T + F0
    
    # Invert the stiffness matrix to get a linear predictor for displacements from forces.
    K_inv = np.linalg.inv(K)
    list = np.matmul(K_inv, F_List.T)

    #Scale the list of displacements to have a magnitude of Disp.
    for i in range(n):
        list[:,i] = list[:,i]/np.linalg.norm(list[:,i])*Disp
    U_list = list[0,:]
    W_list = list[1,:]


    return U_list, W_list

U_list, W_list = build_linear_predictor(cal_instance,Disp, n)


# ------------------------------------------------------------
# Changing the order of indices to be analysed - using bisection sampling to get a better estimate of the failure envelope with fewer analyses.
def generate_bisection_indices(n):
    """Generate indices in bisection order: start with edges, then midpoint, then fill in"""
    indices = []
    to_process = [(0, n-1)]  # (start_idx, end_idx)
    
    while to_process:
        start, end = to_process.pop(0)
        
        if start not in [idx for idx, _ in indices]:
            indices.append((start, 0))  # Add start
        if end not in [idx for idx, _ in indices]:
            indices.append((end, 0))    # Add end
        
        # Add midpoint and prepare next level
        if end - start > 1:
            mid = (start + end) // 2
            if mid not in [idx for idx, _ in indices]:
                indices.append((mid, 0))
                to_process.append((start, mid))
                to_process.append((mid, end))
    
    return indices

# Generate indices in bisection order
bisection_indices = generate_bisection_indices(n)

# Remove duplicates while preserving order
seen = set()
unique_indices = []
for idx, _ in bisection_indices:
    if idx not in seen:
        seen.add(idx)
        unique_indices.append(idx)

analysis_count = 0

for order_index, idx in enumerate(unique_indices, start=1):
    # Run analysis for all 3 FoS maps
    print('\n========================================')
    print('Analysing index: %i' % idx)
    print('========================================')
    
    fos_results = {}
    for fos_mode in fos_modes:
        u = U_list[idx]
        w = W_list[idx]
        
        # Scale by previous results - for specific FoS value fos_mode
        if len(analysis_results) > 1:  # Only scale if we have at least 8 previous analyses to compare to
            analysed_indices = sorted(analysis_results.keys())
            
            # Find closest index less than current idx (left neighbour)
            left_neighbour = None
            for ai in reversed(analysed_indices):
                if ai < idx:
                    left_neighbour = ai
                    break
            
            # Find closest index greater than current idx (right neighbour)
            right_neighbour = None
            for ai in analysed_indices:
                if ai > idx:
                    right_neighbour = ai
                    break
            
            # Collect scale factors from neighbours (magnitude of displacement / Disp)
            neighbour_scalers = []
            
            # Get scaler from left neighbour
            if left_neighbour is not None and fos_mode in analysis_results[left_neighbour]:
                left_disp = analysis_results[left_neighbour][fos_mode]
                left_magnitude = np.linalg.norm([left_disp[0], left_disp[2]])
                left_scaler = left_magnitude / Disp
                neighbour_scalers.append(left_scaler)

            # Get scaler from right neighbour
            if right_neighbour is not None and fos_mode in analysis_results[right_neighbour]:
                right_disp = analysis_results[right_neighbour][fos_mode]
                right_magnitude = np.linalg.norm([right_disp[0], right_disp[2]])
                right_scaler = right_magnitude / Disp
                neighbour_scalers.append(right_scaler)
            
            # Calculate mean scaling factor from neighbours
            if len(neighbour_scalers) > 0:
                mean_scaler = np.mean(neighbour_scalers)
                u = u * mean_scaler
                w = w * mean_scaler
        
        print('FoS map: %i' % fos_mode)
        print('u-displacement [mm]: %2.3f' % u)
        print('w-displacement [mm]: %2.3f' % w)
        
        result = SweepLoad_FailEnvelope(cal_instance, fos_mode, 20, u, w, order_index, idx, st)
        fos_results[fos_mode] = result
    
    # Store results - only for selected FoS maps
    analysis_results[idx] = fos_results
    analysis_count += len(fos_modes)
    
    print('Analysis %i completed' % analysis_count)
    et = time.time()
    elapsed_time = et - st
    print('Execution time: %3.1f sec. (%3.1f min)' % (elapsed_time, elapsed_time/60.0))


# =============================================================================
# Plot Patch Failure Modes
# =============================================================================
import runpy

runpy.run_path('Plot_PatchFailureEnvelope.py')