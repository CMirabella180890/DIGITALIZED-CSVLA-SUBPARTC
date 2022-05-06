# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 07:51:45 2022

@author: claum
"""

from numpy.polynomial.chebyshev import chebfit, chebval
from ambiance import Atmosphere
import shutil
import numpy as np
import json 
from types import SimpleNamespace
from calc_aero_model import initialization
from calc_aero_model import cl_wb_interpolation, cd_wb_interpolation, cm_wb_interpolation, plot_interpolation_model
from calc_chord_distr import initialization_calc_chord, chord_distribution, mean_aerodynamic_chord
from calc_schrenk import calc_schrenk_distribution, plot_schrenk

#CHORD CALCULATOR
chord_calculator = lambda S_wing, b_wing, taper_ratio, y_halfspan: ( (2 * S_wing) / ( (1 + taper_ratio) * b_wing)) * (1 - (((1 - taper_ratio))/b_wing) * np.abs(2 * y_halfspan))

# ============================================================================
# AIRCRAFT DATABASE INITIALIZATION 
# ============================================================================
# =======================================
JSONFileName = "aircraft_input.json"
with open(JSONFileName, "r") as f:
    # ===================================
    # AIRCRAFT DATA INSIDE A DICTIONARY
    # ===================================
    aircraft = json.load(f)
    # ===================================
# ===================================================
#   AIRCRAFT DATA INSIDE A SIMPLE OBJECT VARIABLE
# ===================================================    
aircraft_data = SimpleNamespace(**aircraft) 
# ===================================================    
# AERODYNAMIC MODEL INTERPOLATION 
# ===================================================
CL_max_clean, CL_max_takeoff, CL_max_landing, CL_inverted, CL_star, CL_wb, CD_wb, CM_wb, AOA_wb = initialization(aircraft_data)
# CURVE FITTING FUNCTION
n_Points = 1000
n_Poly1  = 1
n_Poly2  = 2
# CL INTERPOLATION 
index_cl_star, p_alpha_star, ALPHA_star, alpha_interp, CL_fullmodel, p_cl_wb1, p_cl_wb2 = cl_wb_interpolation(n_Points, n_Poly1, n_Poly2, CL_wb, CL_star, AOA_wb)
# CD INTERPOLATION
alfa0L = p_alpha_star[0]
p_cd_wb, CD_fullmodel, CD0 = cd_wb_interpolation(alpha_interp, alfa0L, n_Poly2, CD_wb, AOA_wb)
# CM INTERPOLATION 
p_CM_wb, CM_fullmodel = cm_wb_interpolation(n_Points, n_Poly1, CL_wb, CM_wb, CL_star, CL_fullmodel, index_cl_star)
CM0         = p_CM_wb[0]
CMCL        = p_CM_wb[1]
CM_alfa_deg = CMCL * p_cl_wb1[1]
CM_alfa_rad = CM_alfa_deg / np.deg2rad(1)
CL_wb_max = np.amax(CL_fullmodel)
# ============================================================
# ================ PLOT INTERPOLATION RESULTS ================ 
# ============================================================
fig1, fig2, fig3 = plot_interpolation_model(AOA_wb, CL_wb, CD_wb, CM_wb,\
                                            alpha_interp, CL_fullmodel,\
                                                CD_fullmodel, CM_fullmodel)
# ============================================================
#                      SAVING FIGURES 
# ============================================================
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++ CHORD DISTRIBUTION AND LOAD DISTRIBUTIONS +++++++++++++++++ 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
b_wing, croot, ctip, S_wing, AR_wing, MGC_wing, chord_kink_one,\
    chord_kink_two, wing_type_flag, taper_ratio1, taper_ratio2, taper_ratio3,\
        ht_span, ht_croot, ht_ctip, ht_taper_ratio, ht_S, ht_AR, ht_MGC,\
            vt_span, vt_croot, vt_ctip, vt_taper_ratio, vt_S, vt_AR, vt_MGC = initialization_calc_chord(aircraft_data)
# =============================================================================
#                             CHORD CALCULATION      
# =============================================================================
N           = 1000       
taper_ratio = ctip / croot
chord, y_halfspan   = chord_distribution(b_wing, S_wing, N, croot, ctip, chord_kink_one,\
                           chord_kink_two, wing_type_flag, taper_ratio1,\
                               taper_ratio2, taper_ratio3) 
vt_y_halfspan       = np.linspace(0.0, vt_span, N)
vt_chord            = chord_calculator(vt_S, vt_span, vt_taper_ratio,\
                                       np.flip(vt_y_halfspan) )
ht_y_halfspan       = np.linspace(0.0, 0.5*ht_span, N)
ht_chord            = chord_calculator(ht_S, ht_span, ht_taper_ratio,\
                                       np.flip(ht_y_halfspan)) 
yMAC_wing, MAC_wing = mean_aerodynamic_chord(S_wing, b_wing, taper_ratio,\
                                             chord, y_halfspan)  
yMAC_ht, MAC_ht     = mean_aerodynamic_chord(ht_S, ht_span, ht_taper_ratio,\
                                             ht_chord, ht_y_halfspan) 
yMAC_vt, MAC_vt     = mean_aerodynamic_chord(vt_S, vt_span, vt_taper_ratio,\
                                             vt_chord, vt_y_halfspan)           
# =============================================================================
#                          SCHRENK LOAD CALCULATION 
# =============================================================================
eta, elliptical_load, schrenk_cCL, unit_CL = calc_schrenk_distribution(b_wing, S_wing, chord, y_halfspan)
# DERIVED DRAG AND PITCHING MOMENT DISTRIBUTIONS 
unit_CD = chebval(unit_CL, p_cd_wb)
unit_CM = chebval(unit_CL, p_CM_wb)
# Plotting results
fig4, fig5, fig6, fig7 = plot_schrenk(y_halfspan, chord, elliptical_load, schrenk_cCL, unit_CL, unit_CD, unit_CM)
# ===================================================    
# STORE INSIDE THE SIMPLE NAMESPACE OBJECT
# ===================================================
Interpolation = { "Interpolation": {       
                                           "CL_alfa_deg"    : {"Value": p_cl_wb1[1],               "Unit": "1/deg"},\
                                           "CL0"            : {"Value": p_cl_wb1[0],               "Unit": "Non dimensional"},\
                                           "CL_alfa_rad"    : {"Value": p_cl_wb1[1]/np.deg2rad(1), "Unit":"1/rad"},\
                                           "CL_fullmodel"   : {"Value": CL_fullmodel,              "Unit": "Non dimensional"},\
                                           "alpha_fullmodel": {"Value": alpha_interp,              "Unit": "deg"},\
                                           "CL_wb_max"      : {"Value": CL_wb_max,                 "Unit": "Non dimensional"},\
                                           "alfa_star"      : {"Value": ALPHA_star,                "Unit": "deg"},\
                                           "alfa0L"         : {"Value": alfa0L,                    "Unit": "deg"},\
                                           "index_cl_star"  : {"Value": index_cl_star,             "Unit": "Pure number"},\
                                           "CD_fullmodel"   : {"Value": CD_fullmodel,              "Unit": "Non dimensional"},\
                                           "CD0"            : {"Value": CD0,                       "Unit": "Non dimensional"},\
                                           "p_cl_wb1"       : {"Value": p_cl_wb1,                  "Unit": "CL interp. coeff."},\
                                           "p_cl_wb2"       : {"Value": p_cl_wb2,                  "Unit": "CL interp. coeff."},\
                                           "p_cd_wb"        : {"Value": p_cd_wb,                   "Unit": "CD interp. coeff."},\
                                           "CM0"            : {"Value": CM0,                       "Unit": "Non dimensional"},\
                                           "CMCL"           : {"Value": CMCL,                      "Unit": "Non dimensional"},\
                                           "CM_alfa_deg"    : {"Value": CM_alfa_deg,               "Unit": "1/deg"},\
                                           "CM_alfa_rad"    : {"Value": CM_alfa_rad,               "Unit": "1/rad"},\
                                           "unit_CL"        : {"Value": unit_CL,                   "Unit": "Non dimensional"},\
                                           "unit_CD"        : {"Value": unit_CD,                   "Unit": "Non dimensional"},\
                                           "unit_CM"        : {"Value": unit_CM,                   "Unit": "Non dimensional"},\
                                           "Schrenk_c_CL"   : {"Value": schrenk_cCL,               "Unit": "m"}    
                                           } 
                 }
# UPDATING THE AIRCRAFT DICTIONARY TYPE    
aircraft.update(Interpolation)
# UPDATING THE AIRCRAFT_DATA SIMPLE NAMESPACE OBJECT
aircraft_data = SimpleNamespace(**aircraft)
# =============================================================================
# =============================================================================
#                          UPDATING THE DICTIONARY 
# =============================================================================
Wing = {"Wing": {'c_root': {'Value': 1.4, 'Unit': 'm'},
 'c_tip': {'Value': 1.4, 'Unit': 'm'},
 'c_kink_one': {'Value': 1.4, 'Unit': 'm'},
 'c_kink_two': {'Value': 1.4, 'Unit': 'm'},
 'b': {'Value': 9.62, 'Unit': 'm'},
 'first_sweep': {'Value': 0.0, 'Unit': 'deg'},
 'second_sweep': {'Value': 0.0, 'Unit': 'deg'},
 'third_sweep': {'Value': 0.0, 'Unit': 'deg'},
 'first_sweep_location': {'Value': 0.0, 'Unit': 'Percentage'},
 'second_sweep_loaction': {'Value': 1.0, 'Unit': 'Percentage'},
 'first_dihedral_angle': {'Value': 1.5, 'Unit': 'deg'},
 'second_dihedral_angle': {'Value': 1.5, 'Unit': 'deg'},
 'third_dihedral_angle': {'Value': 1.5, 'Unit': 'deg'},
 'first_airfoil_section': {'Value': 'GOE_398', 'Unit': 'None'},
 'second_airfoil_section': {'Value': 'GOE_398', 'Unit': 'None'},
 'third_airfoil_section': {'Value': 'GOE_398', 'Unit': 'None'},
 'wing_camber': {'Value': 0.049, 'Unit': 'MAC percentage'},
 'wing_camber_location': {'Value': 0.4, 'Unit': 'MAC percentage'},
 'wing_thickness': {'Value': 0.12, 'Unit': 'MAC percentage'},
 'first_twist_angle': {'Value': -1.5, 'Unit': 'deg'},
 'second_twist_angle': {'Value': -1.5, 'Unit': 'deg'},
 'third_twist_angle': {'Value': -1.5, 'Unit': 'deg'},
 'fourth_twist_angle': {'Value': -1.5, 'Unit': 'deg'},
 'twist_location': {'Value': 0.0, 'Unit': 'Chord percentage'},
 'panel_span_one': {'Value': 0.33, 'Unit': 'Semispan percentage'},
 'panel_span_two': {'Value': 0.33, 'Unit': 'Semispan percentage'},
 'panel_span_third': {'Value': 0.33, 'Unit': 'Semispan percentage'},
 'wing_weight_flag': {'Value': 'Strut_braced_wing', 'Unit': 'flag'},
 'wing_type': {'Value': 'Rectangular', 'Unit': 'flag'},
 'wing_xle': {'Value': 1.638, 'Unit': 'm'},
 'wing_yle': {'Value': 0.0, 'Unit': 'm'},
 'wing_zle': {'Value': 0.7, 'Unit': 'm'}, 
 'S': {'Value': S_wing, 'Unit': 'm^2'},
 'AR': {'Value': AR_wing, 'Unit': 'Non dimensional'},
 'MGC': {'Value': MGC_wing, 'Unit': 'm'}, 
 'taper_ratio1': {'Value': taper_ratio1, 'Unit': 'Non dimensional'},
 'taper_ratio2': {'Value': taper_ratio2, 'Unit': 'Non dimensional'},
 'taper_ratio3': {'Value': taper_ratio3, 'Unit': 'Non dimensional'},
 'chord_distribution': {'Value': chord, 'Unit': 'm'},
 'y_halfspan': {'Value': y_halfspan, 'Unit': 'm'},
 'eta':{'Value': eta, 'Unit': 'Non dimensional'},
 'MAC': {'Value': MAC_wing, 'Unit':'m'}, 
 'MAC_location': {'Value': yMAC_wing, 'Unit': 'm'} } }

aircraft.update(Wing)

Horizontal_Tailplane = {'Horizontal_Tailplane': {'ht_croot': {'Value': 0.68, 'Unit': 'm'},
  'ht_ctip': {'Value': 0.68, 'Unit': 'm'},
  'ht_span': {'Value': 2.9, 'Unit': 'm'},
  'ht_twist': {'Value': 0.0, 'Unit': 'deg'},
  'ht_twist_location': {'Value': 0.25, 'Unit': 'Chord percentage'},
  'ht_sweep': {'Value': 0.0, 'Unit': 'deg'},
  'ht_sweep_location': {'Value': 0.0, 'Unit': 'Chord percentage'},
  'ht_sweep_secondary_location': {'Value': 1.0, 'Unit': 'Chord percentage'},
  'ht_dihedral': {'Value': 0.0, 'Unit': 'deg'},
  'ht_camber': {'Value': 0.0, 'Unit': 'Chord percentage'},
  'ht_camber_location': {'Value': 0.2, 'Unit': 'Chord percentage'},
  'ht_thickness_ratio': {'Value': 0.12, 'Unit': 'Chord percentage'},
  'ht_moment_arm': {'Value': 3.78, 'Unit': 'm'},
  'ht_depsilondalfa': {'Value': 0.3, 'Unit': 'Non dimensional'},
  'ht_xloc': {'Value': 6.0, 'Unit': 'm'},
  'ht_yloc': {'Value': 0.0, 'Unit': 'm'},
  'ht_zloc': {'Value': 0.15, 'Unit': 'm'}, 
  'ht_taper_ratio': {'Value': ht_taper_ratio, 'Unit': 'm'}, 
  'ht_S': {'Value': ht_S, 'Unit': 'm^2'},
  'ht_AR': {'Value': ht_AR, 'Unit': 'Non dimensional'},
  'ht_MGC': {'Value': ht_MGC, 'Unit': 'm'},
  'ht_y_halfspan': {'Value': ht_y_halfspan, 'Unit': 'm'},
  'ht_chord': {'Value': ht_chord, 'Unit': 'm'}, 
  'ht_MAC': {'Value': MAC_ht, 'Unit': 'm'}, 
  'ht_MAC_location': {'Value': yMAC_ht, 'Unit':'m'} } }   

aircraft.update(Horizontal_Tailplane)

Vertical_Tailplane = {'Vertical_Tailplane': 
{'vt_flag': {'Value': 'Single fin', 'Unit': 'flag'},
 'vt_camber': {'Value': 0.0, 'Unit': 'Chord percentage'},
 'vt_camber_location': {'Value': 0.2, 'Unit': 'Chord percentage'},
 'vt_thickness_ratio': {'Value': 0.1, 'Unit': 'Chord percentage'},
 'vt_twist': {'Value': 15, 'Unit': 'deg'},
 'vt_twist_location': {'Value': 0.0, 'Unit': 'Chord percentage'},
 'a_vt_rad': {'Value': 3.65, 'Unit': '1/rad'},
 'vt_span': {'Value': 1.29, 'Unit': 'm'},
 'vt_xloc': {'Value': 5.8, 'Unit': 'm'},
 'vt_yloc': {'Value': 0.0, 'Unit': 'm'},
 'vt_zloc': {'Value': 0.7, 'Unit': 'm'},
 'vt_croot': {'Value': 1.0, 'Unit': 'm'},
 'vt_ctip': {'Value': 0.476, 'Unit': 'm'},
 'vt_sweep': {'Value': 30.0, 'Unit': 'deg'},
 'vt_sweep_location': {'Value': 0.0, 'Unit': 'Chord percentage'},
 'vt_sweep_secondary_location': {'Value': 0.0, 'Unit': 'Chord percentage'},
 'vt_dihedral': {'Value': 0.0, 'Unit': 'deg'},
 'vt_taper_ratio': {'Value': vt_taper_ratio, 'Unit': 'm'}, 
 'vt_S': {'Value': vt_S, 'Unit': 'm^2'},
 'vt_AR': {'Value': vt_AR, 'Unit': 'Non dimensional'},
 'vt_MGC': {'Value': vt_MGC, 'Unit': 'm'},
 'vt_y_halfspan': {'Value': vt_y_halfspan, 'Unit': 'm'},
 'vt_chord': {'Value': vt_chord, 'Unit': 'm'}, 
 'vt_MAC': {'Value': MAC_vt, 'Unit': 'm'}, 
 'vt_MAC_location': {'Value': yMAC_vt, 'Unit':'m'} }}
    
aircraft.update(Vertical_Tailplane)    

# MOVING FIGURES INSIDE OUTPUT FOLDER
# FIGURE 1 
original = r'I:\PythonTesiConversion\initialization\figure1.pdf'
target = r'I:\PythonTesiConversion\initialization\Output\Interpolation\figure1.pdf'
shutil.move(original, target)
# FIGURE 2 
original = r'I:\PythonTesiConversion\initialization\figure2.pdf'
target = r'I:\PythonTesiConversion\initialization\Output\Interpolation\figure2.pdf'
shutil.move(original, target)
# FIGURE 3 
original = r'I:\PythonTesiConversion\initialization\figure3.pdf'
target = r'I:\PythonTesiConversion\initialization\Output\Interpolation\figure3.pdf'
shutil.move(original, target)
# FIGURE 4 
original = r'I:\PythonTesiConversion\initialization\figure4.pdf'
target = r'I:\PythonTesiConversion\initialization\Output\Schrenk\figure4.pdf'
shutil.move(original, target)
# FIGURE 5 
original = r'I:\PythonTesiConversion\initialization\figure5.pdf'
target = r'I:\PythonTesiConversion\initialization\Output\Schrenk\figure5.pdf'
shutil.move(original, target)
# FIGURE 6 
original = r'I:\PythonTesiConversion\initialization\figure6.pdf'
target = r'I:\PythonTesiConversion\initialization\Output\Schrenk\figure6.pdf'
shutil.move(original, target)
# FIGURE 7 
original = r'I:\PythonTesiConversion\initialization\figure7.pdf'
target = r'I:\PythonTesiConversion\initialization\Output\Schrenk\figure7.pdf'
shutil.move(original, target)
# =============================================================================
#            +++++++++++++++++++++++++++++++++++++++++++++++++
#            ++++++++++++++ STANDARD ATMOSPHERE ++++++++++++++
#            +++++++++++++++++++++++++++++++++++++++++++++++++
# =============================================================================
sea_level           = aircraft_data.Standard_atmosphere["sea_level"]["h"]["Value"]
operative_ceiling   = aircraft_data.Standard_atmosphere["operative_ceiling"]["h"]["Value"]
theoretical_ceiling = aircraft_data.Standard_atmosphere["theoretical_ceiling"]["h"]["Value"]
ISA_sea_level       = Atmosphere(sea_level)
ISA_operative       = Atmosphere(operative_ceiling)
ISA_theoretical     = Atmosphere(theoretical_ceiling)

Standard_Atmosphere = {'Standard_atmosphere': {
    'sea_level': {
        "h" : {"Value": sea_level, "Unit": "m"},
        'Density' : {'Value': ISA_sea_level.density, 'Unit': 'kg/m^3'},
        'Pressure': {'Value': ISA_sea_level.pressure, 'Unit': 'Pa'}, 
        'Temperature': {'Value': ISA_sea_level.temperature, 'Unit': 'K'} 
        }, 
    'operative_ceiling': {
        "h" : {"Value": operative_ceiling, "Unit": "m"},
        'Density' : {'Value': ISA_operative.density, 'Unit': 'kg/m^3'},
        'Pressure': {'Value': ISA_operative.pressure, 'Unit': 'Pa'}, 
        'Temperature': {'Value': ISA_operative.temperature, 'Unit': 'K'} 
        }, 
    'theoretical_ceiling': {
        "h" : {"Value": theoretical_ceiling, "Unit": "m"},
        'Density' : {'Value': ISA_theoretical.density, 'Unit': 'kg/m^3'},
        'Pressure': {'Value': ISA_theoretical.pressure, 'Unit': 'Pa'}, 
        'Temperature': {'Value': ISA_theoretical.temperature, 'Unit': 'K'} } },
    'Constants': { 
        'g'                 : {'Value': ISA_sea_level.grav_accel, 'Unit': 'm/s^2'},
		"gust_speed_cruise" : {"Value": 15.24,   "Unit": "m/s"},
		"gust_speed_dive"   : {"Value": 7.62,    "Unit": "m/s"} 
                  } }
aircraft.update(Standard_Atmosphere)   

# UPDATING THE SIMPLE NAMESPACE OBJECT
aircraft_data = SimpleNamespace(**aircraft)