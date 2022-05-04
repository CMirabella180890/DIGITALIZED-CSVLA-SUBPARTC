# -*- coding: utf-8 -*-
"""
Created on Sun May  1 07:21:48 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%  GUST ENVELOPE                                                              %
%  NOTE: It's important to remember that in this version of the code the air  %
%        density for all the wind gust calculations are referred to the       %
%        operational altitude of the selected aircraft, whilst the flight     %
%        envelope is referred to Sea Level, Standard Atmosphere air density.  %
%  Vectors with airspeed values                                               %
%  Two vectors with airspeed values within following ranges:                  %
%  1st ---> [0, VC] where VC = MaxCruiseSpeed                                 % 
%  2nd ---> [0, VD] where VD = MaxDiveSpeed                                   %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@author: claum
"""

import numpy as np
import csvla
from numpy.polynomial.chebyshev import chebfit, chebval
obj = csvla.csvla

# AIRCRAFT INPUTS 
N            = 1000 
n_Mass       = aircraft_data.Aircraft_mass["n_Mass"]["Value"]
g            = aircraft_data.Constants["g"]["Value"]
WS           = aircraft_data.Wing_loading["Wing_loading_Pa"]["Value"]
MS           = aircraft_data.Wing_loading["Wing_loading_kg"]["Value"]
rho0         = aircraft_data.Standard_atmosphere["sea_level"]["Density"]["Value"]
rho_op       = aircraft_data.Standard_atmosphere["operative_ceiling"]["Density"]["Value"]
rho_th       = aircraft_data.Standard_atmosphere["theoretical_ceiling"]["Density"]["Value"]
CLalfa_rad   = aircraft_data.Interpolation["CL_alfa_rad"]["Value"]
MGC          = aircraft_data.Wing["MGC"]["Value"]
Ude_cruise   = aircraft_data.Constants["gust_speed_cruise"]["Value"]
Ude_dive     = aircraft_data.Constants["gust_speed_dive"]["Value"]
VC           = aircraft_data.Flight_envelope["Point C"]["VC"]["Value"]
VD           = aircraft_data.Flight_envelope["Point D"]["VD"]["Value"]
VF           = aircraft_data.Flight_envelope["Point F"]["VF"]["Value"]
VE           = aircraft_data.Flight_envelope["Point E"]["VE"]["Value"]
V_fromCtoD   = aircraft_data.Flight_envelope["V_fromCtoD"]["Value"]
V_fromFto0   = aircraft_data.Flight_envelope["V_fromFto0"]["Value"]

# GUST AIRSPEED - CRUISE 
V_gust_cruise = np.ones((N, n_Mass))
V_gust_dive = np.ones((N, n_Mass))
for i in range(0, n_Mass):
    V_gust_cruise[:, i] = np.reshape(np.linspace(0.0, VC[i], N), (N,))
    V_gust_dive[:, i]   = np.reshape(np.linspace(0.0, VD[i], N), (N,))

#############################################
#### GUST PARAMETERS CALCULATIONS - MU_G ####
#############################################
#############################################################################
#################### GUST PARAMETERS CALCULATIONS - K_G #####################
#############################################################################
#############################################################################    
# DENSITY: Sea Level - x = calckg(obj, MassRatio) - GUST ALLEVIATION FACTOR #
# This function calculates the GUST ALLEVIATION FACTOR for the              #
# selected airplane and flight conditions following the CS-VLA              #
# airworthiness prescriprions. To have a complete documentation             #
# check the class fil csvla.m                                               #
#############################################################################
mu_g_op = np.ones((n_Mass,1))
kg_op   = np.ones((n_Mass, 1))
mu_g_sl = np.ones((n_Mass,1))
kg_sl   = np.ones((n_Mass, 1))
for i in range(0, n_Mass): 
    mu_g_op[i] = obj.calcmug(WS[i], MGC, CLalfa_rad, rho_op, g)
    kg_op[i]   = obj.calc_kg(mu_g_op[i])
    mu_g_sl[i] = obj.calcmug(WS[i], MGC, CLalfa_rad, rho0, g)
    kg_sl[i]   = obj.calc_kg(mu_g_sl[i])

#####################################################################
################ GUST PARAMETERS CALCULATIONS - n_G #################
#####################################################################
#####################################################################    
# GUST AIRLOADS CALCULATIONS - OPERATIVE CEILING                    #
# Gust airloads flag                                                #
# --> 'positive'                                                    #
# --> 'negative'                                                    #
# This function is able to calculates in any possible case a vector #
# which contains gust load factors value, following CS-VLA          #
# airworthiness prescription. To have a complete documentation      #
# check the class file csvla.m                                      #
#####################################################################
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++ POSITIVE AND NEGATIVE CRUISE -  CRUISE AND DIVE +++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pos_case_flag    = "positive"
neg_case_flag    = "negative"
ng_pos_cruise_op = np.ones((N, n_Mass))
ng_neg_cruise_op = np.ones((N, n_Mass))
ng_pos_dive_op   = np.ones((N, n_Mass))
ng_neg_dive_op   = np.ones((N, n_Mass))
ng_pos_cruise_sl = np.ones((N, n_Mass))
ng_neg_cruise_sl = np.ones((N, n_Mass))
ng_pos_dive_sl   = np.ones((N, n_Mass))
ng_neg_dive_sl   = np.ones((N, n_Mass))
for i in range(0, n_Mass):
    ng_pos_cruise_op[:,i] = obj.calc_n_gust(rho0, V_gust_cruise[:,i],\
                                            V_gust_dive[:,i], CLalfa_rad,\
                                                kg_op[i], Ude_cruise, WS[i],\
                                                    pos_case_flag)
    ng_neg_cruise_op[:,i] = obj.calc_n_gust(rho0, V_gust_cruise[:,i],\
                                            V_gust_dive[:,i], CLalfa_rad,\
                                                kg_op[i], Ude_cruise, WS[i],\
                                                    neg_case_flag)
    ng_pos_dive_op[:,i] = obj.calc_n_gust(rho0, V_gust_cruise[:,i],\
                                            V_gust_dive[:,i], CLalfa_rad,\
                                                kg_op[i], Ude_dive, WS[i],\
                                                    pos_case_flag)        
    ng_neg_dive_op[:,i] = obj.calc_n_gust(rho0, V_gust_cruise[:,i],\
                                            V_gust_dive[:,i], CLalfa_rad,\
                                                kg_op[i], Ude_dive, WS[i],\
                                                    neg_case_flag)

    ng_pos_cruise_sl[:,i] = obj.calc_n_gust(rho0, V_gust_cruise[:,i],\
                                            V_gust_dive[:,i], CLalfa_rad,\
                                                kg_sl[i], Ude_cruise, WS[i],\
                                                    pos_case_flag)
    ng_neg_cruise_sl[:,i] = obj.calc_n_gust(rho0, V_gust_cruise[:,i],\
                                            V_gust_dive[:,i], CLalfa_rad,\
                                                kg_sl[i], Ude_cruise, WS[i],\
                                                    neg_case_flag)
    ng_pos_dive_sl[:,i] = obj.calc_n_gust(rho0, V_gust_cruise[:,i],\
                                            V_gust_dive[:,i], CLalfa_rad,\
                                                kg_sl[i], Ude_dive, WS[i],\
                                                    pos_case_flag)        
    ng_neg_dive_sl[:,i] = obj.calc_n_gust(rho0, V_gust_cruise[:,i],\
                                            V_gust_dive[:,i], CLalfa_rad,\
                                                kg_sl[i], Ude_dive, WS[i],\
                                                    neg_case_flag)
#############################################################################
######## FINAL STEP - INTERPOLATING BETWEEN CRUISE AND DIVE AIRSPEED ########
#############################################################################
n_gust_fromCtoD = np.ones((N, n_Mass))
n_gust_fromFtoE = np.ones((N, n_Mass))
for i in range(0, n_Mass): 
    x_pos           = np.reshape(np.array([ VC[i], VD[i]] ), (2,))
    y1_pos          = ng_pos_cruise_op[-1, i]
    y2_pos          = ng_pos_dive_op[-1, i]
    y_pos           = np.reshape(np.array([ y1_pos, y2_pos ]), (2,))
    p_gust_positive = chebfit(x_pos, y_pos, deg = 1)
    x_neg           = np.reshape(np.array([ VF[i], VE[i]] ), (2,))
    y1_neg          = ng_neg_cruise_op[-1, i]
    y2_neg          = ng_neg_dive_op[-1, i]
    y_neg           = np.reshape(np.array([ y1_neg, y2_neg ]), (2,))
    p_gust_negative = chebfit(x_neg, y_neg, deg = 1)
    for j in range(0, N):
        n_gust_fromCtoD[j,i] = chebval(V_fromCtoD[j,i], p_gust_positive)
        n_gust_fromFtoE[j,i] = chebval(V_fromFto0[j,i], p_gust_negative)
        
####################################################################
############# PLOTTING RESULTS - GUST ENVELOPE DIAGRAM #############
####################################################################
Reg           = aircraft_data.Aircraft_inputs["Regulation"]["Value"]
Aircraft_name = aircraft_data.Aircraft_inputs["Aircraft_name"]["Value"]
n             = 1
for i in range(0,n_Mass):
    figure_name = "fig" + str(n)
    obj.gust_envelope_diagram(n_from0toS[:,i], n_fromStoA[:,i], n_fromAtoC[:,i], n_fromCtoD[:,i], n_fromDto0[:,i],\
                    V_from0toS[:,i], V_fromStoA[:,i], V_fromAtoC[:,i], V_fromCtoD[:,i], V_fromDto0[:,i],\
                    n_from0toSinv[:,i], n_fromSinvtoG[:,i], n_fromGtoF[:,i], n_fromFto0[:,i],\
                    V_from0toSinv[:,i], V_fromSinvtoG[:,i], V_fromGtoF[:,i], V_fromFto0[:,i],\
                    Reg, Aircraft_name, n,\
                    nS[i], VS[i], nA[i], VA[i], nC[i], VC[i], nD[i], VD[i],\
                    nSinv[i], VSinv[i], nG[i], VG[i], nF[i], VF[i],\
                    ng_pos_cruise_op[:,i], ng_neg_cruise_op[:,i], ng_pos_dive_op[:,i], ng_neg_dive_op[:,i],\
                    n_gust_fromCtoD[:,i], n_gust_fromFtoE[:,i], V_gust_cruise[:,i], V_gust_dive[:,i])
    n = n + 1
# =============================================================================    
# MOVING FIGURES INSIDE OUTPUT FOLDER
# GustEnvelope
n = 1
for i in range(0, n_Mass):
    figure_name = r'\GustEnvelope' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\GustEnvelope'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)
    n = n + 1
# ============================================================================= 
###############################   
###### CLOSE ALL FIGURES ######
###############################
plt.close("all")