# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:29:25 2022

@author: claum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALANCING LOADS                                                               %
% Calculation of the aerodynamic balancing loads; these loads are related to an % 
% unaccelerated flight condition (n=1), straight or inverted. Wing body aero-   %
% dynamic data are used to assess the horizontal tailplane balancing loads.     %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flight_env_case_pos: 
    ---------------> Case1 
    ---------------> Case2 
    ---------------> Case3 
flight_env_case_neg: 
    ---------------> Case1_inverted 
    ---------------> Case2_inverted 
    ---------------> Case3_inverted
    
"""

import numpy as np
import balancing_loads
obj1 = balancing_loads.balancing_loads

flag2 = "deg"
CL_max_clean = max(CL_fullmodel)
CM_gear      = aircraft_data.Wing_body_aerodynamic["CM_landing_gear"]["Value"]


CM_CG_from0toS     = np.ones((N, n_Mass))
CM_CG_fromStoA     = np.ones((N, n_Mass))
CM_CG_fromAtoGust1 = np.ones((N, n_Mass))
CM_CG_fromGust1toC = np.ones((N, n_Mass))
CM_CG_fromCtoGust2 = np.ones((N, n_Mass))
CM_CG_fromGust2toD = np.ones((N, n_Mass))
CM_CG_fromDto0     = np.ones((N, n_Mass))
CM_CG_from0toSinv  = np.ones((N, n_Mass))
CM_CG_fromSinvtoG  = np.ones((N, n_Mass))
CM_CG_fromGtoGust1 = np.ones((N, n_Mass))
CM_CG_fromGust1toF = np.ones((N, n_Mass))
CM_CG_fromFtoE     = np.ones((N, n_Mass))
CM_CG_fromEto0     = np.ones((N, n_Mass))
CM_CG_fromGtoF     = np.ones((N, n_Mass))
CM_CG_fromAtoC     = np.ones((N, n_Mass))
CM_CG_fromCtoD     = np.ones((N, n_Mass))
   
CM_due_to_CD_from0toS     = np.ones((N, n_Mass))
CM_due_to_CD_fromStoA     = np.ones((N, n_Mass))
CM_due_to_CD_fromAtoGust1 = np.ones((N, n_Mass))
CM_due_to_CD_fromGust1toC = np.ones((N, n_Mass))
CM_due_to_CD_fromCtoGust2 = np.ones((N, n_Mass))
CM_due_to_CD_fromGust2toD = np.ones((N, n_Mass))
CM_due_to_CD_fromDto0     = np.ones((N, n_Mass))
CM_due_to_CD_from0toSin   = np.ones((N, n_Mass))
CM_due_to_CD_fromSinvtoG  = np.ones((N, n_Mass))
CM_due_to_CD_fromGtoGust1 = np.ones((N, n_Mass))
CM_due_to_CD_fromGust1toF = np.ones((N, n_Mass))
CM_due_to_CD_fromFtoE     = np.ones((N, n_Mass))
CM_due_to_CD_fromEto0     = np.ones((N, n_Mass))
CM_due_to_CD_fromGtoF     = np.ones((N, n_Mass))
CM_due_to_CD_fromAtoC     = np.ones((N, n_Mass))
CM_due_to_CD_fromCtoD     = np.ones((N, n_Mass))
   
CM_due_to_CL_from0toS     = np.ones((N, n_Mass))
CM_due_to_CL_fromStoA     = np.ones((N, n_Mass))
CM_due_to_CL_fromAtoGust1 = np.ones((N, n_Mass))
CM_due_to_CL_fromGust1toC = np.ones((N, n_Mass))
CM_due_to_CL_fromCtoGust2 = np.ones((N, n_Mass))
CM_due_to_CL_fromGust2toD = np.ones((N, n_Mass))
CM_due_to_CL_fromDto0     = np.ones((N, n_Mass))
CM_due_to_CL_from0toSinv  = np.ones((N, n_Mass))
CM_due_to_CL_fromSinvtoG  = np.ones((N, n_Mass))
CM_due_to_CL_fromGtoGust1 = np.ones((N, n_Mass))
CM_due_to_CL_fromGust1toF = np.ones((N, n_Mass))
CM_due_to_CL_fromFtoE     = np.ones((N, n_Mass))
CM_due_to_CL_fromEto0     = np.ones((N, n_Mass))
CM_due_to_CL_fromGtoF     = np.ones((N, n_Mass))
CM_due_to_CL_fromAtoC     = np.ones((N, n_Mass))
CM_due_to_CL_fromCtoD     = np.ones((N, n_Mass))
   
q_from0toS     = np.ones((N, n_Mass))
q_fromStoA     = np.ones((N, n_Mass))
q_fromAtoGust1 = np.ones((N, n_Mass))
q_fromGust1toC = np.ones((N, n_Mass))
q_fromCtoGust2 = np.ones((N, n_Mass))
q_fromGust2toD = np.ones((N, n_Mass))
q_fromDto0     = np.ones((N, n_Mass))
q_from0toSinv  = np.ones((N, n_Mass))
q_fromSinvtoG  = np.ones((N, n_Mass))
q_fromGtoGust1 = np.ones((N, n_Mass))
q_fromGust1toF = np.ones((N, n_Mass))
q_fromFtoE     = np.ones((N, n_Mass))
q_fromEto0     = np.ones((N, n_Mass))
q_fromGtoF     = np.ones((N, n_Mass))
q_fromAtoC     = np.ones((N, n_Mass))
q_fromCtoD     = np.ones((N, n_Mass))
   
CD_from0toS     = np.ones((N, n_Mass))
CD_fromStoA     = np.ones((N, n_Mass))
CD_fromAtoGust1 = np.ones((N, n_Mass))
CD_fromGust1toC = np.ones((N, n_Mass))
CD_fromCtoGust2 = np.ones((N, n_Mass))
CD_fromGust2toD = np.ones((N, n_Mass))
CD_fromDto0     = np.ones((N, n_Mass))
CD_from0toSinv  = np.ones((N, n_Mass))
CD_fromSinvtoG  = np.ones((N, n_Mass))
CD_fromGtoGust1 = np.ones((N, n_Mass))
CD_fromGust1toF = np.ones((N, n_Mass))
CD_fromFtoE     = np.ones((N, n_Mass))
CD_fromEto0     = np.ones((N, n_Mass))
CD_fromGtoF     = np.ones((N, n_Mass))
CD_fromAtoC     = np.ones((N, n_Mass))
CD_fromCtoD     = np.ones((N, n_Mass))
   
alfa_new_from0toS     = np.ones((N, n_Mass))
alfa_new_fromStoA     = np.ones((N, n_Mass))
alfa_new_fromAtoGust1 = np.ones((N, n_Mass))
alfa_new_fromGust1toC = np.ones((N, n_Mass))
alfa_new_fromCtoGust2 = np.ones((N, n_Mass))
alfa_new_fromGust2toD = np.ones((N, n_Mass))
alfa_new_fromDto0     = np.ones((N, n_Mass))
lfa_new_from0toSinv   = np.ones((N, n_Mass))
alfa_new_fromSinvtoG  = np.ones((N, n_Mass))
alfa_new_fromGtoGust1 = np.ones((N, n_Mass))
alfa_new_fromGust1toF = np.ones((N, n_Mass))
alfa_new_fromFtoE     = np.ones((N, n_Mass))
alfa_new_fromEto0     = np.ones((N, n_Mass))
alfa_new_fromGtoF     = np.ones((N, n_Mass))
alfa_new_fromAtoC     = np.ones((N, n_Mass))
alfa_new_fromCtoD     = np.ones((N, n_Mass))
   
alfa_from0toS     = np.ones((N, n_Mass))
alfa_fromStoA     = np.ones((N, n_Mass))
alfa_fromAtoGust1 = np.ones((N, n_Mass))
alfa_fromGust1toC = np.ones((N, n_Mass))
alfa_fromCtoGust2 = np.ones((N, n_Mass))
alfa_fromGust2toD = np.ones((N, n_Mass))
alfa_fromDto0     = np.ones((N, n_Mass))
alfa_from0toSinv  = np.ones((N, n_Mass))
alfa_fromSinvtoG  = np.ones((N, n_Mass))
alfa_fromGtoGust1 = np.ones((N, n_Mass))
alfa_fromGust1toF = np.ones((N, n_Mass))
alfa_fromFtoE     = np.ones((N, n_Mass))
alfa_fromEto0     = np.ones((N, n_Mass))
alfa_fromGtoF     = np.ones((N, n_Mass))
alfa_fromAtoC     = np.ones((N, n_Mass))
alfa_fromCtoD     = np.ones((N, n_Mass))

CL_wb_from0toS     = np.ones((N, n_Mass))
CL_wb_fromStoA     = np.ones((N, n_Mass))
CL_wb_fromAtoGust1 = np.ones((N, n_Mass))
CL_wb_fromGust1toC = np.ones((N, n_Mass))
CL_wb_fromCtoGust2 = np.ones((N, n_Mass))
CL_wb_fromGust2toD = np.ones((N, n_Mass))
CL_wb_fromDto0     = np.ones((N, n_Mass))
CL_wb_from0toSinv  = np.ones((N, n_Mass))
CL_wb_fromSinvtoG  = np.ones((N, n_Mass))
CL_wb_fromGtoGust1 = np.ones((N, n_Mass))
CL_wb_fromGust1toF = np.ones((N, n_Mass))
CL_wb_fromFtoE     = np.ones((N, n_Mass))
CL_wb_fromEto0     = np.ones((N, n_Mass))
CL_wb_fromGtoF     = np.ones((N, n_Mass))
CL_wb_fromAtoC     = np.ones((N, n_Mass))
CL_wb_fromCtoD     = np.ones((N, n_Mass))

CL_ht_from0toS     = np.ones((N, n_Mass))
CL_ht_fromStoA     = np.ones((N, n_Mass))
CL_ht_fromAtoGust1 = np.ones((N, n_Mass))
CL_ht_fromGust1toC = np.ones((N, n_Mass))
CL_ht_fromCtoGust2 = np.ones((N, n_Mass))
CL_ht_fromGust2toD = np.ones((N, n_Mass))
CL_ht_fromDto0     = np.ones((N, n_Mass))
CL_ht_from0toSinv  = np.ones((N, n_Mass))
CL_ht_fromSinvtoG  = np.ones((N, n_Mass))
CL_ht_fromGtoGust1 = np.ones((N, n_Mass))
CL_ht_fromGust1toF = np.ones((N, n_Mass))
CL_ht_fromFtoE     = np.ones((N, n_Mass))
CL_ht_fromEto0     = np.ones((N, n_Mass))
CL_ht_fromGtoF     = np.ones((N, n_Mass))
CL_ht_fromAtoC     = np.ones((N, n_Mass))
CL_ht_fromCtoD     = np.ones((N, n_Mass))

CL_new_from0toS     = np.ones((N, n_Mass))
CL_new_fromStoA     = np.ones((N, n_Mass))
CL_new_fromAtoGust1 = np.ones((N, n_Mass))
CL_new_fromGust1toC = np.ones((N, n_Mass))
CL_new_fromCtoGust2 = np.ones((N, n_Mass))
CL_new_fromGust2toD = np.ones((N, n_Mass))
CL_new_fromDto0     = np.ones((N, n_Mass))
CL_new_from0toSinv  = np.ones((N, n_Mass))
CL_new_fromSinvtoG  = np.ones((N, n_Mass))
CL_new_fromGtoGust1 = np.ones((N, n_Mass))
CL_new_fromGust1toF = np.ones((N, n_Mass))
CL_new_fromFtoE     = np.ones((N, n_Mass))
CL_new_fromEto0     = np.ones((N, n_Mass))
CL_new_fromGtoF     = np.ones((N, n_Mass))
CL_new_fromAtoC     = np.ones((N, n_Mass))
CL_new_fromCtoD     = np.ones((N, n_Mass))
                           
L_wb_from0toS      = np.ones((N, n_Mass))
L_new_from0toS     = np.ones((N, n_Mass))
L_ht_from0toS      = np.ones((N, n_Mass))
L_wb_fromStoA      = np.ones((N, n_Mass))
L_new_fromStoA     = np.ones((N, n_Mass)) 
L_ht_fromStoA      = np.ones((N, n_Mass)) 
L_new_fromAtoGust1 = np.ones((N, n_Mass))
L_wb_fromAtoGust1  = np.ones((N, n_Mass))
L_ht_fromAtoGust1  = np.ones((N, n_Mass))
L_new_fromGust1toC = np.ones((N, n_Mass))
L_wb_fromGust1toC  = np.ones((N, n_Mass))
L_ht_fromGust1toC  = np.ones((N, n_Mass))
L_new_fromCtoGust2 = np.ones((N, n_Mass))
L_wb_fromCtoGust2  = np.ones((N, n_Mass))
L_ht_fromCtoGust2  = np.ones((N, n_Mass))
L_new_fromGust2toD = np.ones((N, n_Mass))
L_wb_fromGust2toD  = np.ones((N, n_Mass))
L_ht_fromGust2toD  = np.ones((N, n_Mass))
L_new_fromDto0     = np.ones((N, n_Mass))
L_wb_fromDto0      = np.ones((N, n_Mass)) 
L_ht_fromDto0      = np.ones((N, n_Mass)) 
L_new_from0toSinv  = np.ones((N, n_Mass))
L_wb_from0toSinv   = np.ones((N, n_Mass))
L_ht_from0toSinv   = np.ones((N, n_Mass))
L_new_fromSinvtoG  = np.ones((N, n_Mass)) 
L_wb_fromSinvtoG   = np.ones((N, n_Mass))
L_ht_fromSinvtoG   = np.ones((N, n_Mass))
L_new_fromGtoGust1 = np.ones((N, n_Mass))
L_wb_fromGtoGust1  = np.ones((N, n_Mass))
L_ht_fromGtoGust1  = np.ones((N, n_Mass))
L_new_fromGust1toF = np.ones((N, n_Mass))
L_wb_fromGust1toF  = np.ones((N, n_Mass))
L_ht_fromGust1toF  = np.ones((N, n_Mass))
L_new_fromFtoE     = np.ones((N, n_Mass))
L_wb_fromFtoE      = np.ones((N, n_Mass))
L_ht_fromFtoE      = np.ones((N, n_Mass))
L_new_fromEto0     = np.ones((N, n_Mass))
L_wb_fromEto0      = np.ones((N, n_Mass))
L_ht_fromEto0      = np.ones((N, n_Mass))
L_new_fromGtoF     = np.ones((N, n_Mass))
L_wb_fromGtoF      = np.ones((N, n_Mass))
L_ht_fromGtoF      = np.ones((N, n_Mass))
L_new_fromAtoC     = np.ones((N, n_Mass))
L_wb_fromAtoC      = np.ones((N, n_Mass))
L_ht_fromAtoC      = np.ones((N, n_Mass))
L_new_fromCtoD     = np.ones((N, n_Mass))
L_wb_fromCtoD      = np.ones((N, n_Mass))
L_ht_fromCtoD      = np.ones((N, n_Mass))
###########################################
############# UNIT LOAD FACTOR ############
###########################################
V_unit_load_factor        = np.ones((N, n_Mass))
n_unit_load_factor        = np.ones((N, n_Mass))
CL_unit_load_factor       = np.ones((N, n_Mass))
alfa_unit_load_factor     = np.ones((N, n_Mass))
alfa_new_unit_load_factor = np.ones((N, n_Mass))
CD_unit_load_factor       = np.ones((N, n_Mass))
q_unit_load_factor        = np.ones((N, n_Mass))
CM_CL_unit_load_factor    = np.ones((N, n_Mass))
CM_CD_unit_load_factor    = np.ones((N, n_Mass))
CM_CT_unit_load_factor    = np.ones((N, n_Mass))
CM_CG_unit_load_factor    = np.ones((N, n_Mass))
CL_ht_unit_load_factor    = np.ones((N, n_Mass))
CL_new_unit_load_factor   = np.ones((N, n_Mass))
L_wb_unit_load_factor     = np.ones((N, n_Mass))
L_new_unit_load_factor    = np.ones((N, n_Mass))
L_ht_unit_load_factor     = np.ones((N, n_Mass)) 

flag1 = "straight" 
flag2 = "deg"          
# ---------------------------------------------------------------------------------- #
# ------------------------- UNIT LOAD FACTOR CALCULATION --------------------------- #
# ---------------------------------------------------------------------------------- #
for j in range(0, n_Mass):
    V_unit_load_factor[0,j] = VS[j]
    for i in range(1, N):
        V_unit_load_factor[i,j] = V_unit_load_factor[i-1, j] + ( VD[j] - VS[j] ) * ( 1 / len(V_unit_load_factor[:,j]) );
#####################################################################################
################################# CALCULATIONS ######################################
#####################################################################################
for j in range(0, n_Mass):
    for i in range(0, N):
        CL_unit_load_factor[i,j], alfa_unit_load_factor[i,j], alfa_new_unit_load_factor[i,j],\
        CD_unit_load_factor[i,j], q_unit_load_factor[i,j], CM_CL_unit_load_factor[i,j],\
        CM_CD_unit_load_factor[i,j], CM_CT_unit_load_factor[i,j],\
        CM_CG_unit_load_factor[i,j], CL_ht_unit_load_factor[i,j], CL_new_unit_load_factor[i,j],\
        L_wb_unit_load_factor[i,j], L_new_unit_load_factor[i,j], L_ht_unit_load_factor[i,j] = obj1.flight_balancing_parameters(rho0, V_unit_load_factor[i,j], WS[j], n_unit_load_factor[i,j], CL_max_clean, alfa0L,\
                                        p_cd_wb, CL_star, p_cl_wb1,\
                                        p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                        bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                        CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                        flag1, flag2, p_CM_wb[1], obj1)           
###########################################
############# STRAIGHT FLIGHT #############   
########################################### 
if ( flight_env_case_pos == "Case1" ): 
    
    flag1 = "straight"
    
    # V_fe_from0toS     --- n_fe_from0toS
    CL_wb_from0toS        = np.ones((N, n_Mass)) 
    alfa_from0toS         = np.ones((N, n_Mass)) 
    alfa_new_from0toS     = np.ones((N, n_Mass)) 
    CD_from0toS           = np.ones((N, n_Mass)) 
    q_from0toS            = np.ones((N, n_Mass)) 
    CM_due_to_CL_from0toS = np.ones((N, n_Mass))
    CM_due_to_CD_from0toS = np.ones((N, n_Mass))
    CM_due_to_CT_from0toS = np.ones((N, n_Mass))
    CM_CG_from0toS        = np.ones((N, n_Mass))
    CL_ht_from0toS        = np.ones((N, n_Mass))
    CL_new_from0toS       = np.ones((N, n_Mass))
    L_wb_from0toS         = np.ones((N, n_Mass))
    L_new_from0toS        = np.ones((N, n_Mass))
    L_ht_from0toS         = np.ones((N, n_Mass))
    
    # BALANCING LOADS INPUTS PARAMETERS
    #  _________________________________________
    # | rho, V, WS, n, CLmax_aero_model, alfa0l,|
    # |p_drag, CL_star, CL0, CLalfa, p_lift1,   |
    # |p_lift2, xAC, xCG, bCG, MAC, h, CM0,     |
    # |CM_gear, S, l_ht, flag, CM_CL, obj1      |
    # |_________________________________________|

    # BALANCING LOADS CALCULATIONS - FROM 0 TO S
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_from0toS[i,j], alfa_from0toS[i,j], alfa_new_from0toS[i,j],\
            CD_from0toS[i,j], q_from0toS[i,j], CM_due_to_CL_from0toS[i,j],\
            CM_due_to_CD_from0toS[i,j], CM_due_to_CT_from0toS[i,j],\
            CM_CG_from0toS[i,j], CL_ht_from0toS[i,j], CL_new_from0toS[i,j],\
            L_wb_from0toS[i,j], L_new_from0toS[i,j], L_ht_from0toS[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_from0toS[i,j], WS[j], n_fe_from0toS[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1)       
    
    # V_fe_fromStoA     --- n_fe_fromStoA
    CL_wb_fromStoA        = np.ones((N, n_Mass)) 
    alfa_fromStoA         = np.ones((N, n_Mass))
    alfa_new_fromStoA     = np.ones((N, n_Mass)) 
    CD_fromStoA           = np.ones((N, n_Mass)) 
    q_fromStoA            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromStoA = np.ones((N, n_Mass))
    CM_due_to_CD_fromStoA = np.ones((N, n_Mass))
    CM_due_to_CT_fromStoA = np.ones((N, n_Mass))
    CM_CG_fromStoA        = np.ones((N, n_Mass))
    CL_ht_fromStoA        = np.ones((N, n_Mass))
    CL_new_fromStoA       = np.ones((N, n_Mass))
    L_wb_fromStoA         = np.ones((N, n_Mass))
    L_new_fromStoA        = np.ones((N, n_Mass))
    L_ht_fromStoA         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM S TO A
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromStoA[i,j], alfa_fromStoA[i,j], alfa_new_fromStoA[i,j],\
            CD_fromStoA[i,j], q_fromStoA[i,j], CM_due_to_CL_fromStoA[i,j],\
            CM_due_to_CD_fromStoA[i,j], CM_due_to_CT_fromStoA[i,j],\
            CM_CG_fromStoA[i,j], CL_ht_fromStoA[i,j], CL_new_fromStoA[i,j],\
            L_wb_fromStoA[i,j], L_new_fromStoA[i,j], L_ht_fromStoA[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromStoA[i,j], WS[j], n_fe_fromStoA[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1)    

    # V_fe_fromAtoGust1 --- n_fe_fromAtoGust1
    CL_wb_fromAtoGust1        = np.ones((N, n_Mass)) 
    alfa_fromAtoGust1         = np.ones((N, n_Mass)) 
    alfa_new_fromAtoGust1     = np.ones((N, n_Mass)) 
    CD_fromAtoGust1           = np.ones((N, n_Mass)) 
    q_fromAtoGust1            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromAtoGust1 = np.ones((N, n_Mass))
    CM_due_to_CD_fromAtoGust1 = np.ones((N, n_Mass))
    CM_due_to_CT_fromAtoGust1 = np.ones((N, n_Mass))
    CM_CG_fromAtoGust1        = np.ones((N, n_Mass))
    CL_ht_fromAtoGust1        = np.ones((N, n_Mass))
    CL_new_fromAtoGust1       = np.ones((N, n_Mass))
    L_wb_fromAtoGust1         = np.ones((N, n_Mass))
    L_new_fromAtoGust1        = np.ones((N, n_Mass))
    L_ht_fromAtoGust1         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM A TO GUST1
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromAtoGust1[i,j], alfa_fromAtoGust1[i,j], alfa_new_fromAtoGust1[i,j],\
            CD_fromAtoGust1[i,j], q_fromAtoGust1[i,j], CM_due_to_CL_fromAtoGust1[i,j],\
            CM_due_to_CD_fromAtoGust1[i,j], CM_due_to_CT_fromAtoGust1[i,j],\
            CM_CG_fromAtoGust1[i,j], CL_ht_fromAtoGust1[i,j], CL_new_fromAtoGust1[i,j],\
            L_wb_fromAtoGust1[i,j], L_new_fromAtoGust1[i,j], L_ht_fromAtoGust1[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromAtoGust1[i,j], WS[j], n_fe_fromAtoGust1[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 

    # V_fe_fromGust1toC --- n_fe_fromGust1toC
    CL_wb_fromGust1toC        = np.ones((N, n_Mass)) 
    alfa_fromGust1toC         = np.ones((N, n_Mass)) 
    alfa_new_fromGust1toC     = np.ones((N, n_Mass)) 
    CD_fromGust1toC           = np.ones((N, n_Mass)) 
    q_fromGust1toC            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromGust1toC = np.ones((N, n_Mass))
    CM_due_to_CD_fromGust1toC = np.ones((N, n_Mass))
    CM_due_to_CT_fromGust1toC = np.ones((N, n_Mass))
    CM_CG_fromGust1toC        = np.ones((N, n_Mass))
    CL_ht_fromGust1toC        = np.ones((N, n_Mass))
    CL_new_fromGust1toC       = np.ones((N, n_Mass))
    L_wb_fromGust1toC         = np.ones((N, n_Mass))
    L_new_fromGust1toC        = np.ones((N, n_Mass))
    L_ht_fromGust1toC         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM GUST1 TO C
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromGust1toC[i,j], alfa_fromGust1toC[i,j], alfa_new_fromGust1toC[i,j],\
            CD_fromGust1toC[i,j], q_fromGust1toC[i,j], CM_due_to_CL_fromGust1toC[i,j],\
            CM_due_to_CD_fromGust1toC[i,j], CM_due_to_CT_fromGust1toC[i,j],\
            CM_CG_fromGust1toC[i,j], CL_ht_fromGust1toC[i,j], CL_new_fromGust1toC[i,j],\
            L_wb_fromGust1toC[i,j], L_new_fromGust1toC[i,j], L_ht_fromGust1toC[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromGust1toC[i,j], WS[j], n_fe_fromGust1toC[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromCtoGust2 --- n_fe_fromCtoGust2
    CL_wb_fromCtoGust2        = np.ones((N, n_Mass)) 
    alfa_fromCtoGust2         = np.ones((N, n_Mass)) 
    alfa_new_fromCtoGust2     = np.ones((N, n_Mass)) 
    CD_fromCtoGust2           = np.ones((N, n_Mass)) 
    q_fromCtoGust2            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromCtoGust2 = np.ones((N, n_Mass))
    CM_due_to_CD_fromCtoGust2 = np.ones((N, n_Mass))
    CM_due_to_CT_fromCtoGust2 = np.ones((N, n_Mass))
    CM_CG_fromCtoGust2        = np.ones((N, n_Mass))
    CL_ht_fromCtoGust2        = np.ones((N, n_Mass))
    CL_new_fromCtoGust2       = np.ones((N, n_Mass))
    L_wb_fromCtoGust2         = np.ones((N, n_Mass))
    L_new_fromCtoGust2        = np.ones((N, n_Mass))
    L_ht_fromCtoGust2         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM C TO GUST2
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromCtoGust2[i,j], alfa_fromCtoGust2[i,j], alfa_new_fromCtoGust2[i,j],\
            CD_fromCtoGust2[i,j], q_fromCtoGust2[i,j], CM_due_to_CL_fromCtoGust2[i,j],\
            CM_due_to_CD_fromCtoGust2[i,j], CM_due_to_CT_fromCtoGust2[i,j],\
            CM_CG_fromCtoGust2[i,j], CL_ht_fromCtoGust2[i,j], CL_new_fromCtoGust2[i,j],\
            L_wb_fromCtoGust2[i,j], L_new_fromCtoGust2[i,j], L_ht_fromCtoGust2[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromCtoGust2[i,j], WS[j], n_fe_fromCtoGust2[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromGust2toD --- n_fe_fromGust2toD
    CL_wb_fromGust2toD        = np.ones((N, n_Mass)) 
    alfa_fromGust2toD         = np.ones((N, n_Mass))  
    alfa_new_fromGust2toD     = np.ones((N, n_Mass))
    CD_fromGust2toD           = np.ones((N, n_Mass)) 
    q_fromGust2toD            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromGust2toD = np.ones((N, n_Mass))
    CM_due_to_CD_fromGust2toD = np.ones((N, n_Mass))
    CM_due_to_CT_fromGust2toD = np.ones((N, n_Mass))
    CM_CG_fromGust2toD        = np.ones((N, n_Mass))
    CL_ht_fromGust2toD        = np.ones((N, n_Mass))
    CL_new_fromGust2toD       = np.ones((N, n_Mass))
    L_wb_fromGust2toD         = np.ones((N, n_Mass))
    L_new_fromGust2toD        = np.ones((N, n_Mass))
    L_ht_fromGust2toD         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM GUST2 TO D
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromGust2toD[i,j], alfa_fromGust2toD[i,j], alfa_new_fromGust2toD[i,j],\
            CD_fromGust2toD[i,j], q_fromGust2toD[i,j], CM_due_to_CL_fromGust2toD[i,j],\
            CM_due_to_CD_fromGust2toD[i,j], CM_due_to_CT_fromGust2toD[i,j],\
            CM_CG_fromGust2toD[i,j], CL_ht_fromGust2toD[i,j], CL_new_fromGust2toD[i,j],\
            L_wb_fromGust2toD[i,j], L_new_fromGust2toD[i,j], L_ht_fromGust2toD[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromGust2toD[i,j], WS[j], n_fe_fromGust2toD[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromDto0     --- n_fe_fromDto0
    CL_wb_fromDto0        = np.ones((N, n_Mass)) 
    alfa_fromDto0         = np.ones((N, n_Mass)) 
    alfa_new_fromDto0     = np.ones((N, n_Mass))
    CD_fromDto0           = np.ones((N, n_Mass)) 
    q_fromDto0            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromDto0 = np.ones((N, n_Mass))
    CM_due_to_CD_fromDto0 = np.ones((N, n_Mass))
    CM_due_to_CT_fromDto0 = np.ones((N, n_Mass))
    CM_CG_fromDto0        = np.ones((N, n_Mass))
    CL_ht_fromDto0        = np.ones((N, n_Mass))
    CL_new_fromDto0       = np.ones((N, n_Mass))
    L_wb_fromDto0         = np.ones((N, n_Mass))
    L_new_fromDto0        = np.ones((N, n_Mass))
    L_ht_fromDto0         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM D TO 0
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromDto0[i,j], alfa_fromDto0[i,j], alfa_new_fromDto0[i,j],\
            CD_fromDto0[i,j], q_fromDto0[i,j], CM_due_to_CL_fromDto0[i,j],\
            CM_due_to_CD_fromDto0[i,j], CM_due_to_CT_fromDto0[i,j],\
            CM_CG_fromDto0[i,j], CL_ht_fromDto0[i,j], CL_new_fromDto0[i,j],\
            L_wb_fromDto0[i,j], L_new_fromDto0[i,j], L_ht_fromDto0[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromDto0[i,j], WS[j], n_fe_fromDto0[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
elif ( flight_env_case_pos == "Case2" ):
    
    flag1 = "straight"
    
    # V_fe_from0toS --- n_fe_from0toS
    CL_wb_from0toS        = np.ones((N, n_Mass)) 
    alfa_from0toS         = np.ones((N, n_Mass)) 
    alfa_new_from0toS     = np.ones((N, n_Mass)) 
    CD_from0toS           = np.ones((N, n_Mass)) 
    q_from0toS            = np.ones((N, n_Mass)) 
    CM_due_to_CL_from0toS = np.ones((N, n_Mass))
    CM_due_to_CD_from0toS = np.ones((N, n_Mass))
    CM_due_to_CT_from0toS = np.ones((N, n_Mass))
    CM_CG_from0toS        = np.ones((N, n_Mass))
    CL_ht_from0toS        = np.ones((N, n_Mass))
    CL_new_from0toS       = np.ones((N, n_Mass))
    L_wb_from0toS         = np.ones((N, n_Mass))
    L_new_from0toS        = np.ones((N, n_Mass))
    L_ht_from0toS         = np.ones((N, n_Mass))

    # BALANCING LOADS CALCULATIONS - FROM 0 TO S
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_from0toS[i,j], alfa_from0toS[i,j], alfa_new_from0toS[i,j],\
            CD_from0toS[i,j], q_from0toS[i,j], CM_due_to_CL_from0toS[i,j],\
            CM_due_to_CD_from0toS[i,j], CM_due_to_CT_from0toS[i,j],\
            CM_CG_from0toS[i,j], CL_ht_from0toS[i,j], CL_new_from0toS[i,j],\
            L_wb_from0toS[i,j], L_new_from0toS[i,j], L_ht_from0toS[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_from0toS[i,j], WS[j], n_fe_from0toS[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1)   
    
    # V_fe_fromStoA --- n_fe_fromStoA
    CL_wb_fromStoA        = np.ones((N, n_Mass)) 
    alfa_fromStoA         = np.ones((N, n_Mass)) 
    alfa_new_fromStoA     = np.ones((N, n_Mass)) 
    CD_fromStoA           = np.ones((N, n_Mass)) 
    q_fromStoA            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromStoA = np.ones((N, n_Mass))
    CM_due_to_CD_fromStoA = np.ones((N, n_Mass))
    CM_due_to_CT_fromStoA = np.ones((N, n_Mass))
    CM_CG_fromStoA        = np.ones((N, n_Mass))
    CL_ht_fromStoA        = np.ones((N, n_Mass))
    CL_new_fromStoA       = np.ones((N, n_Mass))
    L_wb_fromStoA         = np.ones((N, n_Mass))
    L_new_fromStoA        = np.ones((N, n_Mass))
    L_ht_fromStoA         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM S TO A
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromStoA[i,j], alfa_fromStoA[i,j], alfa_new_fromStoA[i,j],\
            CD_fromStoA[i,j], q_fromStoA[i,j], CM_due_to_CL_fromStoA[i,j],\
            CM_due_to_CD_fromStoA[i,j], CM_due_to_CT_fromStoA[i,j],\
            CM_CG_fromStoA[i,j], CL_ht_fromStoA[i,j], CL_new_fromStoA[i,j],\
            L_wb_fromStoA[i,j], L_new_fromStoA[i,j], L_ht_fromStoA[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromStoA[i,j], WS[j], n_fe_fromStoA[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1)    
    
    # V_fe_fromAtoC --- n_fe_fromAtoC
    CL_wb_fromAtoC        = np.ones((N, n_Mass)) 
    alfa_fromAtoC         = np.ones((N, n_Mass)) 
    alfa_new_fromAtoC     = np.ones((N, n_Mass)) 
    CD_fromAtoC           = np.ones((N, n_Mass)) 
    q_fromAtoC            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromAtoC = np.ones((N, n_Mass))
    CM_due_to_CD_fromAtoC = np.ones((N, n_Mass))
    CM_due_to_CT_fromAtoC = np.ones((N, n_Mass))
    CM_CG_fromAtoC        = np.ones((N, n_Mass))
    CL_ht_fromAtoC        = np.ones((N, n_Mass))
    CL_new_fromAtoC       = np.ones((N, n_Mass))
    L_wb_fromAtoC         = np.ones((N, n_Mass))
    L_new_fromAtoC        = np.ones((N, n_Mass))
    L_ht_fromAtoC         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM A TO C
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromAtoC[i,j], alfa_fromAtoC[i,j], alfa_new_fromAtoC[i,j],\
            CD_fromAtoC[i,j], q_fromAtoC[i,j], CM_due_to_CL_fromAtoC[i,j],\
            CM_due_to_CD_fromAtoC[i,j], CM_due_to_CT_fromAtoC[i,j],\
            CM_CG_fromAtoC[i,j], CL_ht_fromAtoC[i,j], CL_new_fromAtoC[i,j],\
            L_wb_fromAtoC[i,j], L_new_fromAtoC[i,j], L_ht_fromAtoC[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromAtoC[i,j], WS[j], n_fe_fromAtoC[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1)   
    
    # V_fe_fromCtoD --- n_fe_fromCtoD
    CL_wb_fromCtoD        = np.ones((N, n_Mass)) 
    alfa_fromCtoD         = np.ones((N, n_Mass)) 
    alfa_new_fromCtoD     = np.ones((N, n_Mass)) 
    CD_fromCtoD           = np.ones((N, n_Mass)) 
    q_fromCtoD            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromCtoD = np.ones((N, n_Mass))
    CM_due_to_CD_fromCtoD = np.ones((N, n_Mass))
    CM_due_to_CT_fromCtoD = np.ones((N, n_Mass))
    CM_CG_fromCtoD        = np.ones((N, n_Mass))
    CL_ht_fromCtoD        = np.ones((N, n_Mass))
    CL_new_fromCtoD       = np.ones((N, n_Mass))
    L_wb_fromCtoD         = np.ones((N, n_Mass))
    L_new_fromCtoD        = np.ones((N, n_Mass))
    L_ht_fromCtoD         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM C TO D
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromCtoD[i,j], alfa_fromCtoD[i,j], alfa_new_fromCtoD[i,j],\
            CD_fromCtoD[i,j], q_fromCtoD[i,j], CM_due_to_CL_fromCtoD[i,j],\
            CM_due_to_CD_fromCtoD[i,j], CM_due_to_CT_fromCtoD[i,j],\
            CM_CG_fromCtoD[i,j], CL_ht_fromCtoD[i,j], CL_new_fromCtoD[i,j],\
            L_wb_fromCtoD[i,j], L_new_fromCtoD[i,j], L_ht_fromCtoD[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromCtoD[i,j], WS[j], n_fe_fromCtoD[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1)   
    
    # V_fe_fromDto0 --- n_fe_fromDto0
    CL_wb_fromDto0        = np.ones((N, n_Mass)) 
    alfa_fromDto0         = np.ones((N, n_Mass)) 
    alfa_new_fromDto0     = np.ones((N, n_Mass)) 
    CD_fromDto0           = np.ones((N, n_Mass)) 
    q_fromDto0            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromDto0 = np.ones((N, n_Mass))
    CM_due_to_CD_fromDto0 = np.ones((N, n_Mass))
    CM_due_to_CT_fromDto0 = np.ones((N, n_Mass))
    CM_CG_fromDto0        = np.ones((N, n_Mass))
    CL_ht_fromDto0        = np.ones((N, n_Mass))
    CL_new_fromDto0       = np.ones((N, n_Mass))
    L_wb_fromDto0         = np.ones((N, n_Mass))
    L_new_fromDto0        = np.ones((N, n_Mass))
    L_ht_fromDto0         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM D TO 0
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromDto0[i,j], alfa_fromDto0[i,j], alfa_new_fromDto0[i,j],\
            CD_fromDto0[i,j], q_fromDto0[i,j], CM_due_to_CL_fromDto0[i,j],\
            CM_due_to_CD_fromDto0[i,j], CM_due_to_CT_fromDto0[i,j],\
            CM_CG_fromDto0[i,j], CL_ht_fromDto0[i,j], CL_new_fromDto0[i,j],\
            L_wb_fromDto0[i,j], L_new_fromDto0[i,j], L_ht_fromDto0[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromDto0[i,j], WS[j], n_fe_fromDto0[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
elif ( flight_env_case_pos == "Case3" ):
    
    flag1 = "straight"
    
    # V_fe_from0toS     --- n_fe_from0toS  
    CL_wb_from0toS        = np.ones((N, n_Mass)) 
    alfa_from0toS         = np.ones((N, n_Mass)) 
    alfa_new_from0toS     = np.ones((N, n_Mass)) 
    CD_from0toS           = np.ones((N, n_Mass)) 
    q_from0toS            = np.ones((N, n_Mass)) 
    CM_due_to_CL_from0toS = np.ones((N, n_Mass))
    CM_due_to_CD_from0toS = np.ones((N, n_Mass))
    CM_due_to_CT_from0toS = np.ones((N, n_Mass))
    CM_CG_from0toS        = np.ones((N, n_Mass))
    CL_ht_from0toS        = np.ones((N, n_Mass))
    CL_new_from0toS       = np.ones((N, n_Mass))
    L_wb_from0toS         = np.ones((N, n_Mass))
    L_new_from0toS        = np.ones((N, n_Mass))
    L_ht_from0toS         = np.ones((N, n_Mass))

    # BALANCING LOADS CALCULATIONS - FROM 0 TO S
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_from0toS[i,j], alfa_from0toS[i,j], alfa_new_from0toS[i,j],\
            CD_from0toS[i,j], q_from0toS[i,j], CM_due_to_CL_from0toS[i,j],\
            CM_due_to_CD_from0toS[i,j], CM_due_to_CT_from0toS[i,j],\
            CM_CG_from0toS[i,j], CL_ht_from0toS[i,j], CL_new_from0toS[i,j],\
            L_wb_from0toS[i,j], L_new_from0toS[i,j], L_ht_from0toS[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_from0toS[i,j], WS[j], n_fe_from0toS[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1)   
    
    # V_fe_fromStoA     --- n_fe_fromStoA   
    CL_wb_fromStoA        = np.ones((N, n_Mass)) 
    alfa_fromStoA         = np.ones((N, n_Mass)) 
    alfa_new_fromStoA     = np.ones((N, n_Mass)) 
    CD_fromStoA           = np.ones((N, n_Mass)) 
    q_fromStoA            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromStoA = np.ones((N, n_Mass))
    CM_due_to_CD_fromStoA = np.ones((N, n_Mass))
    CM_due_to_CT_fromStoA = np.ones((N, n_Mass))
    CM_CG_fromStoA        = np.ones((N, n_Mass))
    CL_ht_fromStoA        = np.ones((N, n_Mass))
    CL_new_fromStoA       = np.ones((N, n_Mass))
    L_wb_fromStoA         = np.ones((N, n_Mass))
    L_new_fromStoA        = np.ones((N, n_Mass))
    L_ht_fromStoA         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM S TO A
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromStoA[i,j], alfa_fromStoA[i,j], alfa_new_fromStoA[i,j],\
            CD_fromStoA[i,j], q_fromStoA[i,j], CM_due_to_CL_fromStoA[i,j],\
            CM_due_to_CD_fromStoA[i,j], CM_due_to_CT_fromStoA[i,j],\
            CM_CG_fromStoA[i,j], CL_ht_fromStoA[i,j], CL_new_fromStoA[i,j],\
            L_wb_fromStoA[i,j], L_new_fromStoA[i,j], L_ht_fromStoA[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromStoA[i,j], WS[j], n_fe_fromStoA[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1)    
    
    # V_fe_fromAtoGust1 --- n_fe_fromAtoGust1
    CL_wb_fromAtoGust1        = np.ones((N, n_Mass)) 
    alfa_fromAtoGust1         = np.ones((N, n_Mass)) 
    alfa_new_fromAtoGust1     = np.ones((N, n_Mass)) 
    CD_fromAtoGust1           = np.ones((N, n_Mass)) 
    q_fromAtoGust1            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromAtoGust1 = np.ones((N, n_Mass))
    CM_due_to_CD_fromAtoGust1 = np.ones((N, n_Mass))
    CM_due_to_CT_fromAtoGust1 = np.ones((N, n_Mass))
    CM_CG_fromAtoGust1        = np.ones((N, n_Mass))
    CL_ht_fromAtoGust1        = np.ones((N, n_Mass))
    CL_new_fromAtoGust1       = np.ones((N, n_Mass))
    L_wb_fromAtoGust1         = np.ones((N, n_Mass))
    L_new_fromAtoGust1        = np.ones((N, n_Mass))
    L_ht_fromAtoGust1         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM A TO GUST1
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromAtoGust1[i,j], alfa_fromAtoGust1[i,j], alfa_new_fromAtoGust1[i,j],\
            CD_fromAtoGust1[i,j], q_fromAtoGust1[i,j], CM_due_to_CL_fromAtoGust1[i,j],\
            CM_due_to_CD_fromAtoGust1[i,j], CM_due_to_CT_fromAtoGust1[i,j],\
            CM_CG_fromAtoGust1[i,j], CL_ht_fromAtoGust1[i,j], CL_new_fromAtoGust1[i,j],\
            L_wb_fromAtoGust1[i,j], L_new_fromAtoGust1[i,j], L_ht_fromAtoGust1[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromAtoGust1[i,j], WS[j], n_fe_fromAtoGust1[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromGust1toC --- n_fe_fromGust1toC
    CL_wb_fromGust1toC        = np.ones((N, n_Mass)) 
    alfa_fromGust1toC         = np.ones((N, n_Mass)) 
    alfa_new_fromGust1toC     = np.ones((N, n_Mass)) 
    CD_fromGust1toC           = np.ones((N, n_Mass)) 
    q_fromGust1toC            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromGust1toC = np.ones((N, n_Mass))
    CM_due_to_CD_fromGust1toC = np.ones((N, n_Mass))
    CM_due_to_CT_fromGust1toC = np.ones((N, n_Mass))
    CM_CG_fromGust1toC        = np.ones((N, n_Mass))
    CL_ht_fromGust1toC        = np.ones((N, n_Mass))
    CL_new_fromGust1toC       = np.ones((N, n_Mass))
    L_wb_fromGust1toC         = np.ones((N, n_Mass))
    L_new_fromGust1toC        = np.ones((N, n_Mass))
    L_ht_fromGust1toC         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM GUST1 TO C
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromGust1toC[i,j], alfa_fromGust1toC[i,j], alfa_new_fromGust1toC[i,j],\
            CD_fromGust1toC[i,j], q_fromGust1toC[i,j], CM_due_to_CL_fromGust1toC[i,j],\
            CM_due_to_CD_fromGust1toC[i,j], CM_due_to_CT_fromGust1toC[i,j],\
            CM_CG_fromGust1toC[i,j], CL_ht_fromGust1toC[i,j], CL_new_fromGust1toC[i,j],\
            L_wb_fromGust1toC[i,j], L_new_fromGust1toC[i,j], L_ht_fromGust1toC[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromGust1toC[i,j], WS[j], n_fe_fromGust1toC[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromCtoGust2 --- n_fe_fromCtoGust2
    CL_wb_fromCtoGust2        = np.ones((N, n_Mass)) 
    alfa_fromCtoGust2         = np.ones((N, n_Mass)) 
    alfa_new_fromCtoGust2     = np.ones((N, n_Mass)) 
    CD_fromCtoGust2           = np.ones((N, n_Mass)) 
    q_fromCtoGust2            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromCtoGust2 = np.ones((N, n_Mass))
    CM_due_to_CD_fromCtoGust2 = np.ones((N, n_Mass))
    CM_due_to_CT_fromCtoGust2 = np.ones((N, n_Mass))
    CM_CG_fromCtoGust2        = np.ones((N, n_Mass))
    CL_ht_fromCtoGust2        = np.ones((N, n_Mass))
    CL_new_fromCtoGust2       = np.ones((N, n_Mass))
    L_wb_fromCtoGust2         = np.ones((N, n_Mass))
    L_new_fromCtoGust2        = np.ones((N, n_Mass))
    L_ht_fromCtoGust2         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM C TO GUST2
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromCtoGust2[i,j], alfa_fromCtoGust2[i,j], alfa_new_fromCtoGust2[i,j],\
            CD_fromCtoGust2[i,j], q_fromCtoGust2[i,j], CM_due_to_CL_fromCtoGust2[i,j],\
            CM_due_to_CD_fromCtoGust2[i,j], CM_due_to_CT_fromCtoGust2[i,j],\
            CM_CG_fromCtoGust2[i,j], CL_ht_fromCtoGust2[i,j], CL_new_fromCtoGust2[i,j],\
            L_wb_fromCtoGust2[i,j], L_new_fromCtoGust2[i,j], L_ht_fromCtoGust2[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromCtoGust2[i,j], WS[j], n_fe_fromCtoGust2[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromGust2toD --- n_fe_fromGust2toD
    CL_wb_fromGust2toD        = np.ones((N, n_Mass)) 
    alfa_fromGust2toD         = np.ones((N, n_Mass)) 
    alfa_new_fromGust2toD     = np.ones((N, n_Mass)) 
    CD_fromGust2toD           = np.ones((N, n_Mass)) 
    q_fromGust2toD            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromGust2toD = np.ones((N, n_Mass))
    CM_due_to_CD_fromGust2toD = np.ones((N, n_Mass))
    CM_due_to_CT_fromGust2toD = np.ones((N, n_Mass))
    CM_CG_fromGust2toD        = np.ones((N, n_Mass))
    CL_ht_fromGust2toD        = np.ones((N, n_Mass))
    CL_new_fromGust2toD       = np.ones((N, n_Mass))
    L_wb_fromGust2toD         = np.ones((N, n_Mass))
    L_new_fromGust2toD        = np.ones((N, n_Mass))
    L_ht_fromGust2toD         = np.ones((N, n_Mass))

    # BALANCING LOADS CALCULATIONS - FROM GUST2 TO D
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromGust2toD[i,j], alfa_fromGust2toD[i,j], alfa_new_fromGust2toD[i,j], CD_fromGust2toD[i,j],\
            q_fromGust2toD[i,j], CM_due_to_CL_fromGust2toD[i,j],\
            CM_due_to_CD_fromGust2toD[i,j], CM_due_to_CT_fromGust2toD[i,j],\
            CM_CG_fromGust2toD[i,j], CL_ht_fromGust2toD[i,j], CL_new_fromGust2toD[i,j],\
            L_wb_fromGust2toD[i,j], L_new_fromGust2toD[i,j], L_ht_fromGust2toD[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromGust2toD[i,j], WS[j], n_fe_fromGust2toD[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromDto0     --- n_fe_fromDto0  
    CL_wb_fromDto0        = np.ones((N, n_Mass)) 
    alfa_fromDto0         = np.ones((N, n_Mass)) 
    alfa_new_fromDto0     = np.ones((N, n_Mass)) 
    CD_fromDto0           = np.ones((N, n_Mass)) 
    q_fromDto0            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromDto0 = np.ones((N, n_Mass))
    CM_due_to_CD_fromDto0 = np.ones((N, n_Mass))
    CM_due_to_CT_fromDto0 = np.ones((N, n_Mass))
    CM_CG_fromDto0        = np.ones((N, n_Mass))
    CL_ht_fromDto0        = np.ones((N, n_Mass))
    CL_new_fromDto0       = np.ones((N, n_Mass))
    L_wb_fromDto0         = np.ones((N, n_Mass))
    L_new_fromDto0        = np.ones((N, n_Mass))
    L_ht_fromDto0         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM D TO 0
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromDto0[i,j], alfa_fromDto0[i,j], alfa_new_fromDto0[i,j], CD_fromDto0[i,j],\
            q_fromDto0[i,j], CM_due_to_CL_fromDto0[i,j],\
            CM_due_to_CD_fromDto0[i,j], CM_due_to_CT_fromDto0[i,j],\
            CM_CG_fromDto0[i,j], CL_ht_fromDto0[i,j], CL_new_fromDto0[i,j],\
            L_wb_fromDto0[i,j], L_new_fromDto0[i,j], L_ht_fromDto0[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromDto0[i,j], WS[j], n_fe_fromDto0[i,j], CL_max_clean, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
###########################################
############# INVERTED FLIGHT #############   
########################################### 
if ( flight_env_case_neg == "Case1_inverted" ): 
    
    flag1 = "inverted"
    
    # V_fe_from0toSinv  --- n_fe_from0toSinv 
    CL_wb_from0toSinv        = np.ones((N, n_Mass)) 
    alfa_from0toSinv         = np.ones((N, n_Mass)) 
    alfa_new_from0toSinv     = np.ones((N, n_Mass)) 
    CD_from0toSinv           = np.ones((N, n_Mass)) 
    q_from0toSinv            = np.ones((N, n_Mass)) 
    CM_due_to_CL_from0toSinv = np.ones((N, n_Mass))
    CM_due_to_CD_from0toSinv = np.ones((N, n_Mass))
    CM_due_to_CT_from0toSinv = np.ones((N, n_Mass))
    CM_CG_from0toSinv        = np.ones((N, n_Mass))
    CL_ht_from0toSinv        = np.ones((N, n_Mass))
    CL_new_from0toSinv       = np.ones((N, n_Mass))
    L_wb_from0toSinv         = np.ones((N, n_Mass))
    L_new_from0toSinv        = np.ones((N, n_Mass))
    L_ht_from0toSinv         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO 0 TO S INVERTED 
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_from0toSinv[i,j], alfa_from0toSinv[i,j], alfa_new_from0toSinv[i,j],\
            CD_from0toSinv[i,j], q_from0toSinv[i,j], CM_due_to_CL_from0toSinv[i,j],\
            CM_due_to_CD_from0toSinv[i,j], CM_due_to_CT_from0toSinv[i,j],\
            CM_CG_from0toSinv[i,j], CL_ht_from0toSinv[i,j], CL_new_from0toSinv[i,j],\
            L_wb_from0toSinv[i,j], L_new_from0toSinv[i,j], L_ht_from0toSinv[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_from0toSinv[i,j], WS[j], n_fe_from0toSinv[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromSinvtoG  --- n_fe_fromSinvtoG
    CL_wb_fromSinvtoG        = np.ones((N, n_Mass)) 
    alfa_fromSinvtoG         = np.ones((N, n_Mass)) 
    alfa_new_fromSinvtoG     = np.ones((N, n_Mass)) 
    CD_fromSinvtoG           = np.ones((N, n_Mass)) 
    q_fromSinvtoG            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromSinvtoG = np.ones((N, n_Mass))
    CM_due_to_CD_fromSinvtoG = np.ones((N, n_Mass))
    CM_due_to_CT_fromSinvtoG = np.ones((N, n_Mass))
    CM_CG_fromSinvtoG        = np.ones((N, n_Mass))
    CL_ht_fromSinvtoG        = np.ones((N, n_Mass))
    CL_new_fromSinvtoG       = np.ones((N, n_Mass))
    L_wb_fromSinvtoG         = np.ones((N, n_Mass))
    L_new_fromSinvtoG        = np.ones((N, n_Mass))
    L_ht_fromSinvtoG         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO S INVERTED TO G
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromSinvtoG[i,j], alfa_fromSinvtoG[i,j], alfa_new_fromSinvtoG[i,j],\
            CD_fromSinvtoG[i,j], q_fromSinvtoG[i,j], CM_due_to_CL_fromSinvtoG[i,j],\
            CM_due_to_CD_fromSinvtoG[i,j], CM_due_to_CT_fromSinvtoG[i,j],\
            CM_CG_fromSinvtoG[i,j], CL_ht_fromSinvtoG[i,j], CL_new_fromSinvtoG[i,j],\
            L_wb_fromSinvtoG[i,j], L_new_fromSinvtoG[i,j], L_ht_fromSinvtoG[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromSinvtoG[i,j], WS[j], n_fe_fromSinvtoG[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromGtoGust1 --- n_fe_fromGtoGust1
    CL_wb_fromGtoGust1        = np.ones((N, n_Mass)) 
    alfa_fromGtoGust1         = np.ones((N, n_Mass)) 
    alfa_new_fromGtoGust1     = np.ones((N, n_Mass)) 
    CD_fromGtoGust1           = np.ones((N, n_Mass)) 
    q_fromGtoGust1            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromGtoGust1 = np.ones((N, n_Mass))
    CM_due_to_CD_fromGtoGust1 = np.ones((N, n_Mass))
    CM_due_to_CT_fromGtoGust1 = np.ones((N, n_Mass))
    CM_CG_fromGtoGust1        = np.ones((N, n_Mass))
    CL_ht_fromGtoGust1        = np.ones((N, n_Mass))
    CL_new_fromGtoGust1       = np.ones((N, n_Mass))
    L_wb_fromGtoGust1         = np.ones((N, n_Mass))
    L_new_fromGtoGust1        = np.ones((N, n_Mass))
    L_ht_fromGtoGust1         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM G TO GUST1
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromGtoGust1[i,j], alfa_fromGtoGust1[i,j], alfa_new_fromGtoGust1[i,j],\
            CD_fromGtoGust1[i,j], q_fromGtoGust1[i,j], CM_due_to_CL_fromGtoGust1[i,j],\
            CM_due_to_CD_fromGtoGust1[i,j], CM_due_to_CT_fromGtoGust1[i,j],\
            CM_CG_fromGtoGust1[i,j], CL_ht_fromGtoGust1[i,j], CL_new_fromGtoGust1[i,j],\
            L_wb_fromGtoGust1[i,j], L_new_fromGtoGust1[i,j], L_ht_fromGtoGust1[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromGtoGust1[i,j], WS[j], n_fe_fromGtoGust1[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromGust1toF --- n_fe_fromGust1toF
    CL_wb_fromGust1toF        = np.ones((N, n_Mass)) 
    alfa_fromGust1toF         = np.ones((N, n_Mass)) 
    alfa_new_fromGust1toF     = np.ones((N, n_Mass)) 
    CD_fromGust1toF           = np.ones((N, n_Mass)) 
    q_fromGust1toF            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromGust1toF = np.ones((N, n_Mass))
    CM_due_to_CD_fromGust1toF = np.ones((N, n_Mass))
    CM_due_to_CT_fromGust1toF = np.ones((N, n_Mass))
    CM_CG_fromGust1toF        = np.ones((N, n_Mass))
    CL_ht_fromGust1toF        = np.ones((N, n_Mass))
    CL_new_fromGust1toF       = np.ones((N, n_Mass))
    L_wb_fromGust1toF         = np.ones((N, n_Mass))
    L_new_fromGust1toF        = np.ones((N, n_Mass))
    L_ht_fromGust1toF         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROM GUST1 TO F
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromGust1toF[i,j], alfa_fromGust1toF[i,j], alfa_new_fromGust1toF[i,j],\
            CD_fromGust1toF[i,j], q_fromGust1toF[i,j], CM_due_to_CL_fromGust1toF[i,j],\
            CM_due_to_CD_fromGust1toF[i,j], CM_due_to_CT_fromGust1toF[i,j],\
            CM_CG_fromGust1toF[i,j], CL_ht_fromGust1toF[i,j], CL_new_fromGust1toF[i,j],\
            L_wb_fromGust1toF[i,j], L_new_fromGust1toF[i,j], L_ht_fromGust1toF[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromGust1toF[i,j], WS[j], n_fe_fromGust1toF[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromFtoE     --- n_fe_fromFtoE  
    CL_wb_fromFtoE        = np.ones((N, n_Mass)) 
    alfa_fromFtoE         = np.ones((N, n_Mass)) 
    alfa_new_fromFtoE     = np.ones((N, n_Mass)) 
    CD_fromFtoE           = np.ones((N, n_Mass)) 
    q_fromFtoE            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromFtoE = np.ones((N, n_Mass))
    CM_due_to_CD_fromFtoE = np.ones((N, n_Mass))
    CM_due_to_CT_fromFtoE = np.ones((N, n_Mass))
    CM_CG_fromFtoE        = np.ones((N, n_Mass))
    CL_ht_fromFtoE        = np.ones((N, n_Mass))
    CL_new_fromFtoE       = np.ones((N, n_Mass))
    L_wb_fromFtoE         = np.ones((N, n_Mass))
    L_new_fromFtoE        = np.ones((N, n_Mass))
    L_ht_fromFtoE         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO F TO E
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromFtoE[i,j], alfa_fromFtoE[i,j], alfa_new_fromFtoE[i,j],\
            CD_fromFtoE[i,j], q_fromFtoE[i,j], CM_due_to_CL_fromFtoE[i,j],\
            CM_due_to_CD_fromFtoE[i,j], CM_due_to_CT_fromFtoE[i,j],\
            CM_CG_fromFtoE[i,j], CL_ht_fromFtoE[i,j], CL_new_fromFtoE[i,j],\
            L_wb_fromFtoE[i,j], L_new_fromFtoE[i,j], L_ht_fromFtoE[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromFtoE[i,j], WS[j], n_fe_fromFtoE[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromEto0     --- n_fe_fromEto0    
    CL_wb_fromEto0        = np.ones((N, n_Mass)) 
    alfa_fromEto0         = np.ones((N, n_Mass)) 
    alfa_new_fromEto0     = np.ones((N, n_Mass)) 
    CD_fromEto0           = np.ones((N, n_Mass)) 
    q_fromEto0            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromEto0 = np.ones((N, n_Mass))
    CM_due_to_CD_fromEto0 = np.ones((N, n_Mass))
    CM_due_to_CT_fromEto0 = np.ones((N, n_Mass))
    CM_CG_fromEto0        = np.ones((N, n_Mass))
    CL_ht_fromEto0        = np.ones((N, n_Mass))
    CL_new_fromEto0       = np.ones((N, n_Mass))
    L_wb_fromEto0         = np.ones((N, n_Mass))
    L_new_fromEto0        = np.ones((N, n_Mass))
    L_ht_fromEto0         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO E TO 0
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromEto0[i,j], alfa_fromEto0[i,j], alfa_new_fromEto0[i,j],\
            CD_fromEto0[i,j], q_fromEto0[i,j], CM_due_to_CL_fromEto0[i,j],\
            CM_due_to_CD_fromEto0[i,j], CM_due_to_CT_fromEto0[i,j],\
            CM_CG_fromEto0[i,j], CL_ht_fromEto0[i,j], CL_new_fromEto0[i,j],\
            L_wb_fromEto0[i,j], L_new_fromEto0[i,j], L_ht_fromEto0[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromEto0[i,j], WS[j], n_fe_fromEto0[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
elif ( flight_env_case_neg == "Case2_inverted" ):
    
    flag1 = "inverted"
        
    # V_fe_from0toSinv --- n_fe_from0toSinv
    CL_wb_from0toSinv        = np.ones((N, n_Mass)) 
    alfa_from0toSinv         = np.ones((N, n_Mass)) 
    alfa_new_from0toSinv     = np.ones((N, n_Mass)) 
    CD_from0toSinv           = np.ones((N, n_Mass)) 
    q_from0toSinv            = np.ones((N, n_Mass)) 
    CM_due_to_CL_from0toSinv = np.ones((N, n_Mass))
    CM_due_to_CD_from0toSinv = np.ones((N, n_Mass))
    CM_due_to_CT_from0toSinv = np.ones((N, n_Mass))
    CM_CG_from0toSinv        = np.ones((N, n_Mass))
    CL_ht_from0toSinv        = np.ones((N, n_Mass))
    CL_new_from0toSinv       = np.ones((N, n_Mass))
    L_wb_from0toSinv         = np.ones((N, n_Mass))
    L_new_from0toSinv        = np.ones((N, n_Mass))
    L_ht_from0toSinv         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO 0 TO S INVERTED 
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_from0toSinv[i,j], alfa_from0toSinv[i,j], alfa_new_from0toSinv[i,j],\
            CD_from0toSinv[i,j], q_from0toSinv[i,j], CM_due_to_CL_from0toSinv[i,j],\
            CM_due_to_CD_from0toSinv[i,j], CM_due_to_CT_from0toSinv[i,j],\
            CM_CG_from0toSinv[i,j], CL_ht_from0toSinv[i,j], CL_new_from0toSinv[i,j],\
            L_wb_from0toSinv[i,j], L_new_from0toSinv[i,j], L_ht_from0toSinv[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_from0toSinv[i,j], WS[j], n_fe_from0toSinv[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromSinvtoG --- n_fe_fromSinvtoG
    CL_wb_fromSinvtoG        = np.ones((N, n_Mass)) 
    alfa_fromSinvtoG         = np.ones((N, n_Mass)) 
    alfa_new_fromSinvtoG     = np.ones((N, n_Mass)) 
    CD_fromSinvtoG           = np.ones((N, n_Mass)) 
    q_fromSinvtoG            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromSinvtoG = np.ones((N, n_Mass))
    CM_due_to_CD_fromSinvtoG = np.ones((N, n_Mass))
    CM_due_to_CT_fromSinvtoG = np.ones((N, n_Mass))
    CM_CG_fromSinvtoG        = np.ones((N, n_Mass))
    CL_ht_fromSinvtoG        = np.ones((N, n_Mass))
    CL_new_fromSinvtoG       = np.ones((N, n_Mass))
    L_wb_fromSinvtoG         = np.ones((N, n_Mass))
    L_new_fromSinvtoG        = np.ones((N, n_Mass))
    L_ht_fromSinvtoG         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO S INVERTED TO G
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromSinvtoG[i,j], alfa_fromSinvtoG[i,j], alfa_new_fromSinvtoG[i,j],\
            CD_fromSinvtoG[i,j], q_fromSinvtoG[i,j], CM_due_to_CL_fromSinvtoG[i,j],\
            CM_due_to_CD_fromSinvtoG[i,j], CM_due_to_CT_fromSinvtoG[i,j],\
            CM_CG_fromSinvtoG[i,j], CL_ht_fromSinvtoG[i,j], CL_new_fromSinvtoG[i,j],\
            L_wb_fromSinvtoG[i,j], L_new_fromSinvtoG[i,j], L_ht_fromSinvtoG[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromSinvtoG[i,j], WS[j], n_fe_fromSinvtoG[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromGtoF    --- n_fe_fromGtoF  
    CL_wb_fromGtoF        = np.ones((N, n_Mass)) 
    alfa_fromGtoF         = np.ones((N, n_Mass)) 
    alfa_new_fromGtoF     = np.ones((N, n_Mass)) 
    CD_fromGtoF           = np.ones((N, n_Mass)) 
    q_fromGtoF            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromGtoF = np.ones((N, n_Mass))
    CM_due_to_CD_fromGtoF = np.ones((N, n_Mass))
    CM_due_to_CT_fromGtoF = np.ones((N, n_Mass))
    CM_CG_fromGtoF        = np.ones((N, n_Mass))
    CL_ht_fromGtoF        = np.ones((N, n_Mass))
    CL_new_fromGtoF       = np.ones((N, n_Mass))
    L_wb_fromGtoF         = np.ones((N, n_Mass))
    L_new_fromGtoF        = np.ones((N, n_Mass))
    L_ht_fromGtoF         = np.ones((N, n_Mass)) 
    
    # BALANCING LOADS CALCULATIONS - FROMTO G TO F
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromGtoF[i,j], alfa_fromGtoF[i,j], alfa_new_fromGtoF[i,j],\
            CD_fromGtoF[i,j], q_fromGtoF[i,j], CM_due_to_CL_fromGtoF[i,j],\
            CM_due_to_CD_fromGtoF[i,j], CM_due_to_CT_fromGtoF[i,j],\
            CM_CG_fromGtoF[i,j], CL_ht_fromGtoF[i,j], CL_new_fromGtoF[i,j],\
            L_wb_fromGtoF[i,j], L_new_fromGtoF[i,j], L_ht_fromGtoF[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromGtoF[i,j], WS[j], n_fe_fromGtoF[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromFtoE    --- n_fe_fromFtoE  
    CL_wb_fromFtoE        = np.ones((N, n_Mass)) 
    alfa_fromFtoE         = np.ones((N, n_Mass)) 
    alfa_new_fromFtoE     = np.ones((N, n_Mass)) 
    CD_fromFtoE           = np.ones((N, n_Mass)) 
    q_fromFtoE            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromFtoE = np.ones((N, n_Mass))
    CM_due_to_CD_fromFtoE = np.ones((N, n_Mass))
    CM_due_to_CT_fromFtoE = np.ones((N, n_Mass))
    CM_CG_fromFtoE        = np.ones((N, n_Mass))
    CL_ht_fromFtoE        = np.ones((N, n_Mass))
    CL_new_fromFtoE       = np.ones((N, n_Mass))
    L_wb_fromFtoE         = np.ones((N, n_Mass))
    L_new_fromFtoE        = np.ones((N, n_Mass))
    L_ht_fromFtoE         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO F TO E
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromFtoE[i,j], alfa_fromFtoE[i,j], alfa_new_fromFtoE[i,j],\
            CD_fromFtoE[i,j], q_fromFtoE[i,j], CM_due_to_CL_fromFtoE[i,j],\
            CM_due_to_CD_fromFtoE[i,j], CM_due_to_CT_fromFtoE[i,j],\
            CM_CG_fromFtoE[i,j], CL_ht_fromFtoE[i,j], CL_new_fromFtoE[i,j],\
            L_wb_fromFtoE[i,j], L_new_fromFtoE[i,j], L_ht_fromFtoE[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromFtoE[i,j], WS[j], n_fe_fromFtoE[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromEto0    --- n_fe_fromEto0 
    CL_wb_fromEto0        = np.ones((N, n_Mass)) 
    alfa_fromEto0         = np.ones((N, n_Mass)) 
    alfa_new_fromEto0     = np.ones((N, n_Mass)) 
    CD_fromEto0           = np.ones((N, n_Mass)) 
    q_fromEto0            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromEto0 = np.ones((N, n_Mass))
    CM_due_to_CD_fromEto0 = np.ones((N, n_Mass))
    CM_due_to_CT_fromEto0 = np.ones((N, n_Mass))
    CM_CG_fromEto0        = np.ones((N, n_Mass))
    CL_ht_fromEto0        = np.ones((N, n_Mass))
    CL_new_fromEto0       = np.ones((N, n_Mass))
    L_wb_fromEto0         = np.ones((N, n_Mass))
    L_new_fromEto0        = np.ones((N, n_Mass))
    L_ht_fromEto0         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO E TO 0
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromEto0[i,j], alfa_fromEto0[i,j], alfa_new_fromEto0[i,j],\
            CD_fromEto0[i,j], q_fromEto0[i,j], CM_due_to_CL_fromEto0[i,j],\
            CM_due_to_CD_fromEto0[i,j], CM_due_to_CT_fromEto0[i,j],\
            CM_CG_fromEto0[i,j], CL_ht_fromEto0[i,j], CL_new_fromEto0[i,j],\
            L_wb_fromEto0[i,j], L_new_fromEto0[i,j], L_ht_fromEto0[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromEto0[i,j], WS[j], n_fe_fromEto0[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
elif ( flight_env_case_neg == "Case3_inverted" ):
    
    flag1 = "inverted"
        
    # V_fe_from0toSinv --- n_fe_from0toSinv
    CL_wb_from0toSinv        = np.ones((N, n_Mass)) 
    alfa_from0toSinv         = np.ones((N, n_Mass)) 
    alfa_new_from0toSinv     = np.ones((N, n_Mass)) 
    CD_from0toSinv           = np.ones((N, n_Mass)) 
    q_from0toSinv            = np.ones((N, n_Mass)) 
    CM_due_to_CL_from0toSinv = np.ones((N, n_Mass))
    CM_due_to_CD_from0toSinv = np.ones((N, n_Mass))
    CM_due_to_CT_from0toSinv = np.ones((N, n_Mass))
    CM_CG_from0toSinv        = np.ones((N, n_Mass))
    CL_ht_from0toSinv        = np.ones((N, n_Mass))
    CL_new_from0toSinv       = np.ones((N, n_Mass))
    L_wb_from0toSinv         = np.ones((N, n_Mass))
    L_new_from0toSinv        = np.ones((N, n_Mass))
    L_ht_from0toSinv         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO 0 TO S INVERTED 
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_from0toSinv[i,j], alfa_from0toSinv[i,j], alfa_new_from0toSinv[i,j],\
            CD_from0toSinv[i,j], q_from0toSinv[i,j], CM_due_to_CL_from0toSinv[i,j],\
            CM_due_to_CD_from0toSinv[i,j], CM_due_to_CT_from0toSinv[i,j],\
            CM_CG_from0toSinv[i,j], CL_ht_from0toSinv[i,j], CL_new_from0toSinv[i,j],\
            L_wb_from0toSinv[i,j], L_new_from0toSinv[i,j], L_ht_from0toSinv[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_from0toSinv[i,j], WS[j], n_fe_from0toSinv[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromSinvtoG --- n_fe_fromSinvtoG
    CL_wb_fromSinvtoG        = np.ones((N, n_Mass)) 
    alfa_fromSinvtoG         = np.ones((N, n_Mass)) 
    alfa_new_fromSinvtoG     = np.ones((N, n_Mass)) 
    CD_fromSinvtoG           = np.ones((N, n_Mass)) 
    q_fromSinvtoG            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromSinvtoG = np.ones((N, n_Mass))
    CM_due_to_CD_fromSinvtoG = np.ones((N, n_Mass))
    CM_due_to_CT_fromSinvtoG = np.ones((N, n_Mass))
    CM_CG_fromSinvtoG        = np.ones((N, n_Mass))
    CL_ht_fromSinvtoG        = np.ones((N, n_Mass))
    CL_new_fromSinvtoG       = np.ones((N, n_Mass))
    L_wb_fromSinvtoG         = np.ones((N, n_Mass))
    L_new_fromSinvtoG        = np.ones((N, n_Mass))
    L_ht_fromSinvtoG         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO S INVERTED TO G
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromSinvtoG[i,j], alfa_fromSinvtoG[i,j], alfa_new_fromSinvtoG[i,j],\
            CD_fromSinvtoG[i,j], q_fromSinvtoG[i,j], CM_due_to_CL_fromSinvtoG[i,j],\
            CM_due_to_CD_fromSinvtoG[i,j], CM_due_to_CT_fromSinvtoG[i,j],\
            CM_CG_fromSinvtoG[i,j], CL_ht_fromSinvtoG[i,j], CL_new_fromSinvtoG[i,j],\
            L_wb_fromSinvtoG[i,j], L_new_fromSinvtoG[i,j], L_ht_fromSinvtoG[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromSinvtoG[i,j], WS[j], n_fe_fromSinvtoG[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromGtoF    --- n_fe_fromGtoF  
    CL_wb_fromGtoF        = np.ones((N, n_Mass)) 
    alfa_fromGtoF         = np.ones((N, n_Mass)) 
    alfa_new_fromGtoF     = np.ones((N, n_Mass)) 
    CD_fromGtoF           = np.ones((N, n_Mass)) 
    q_fromGtoF            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromGtoF = np.ones((N, n_Mass))
    CM_due_to_CD_fromGtoF = np.ones((N, n_Mass))
    CM_due_to_CT_fromGtoF = np.ones((N, n_Mass))
    CM_CG_fromGtoF        = np.ones((N, n_Mass))
    CL_ht_fromGtoF        = np.ones((N, n_Mass))
    CL_new_fromGtoF       = np.ones((N, n_Mass))
    L_wb_fromGtoF         = np.ones((N, n_Mass))
    L_new_fromGtoF        = np.ones((N, n_Mass))
    L_ht_fromGtoF         = np.ones((N, n_Mass)) 
    
    # BALANCING LOADS CALCULATIONS - FROMTO G TO F
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromGtoF[i,j], alfa_fromGtoF[i,j], alfa_new_fromGtoF[i,j],\
            CD_fromGtoF[i,j], q_fromGtoF[i,j], CM_due_to_CL_fromGtoF[i,j],\
            CM_due_to_CD_fromGtoF[i,j], CM_due_to_CT_fromGtoF[i,j],\
            CM_CG_fromGtoF[i,j], CL_ht_fromGtoF[i,j], CL_new_fromGtoF[i,j],\
            L_wb_fromGtoF[i,j], L_new_fromGtoF[i,j], L_ht_fromGtoF[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromGtoF[i,j], WS[j], n_fe_fromGtoF[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromFtoE    --- n_fe_fromFtoE   
    CL_wb_fromFtoE        = np.ones((N, n_Mass)) 
    alfa_fromFtoE         = np.ones((N, n_Mass)) 
    alfa_new_fromFtoE     = np.ones((N, n_Mass)) 
    CD_fromFtoE           = np.ones((N, n_Mass)) 
    q_fromFtoE            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromFtoE = np.ones((N, n_Mass))
    CM_due_to_CD_fromFtoE = np.ones((N, n_Mass))
    CM_due_to_CT_fromFtoE = np.ones((N, n_Mass))
    CM_CG_fromFtoE        = np.ones((N, n_Mass))
    CL_ht_fromFtoE        = np.ones((N, n_Mass))
    CL_new_fromFtoE       = np.ones((N, n_Mass))
    L_wb_fromFtoE         = np.ones((N, n_Mass))
    L_new_fromFtoE        = np.ones((N, n_Mass))
    L_ht_fromFtoE         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO F TO E
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromFtoE[i,j], alfa_fromFtoE[i,j], alfa_new_fromFtoE[i,j],\
            CD_fromFtoE[i,j], q_fromFtoE[i,j], CM_due_to_CL_fromFtoE[i,j],\
            CM_due_to_CD_fromFtoE[i,j], CM_due_to_CT_fromFtoE[i,j],\
            CM_CG_fromFtoE[i,j], CL_ht_fromFtoE[i,j], CL_new_fromFtoE[i,j],\
            L_wb_fromFtoE[i,j], L_new_fromFtoE[i,j], L_ht_fromFtoE[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromFtoE[i,j], WS[j], n_fe_fromFtoE[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
    
    # V_fe_fromEto0    --- n_fe_fromEto0   
    CL_wb_fromEto0        = np.ones((N, n_Mass)) 
    alfa_fromEto0         = np.ones((N, n_Mass)) 
    alfa_new_fromEto0     = np.ones((N, n_Mass)) 
    CD_fromEto0           = np.ones((N, n_Mass)) 
    q_fromEto0            = np.ones((N, n_Mass)) 
    CM_due_to_CL_fromEto0 = np.ones((N, n_Mass))
    CM_due_to_CD_fromEto0 = np.ones((N, n_Mass))
    CM_due_to_CT_fromEto0 = np.ones((N, n_Mass))
    CM_CG_fromEto0        = np.ones((N, n_Mass))
    CL_ht_fromEto0        = np.ones((N, n_Mass))
    CL_new_fromEto0       = np.ones((N, n_Mass))
    L_wb_fromEto0         = np.ones((N, n_Mass))
    L_new_fromEto0        = np.ones((N, n_Mass))
    L_ht_fromEto0         = np.ones((N, n_Mass))
    
    # BALANCING LOADS CALCULATIONS - FROMTO F TO E
    for j in range(0, n_Mass):
        for i in range(0, N):
            CL_wb_fromEto0[i,j], alfa_fromEto0[i,j], alfa_new_fromEto0[i,j],\
            CD_fromEto0[i,j], q_fromEto0[i,j], CM_due_to_CL_fromEto0[i,j],\
            CM_due_to_CD_fromEto0[i,j], CM_due_to_CT_fromEto0[i,j],\
            CM_CG_fromEto0[i,j], CL_ht_fromEto0[i,j], CL_new_fromEto0[i,j],\
            L_wb_fromEto0[i,j], L_new_fromEto0[i,j], L_ht_fromEto0[i,j] = obj1.flight_balancing_parameters(rho0, V_fe_fromEto0[i,j], WS[j], n_fe_fromEto0[i,j], CL_inverted, alfa0L,\
                                            p_cd_wb, CL_star, p_cl_wb1,\
                                            p_cl_wb2, x_ac_envelope[j], xcg_envelope[j],\
                                            bcg_envelope[j], MAC, thrust_line_envelope[j],\
                                            CM0, CM_gear, S_wing, l_ht_envelope[j],\
                                            flag1, flag2, p_CM_wb[1], obj1) 
                
###############################################################################
# PLOTTING RESULTS - MAIN WING LIFT, BOTH WING BODY AND FULLVEHICLE
###############################################################################
n1 = 1   
n2 = 1         
for i in range(0, n_Mass):    
    figure_name1 = "fig" + str(n1)  
    n1 = n1 + 1 
    figure_name2 = 'fig' + str(n1)    
    n1 = n1 + 1
    figure_name3 = 'fig' + str(n1)             
    figure_name1, figure_name2, figure_name3 = obj1.main_wing_lift(V_fe_from0toS[:,i], V_fe_fromStoA[:,i], V_fe_fromAtoGust1[:,i],\
                       V_fe_fromGust1toC[:,i], V_fe_fromCtoGust2[:,i], V_fe_fromGust2toD[:,i],\
                       V_fe_fromDto0[:,i], V_fe_from0toSinv[:,i], V_fe_fromSinvtoG[:,i],\
                       V_fe_fromGtoGust1[:,i], V_fe_fromGust1toF[:,i], V_fe_fromFtoE[:,i],\
                       V_fe_fromEto0[:,i], V_fe_fromGtoF[:,i], V_fe_fromAtoC[:,i],\
                       V_fe_fromCtoD[:,i], L_wb_from0toS[:,i], L_new_from0toS[:,i], L_wb_fromStoA[:,i],\
                       L_new_fromStoA[:,i], L_new_fromAtoGust1[:,i], L_wb_fromAtoGust1[:,i],\
                       L_new_fromGust1toC[:,i], L_wb_fromGust1toC[:,i], L_new_fromCtoGust2[:,i],\
                       L_wb_fromCtoGust2[:,i], L_new_fromGust2toD[:,i], L_wb_fromGust2toD[:,i],\
                       L_new_fromDto0[:,i], L_wb_fromDto0[:,i], L_new_from0toSinv[:,i],\
                       L_wb_from0toSinv[:,i], L_new_fromSinvtoG[:,i], L_wb_fromSinvtoG[:,i],\
                       L_new_fromGtoGust1[:,i], L_wb_fromGtoGust1[:,i], L_new_fromGust1toF[:,i],\
                       L_wb_fromGust1toF[:,i], L_new_fromFtoE[:,i], L_wb_fromFtoE[:,i],\
                       L_new_fromEto0[:,i], L_wb_fromEto0[:,i], L_new_fromGtoF[:,i],\
                       L_wb_fromGtoF[:,i], L_new_fromAtoC[:,i], L_wb_fromAtoC[:,i],\
                       L_new_fromCtoD[:,i], L_wb_fromCtoD[:,i], L_ht_from0toS[:,i],\
                       L_ht_fromStoA[:,i], L_ht_fromAtoGust1[:,i], L_ht_fromGust1toC[:,i],\
                       L_ht_fromCtoGust2[:,i], L_ht_fromGust2toD[:,i], L_ht_fromDto0[:,i],\
                       L_ht_from0toSinv[:,i], L_ht_fromSinvtoG[:,i], L_ht_fromGtoGust1[:,i],\
                       L_ht_fromGust1toF[:,i], L_ht_fromFtoE[:,i],  L_ht_fromEto0[:,i],\
                       L_ht_fromGtoF[:,i], L_ht_fromAtoC[:,i], L_ht_fromCtoD[:,i], L_ht_unit_load_factor[:,i],\
                       V_unit_load_factor[:,i], Aircraft_name, n2, flight_env_case_pos, flight_env_case_neg)   
    n1 = n1 + 1    
    n2 = n2 + 1
    
# =============================================================================    
# MOVING FIGURES INSIDE OUTPUT FOLDER
# FinalEnvelope
n = 1
for i in range(0, n_Mass):
    
    figure_name = r'\MainwingLiftFullVehicle' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\BalancingLoads'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)
    
    figure_name = r'\MainwingLiftWingBody' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\BalancingLoads'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)
    
    figure_name = r'\HorizontalTailplaneLift' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\BalancingLoads'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)

    n = n + 1
# ============================================================================= 
###############################   
###### CLOSE ALL FIGURES ######
###############################
plt.close("all")    
###############################################################################
################## STORE INSIDE THE SIMPLE NAMESPACE OBJECT ################### 
############################################################################### 
Balancin_loads  = {}
Balancing_loads = obj1.balancing_values_storage(   
                       Balancin_loads,\
                       CM_CG_from0toS, CM_CG_fromStoA, CM_CG_fromAtoGust1,\
                       CM_CG_fromGust1toC, CM_CG_fromCtoGust2, CM_CG_fromGust2toD,\
                       CM_CG_fromDto0, CM_CG_from0toSinv, CM_CG_fromSinvtoG,\
                       CM_CG_fromGtoGust1, CM_CG_fromGust1toF, CM_CG_fromFtoE,\
                       CM_CG_fromEto0, CM_CG_fromGtoF, CM_CG_fromAtoC, CM_CG_fromCtoD,\
                           
                       CM_due_to_CD_from0toS, CM_due_to_CD_fromStoA, CM_due_to_CD_fromAtoGust1,\
                       CM_due_to_CD_fromGust1toC, CM_due_to_CD_fromCtoGust2, CM_due_to_CD_fromGust2toD,\
                       CM_due_to_CD_fromDto0, CM_due_to_CD_from0toSinv, CM_due_to_CD_fromSinvtoG,\
                       CM_due_to_CD_fromGtoGust1, CM_due_to_CD_fromGust1toF, CM_due_to_CD_fromFtoE,\
                       CM_due_to_CD_fromEto0, CM_due_to_CD_fromGtoF, CM_due_to_CD_fromAtoC, CM_due_to_CD_fromCtoD,\
                           
                       CM_due_to_CL_from0toS, CM_due_to_CL_fromStoA, CM_due_to_CL_fromAtoGust1,\
                       CM_due_to_CL_fromGust1toC, CM_due_to_CL_fromCtoGust2, CM_due_to_CL_fromGust2toD,\
                       CM_due_to_CL_fromDto0, CM_due_to_CL_from0toSinv, CM_due_to_CL_fromSinvtoG,\
                       CM_due_to_CL_fromGtoGust1, CM_due_to_CL_fromGust1toF, CM_due_to_CL_fromFtoE,\
                       CM_due_to_CL_fromEto0, CM_due_to_CL_fromGtoF, CM_due_to_CL_fromAtoC, CM_due_to_CL_fromCtoD,\
                           
                       q_from0toS, q_fromStoA, q_fromAtoGust1,\
                       q_fromGust1toC, q_fromCtoGust2, q_fromGust2toD,\
                       q_fromDto0, q_from0toSinv, q_fromSinvtoG,\
                       q_fromGtoGust1, q_fromGust1toF, q_fromFtoE,\
                       q_fromEto0, q_fromGtoF, q_fromAtoC, q_fromCtoD,\
                           
                       CD_from0toS, CD_fromStoA, CD_fromAtoGust1,\
                       CD_fromGust1toC, CD_fromCtoGust2, CD_fromGust2toD,\
                       CD_fromDto0, CD_from0toSinv, CD_fromSinvtoG,\
                       CD_fromGtoGust1, CD_fromGust1toF, CD_fromFtoE,\
                       CD_fromEto0, CD_fromGtoF, CD_fromAtoC, CD_fromCtoD,\
                           
                       alfa_new_from0toS, alfa_new_fromStoA, alfa_new_fromAtoGust1,\
                       alfa_new_fromGust1toC, alfa_new_fromCtoGust2, alfa_new_fromGust2toD,\
                       alfa_new_fromDto0, alfa_new_from0toSinv, alfa_new_fromSinvtoG,\
                       alfa_new_fromGtoGust1, alfa_new_fromGust1toF, alfa_new_fromFtoE,\
                       alfa_new_fromEto0, alfa_new_fromGtoF, alfa_new_fromAtoC, alfa_new_fromCtoD,\
                           
                       alfa_from0toS, alfa_fromStoA, alfa_fromAtoGust1,\
                       alfa_fromGust1toC, alfa_fromCtoGust2, alfa_fromGust2toD,\
                       alfa_fromDto0, alfa_from0toSinv, alfa_fromSinvtoG,\
                       alfa_fromGtoGust1, alfa_fromGust1toF, alfa_fromFtoE,\
                       alfa_fromEto0, alfa_fromGtoF, alfa_fromAtoC, alfa_fromCtoD,\
                      
                       CL_wb_from0toS, CL_wb_fromStoA, CL_wb_fromAtoGust1,\
                       CL_wb_fromGust1toC, CL_wb_fromCtoGust2, CL_wb_fromGust2toD,\
                       CL_wb_fromDto0, CL_wb_from0toSinv, CL_wb_fromSinvtoG,\
                       CL_wb_fromGtoGust1, CL_wb_fromGust1toF, CL_wb_fromFtoE,\
                       CL_wb_fromEto0, CL_wb_fromGtoF, CL_wb_fromAtoC, CL_wb_fromCtoD,\
                      
                       CL_ht_from0toS, CL_ht_fromStoA, CL_ht_fromAtoGust1,\
                       CL_ht_fromGust1toC, CL_ht_fromCtoGust2, CL_ht_fromGust2toD,\
                       CL_ht_fromDto0, CL_ht_from0toSinv, CL_ht_fromSinvtoG,\
                       CL_ht_fromGtoGust1, CL_ht_fromGust1toF, CL_ht_fromFtoE,\
                       CL_ht_fromEto0, CL_ht_fromGtoF, CL_ht_fromAtoC, CL_ht_fromCtoD,\
                      
                       CL_new_from0toS, CL_new_fromStoA, CL_new_fromAtoGust1,\
                       CL_new_fromGust1toC, CL_new_fromCtoGust2, CL_new_fromGust2toD,\
                       CL_new_fromDto0, CL_new_from0toSinv, CL_new_fromSinvtoG,\
                       CL_new_fromGtoGust1, CL_new_fromGust1toF, CL_new_fromFtoE,\
                       CL_new_fromEto0, CL_new_fromGtoF, CL_new_fromAtoC, CL_new_fromCtoD,\
                           
                       n_fe_from0toS, n_fe_fromStoA, n_fe_fromAtoGust1,\
                       n_fe_fromGust1toC, n_fe_fromCtoGust2, n_fe_fromGust2toD,\
                       n_fe_fromDto0, n_fe_from0toSinv, n_fe_fromSinvtoG,\
                       n_fe_fromGtoGust1, n_fe_fromGust1toF, n_fe_fromFtoE,\
                       n_fe_fromEto0, n_fe_fromGtoF, n_fe_fromAtoC, n_fe_fromCtoD,\
                           
                       V_fe_from0toS, V_fe_fromStoA, V_fe_fromAtoGust1,\
                       V_fe_fromGust1toC, V_fe_fromCtoGust2, V_fe_fromGust2toD,\
                       V_fe_fromDto0, V_fe_from0toSinv, V_fe_fromSinvtoG,\
                       V_fe_fromGtoGust1, V_fe_fromGust1toF, V_fe_fromFtoE,\
                       V_fe_fromEto0, V_fe_fromGtoF, V_fe_fromAtoC, V_fe_fromCtoD,\
                           
                       L_wb_from0toS, L_new_from0toS, L_wb_fromStoA,\
                       L_new_fromStoA, L_new_fromAtoGust1, L_wb_fromAtoGust1,\
                       L_new_fromGust1toC, L_wb_fromGust1toC, L_new_fromCtoGust2,\
                       L_wb_fromCtoGust2, L_new_fromGust2toD, L_wb_fromGust2toD,\
                       L_new_fromDto0, L_wb_fromDto0, L_new_from0toSinv,\
                       L_wb_from0toSinv, L_new_fromSinvtoG, L_wb_fromSinvtoG,\
                       L_new_fromGtoGust1, L_wb_fromGtoGust1, L_new_fromGust1toF,\
                       L_wb_fromGust1toF, L_new_fromFtoE, L_wb_fromFtoE,\
                       L_new_fromEto0, L_wb_fromEto0, L_new_fromGtoF,\
                       L_wb_fromGtoF, L_new_fromAtoC, L_wb_fromAtoC,\
                       L_new_fromCtoD, L_wb_fromCtoD, L_ht_from0toS,\
                       L_ht_fromStoA, L_ht_fromAtoGust1, L_ht_fromGust1toC,\
                       L_ht_fromCtoGust2, L_ht_fromGust2toD, L_ht_fromDto0,\
                       L_ht_from0toSinv, L_ht_fromSinvtoG, L_ht_fromGtoGust1,\
                       L_ht_fromGust1toF, L_ht_fromFtoE,  L_ht_fromEto0,\
                       L_ht_fromGtoF, L_ht_fromAtoC, L_ht_fromCtoD,\
                       flight_env_case_pos, flight_env_case_neg)
# UPDATING DICTIONARY
aircraft.update(Balancing_loads)
 # UPDATING THE SIMPLE NAMESPACE OBJECT
aircraft_data = SimpleNamespace(**aircraft)      
         