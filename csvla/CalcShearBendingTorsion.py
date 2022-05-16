# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:22:59 2022

@author: claum
"""

import numpy as np
import ShearBendingTorsionClass
obj2 = ShearBendingTorsionClass.ShearBendingTorsionClass

# INITIALIZATION 
n_Station   = len(y_halfspan)
twist_angle = aircraft_data.Wing['first_twist_angle']['Value']

# Preallocation - Point S
alfaS         = np.ones((n_Mass, 1))
qS            = np.ones((n_Mass, 1))
cCl_S         = np.ones((n_Station, n_Mass))
cCd_S         = np.ones((n_Station, n_Mass))
AoA_tot_deg_S = np.ones((n_Mass, 1))
AoA_tot_rad_S = np.ones((n_Mass, 1))
cCz_S         = np.ones((n_Station, n_Mass))
cCa_S         = np.ones((n_Station, n_Mass))
Norm_force_S  = np.ones((n_Station, n_Mass))
Axial_force_S = np.ones((n_Station, n_Mass))
Shear_S       = np.ones((n_Station, n_Mass))
Bending_S     = np.ones((n_Station, n_Mass))
m_distr_S     = np.ones((n_Station, n_Mass))
Torsion_S     = np.ones((n_Station, n_Mass))

# Preallocation - Point A
alfaA         = np.ones((n_Mass, 1))
qA            = np.ones((n_Mass, 1))
cCl_A         = np.ones((n_Station, n_Mass))
cCd_A         = np.ones((n_Station, n_Mass))
AoA_tot_deg_A = np.ones((n_Mass, 1))
AoA_tot_rad_A = np.ones((n_Mass, 1))
cCz_A         = np.ones((n_Station, n_Mass))
cCa_A         = np.ones((n_Station, n_Mass))
Norm_force_A  = np.ones((n_Station, n_Mass))
Axial_force_A = np.ones((n_Station, n_Mass))
Shear_A       = np.ones((n_Station, n_Mass))
Bending_A     = np.ones((n_Station, n_Mass))
m_distr_A     = np.ones((n_Station, n_Mass))
Torsion_A     = np.ones((n_Station, n_Mass))

# Preallocation - Point C
alfaC         = np.ones((n_Mass, 1))
qC            = np.ones((n_Mass, 1))
cCl_C         = np.ones((n_Station, n_Mass))
cCd_C         = np.ones((n_Station, n_Mass))
AoA_tot_deg_C = np.ones((n_Mass, 1))
AoA_tot_rad_C = np.ones((n_Mass, 1))
cCz_C         = np.ones((n_Station, n_Mass))
cCa_C         = np.ones((n_Station, n_Mass))
Norm_force_C  = np.ones((n_Station, n_Mass))
Axial_force_C = np.ones((n_Station, n_Mass))
Shear_C       = np.ones((n_Station, n_Mass))
Bending_C     = np.ones((n_Station, n_Mass))
m_distr_C     = np.ones((n_Station, n_Mass))
Torsion_C     = np.ones((n_Station, n_Mass))

# Preallocation - Point D
alfaD         = np.ones((n_Mass, 1))
qD            = np.ones((n_Mass, 1))
cCl_D         = np.ones((n_Station, n_Mass))
cCd_D         = np.ones((n_Station, n_Mass))
AoA_tot_deg_D = np.ones((n_Mass, 1))
AoA_tot_rad_D = np.ones((n_Mass, 1))
cCz_D         = np.ones((n_Station, n_Mass))
cCa_D         = np.ones((n_Station, n_Mass))
Norm_force_D  = np.ones((n_Station, n_Mass))
Axial_force_D = np.ones((n_Station, n_Mass))
Shear_D       = np.ones((n_Station, n_Mass))
Bending_D     = np.ones((n_Station, n_Mass))
m_distr_D     = np.ones((n_Station, n_Mass))
Torsion_D     = np.ones((n_Station, n_Mass))

# Preallocation - Point S inverted
alfaSinv         = np.ones((n_Mass, 1))
qSinv            = np.ones((n_Mass, 1))
cCl_Sinv         = np.ones((n_Station, n_Mass))
cCd_Sinv         = np.ones((n_Station, n_Mass))
AoA_tot_deg_Sinv = np.ones((n_Mass, 1))
AoA_tot_rad_Sinv = np.ones((n_Mass, 1))
cCz_Sinv         = np.ones((n_Station, n_Mass))
cCa_Sinv         = np.ones((n_Station, n_Mass))
Norm_force_Sinv  = np.ones((n_Station, n_Mass))
Axial_force_Sinv = np.ones((n_Station, n_Mass))
Shear_Sinv       = np.ones((n_Station, n_Mass))
Bending_Sinv     = np.ones((n_Station, n_Mass))
m_distr_Sinv     = np.ones((n_Station, n_Mass))
Torsion_Sinv     = np.ones((n_Station, n_Mass))

# Preallocation - Point G
alfaG         = np.ones((n_Mass, 1))
qG            = np.ones((n_Mass, 1))
cCl_G         = np.ones((n_Station, n_Mass))
cCd_G         = np.ones((n_Station, n_Mass))
AoA_tot_deg_G = np.ones((n_Mass, 1))
AoA_tot_rad_G = np.ones((n_Mass, 1))
cCz_G         = np.ones((n_Station, n_Mass))
cCa_G         = np.ones((n_Station, n_Mass))
Norm_force_G  = np.ones((n_Station, n_Mass))
Axial_force_G = np.ones((n_Station, n_Mass))
Shear_G       = np.ones((n_Station, n_Mass))
Bending_G     = np.ones((n_Station, n_Mass))
m_distr_G     = np.ones((n_Station, n_Mass))
Torsion_G     = np.ones((n_Station, n_Mass))

# Preallocation - Point F
alfaF         = np.ones((n_Mass, 1))
qF            = np.ones((n_Mass, 1))
cCl_F         = np.ones((n_Station, n_Mass))
cCd_F         = np.ones((n_Station, n_Mass))
AoA_tot_deg_F = np.ones((n_Mass, 1))
AoA_tot_rad_F = np.ones((n_Mass, 1))
cCz_F         = np.ones((n_Station, n_Mass))
cCa_F         = np.ones((n_Station, n_Mass))
Norm_force_F  = np.ones((n_Station, n_Mass))
Axial_force_F = np.ones((n_Station, n_Mass))
Shear_F       = np.ones((n_Station, n_Mass))
Bending_F     = np.ones((n_Station, n_Mass))
m_distr_F     = np.ones((n_Station, n_Mass))
Torsion_F     = np.ones((n_Station, n_Mass))

# Preallocation - Point E
alfaE         = np.ones((n_Mass, 1))
qE            = np.ones((n_Mass, 1))
cCl_E         = np.ones((n_Station, n_Mass))
cCd_E         = np.ones((n_Station, n_Mass))
AoA_tot_deg_E = np.ones((n_Mass, 1))
AoA_tot_rad_E = np.ones((n_Mass, 1))
cCz_E         = np.ones((n_Station, n_Mass))
cCa_E         = np.ones((n_Station, n_Mass))
Norm_force_E  = np.ones((n_Station, n_Mass))
Axial_force_E = np.ones((n_Station, n_Mass))
Shear_E       = np.ones((n_Station, n_Mass))
Bending_E     = np.ones((n_Station, n_Mass))
m_distr_E     = np.ones((n_Station, n_Mass))
Torsion_E     = np.ones((n_Station, n_Mass))

###########################################################################
################### STRAIGHT AND INVERTED FLIGHT POINTS ###################
###########################################################################

# Figure number - n_fig
n_fig = 1

for j in range(0, n_Mass): 
    
    #######################
    ####### POINT S #######
    #######################
    alfaS[j] = alfa_from0toS[-1,j]
    qS[j]    = q_from0toS[-1,j]
    Point_S  = "PointS"
    cCl_S[:,j], cCd_S[:,j], AoA_tot_deg_S[j], AoA_tot_rad_S[j], cCz_S[:,j],\
    cCa_S[:,j], Norm_force_S[:,j], Axial_force_S[:,j], Shear_S[:,j],\
    Bending_S[:,j], m_distr_S[:,j], Torsion_S[:,j],\
    fig_S = obj2.CalcShearBendingTorsion_func(chord[:], y_halfspan[:],\
                                                     cl_S[:,j], cd_S[:,j],\
                                                     cm_S[:,j], alfaS[j],\
                                                     twist_angle, qS[j],\
                                                     Point_S, n_fig, obj2)
    #######################
    ####### POINT A #######
    #######################
    alfaA[j] = alfa_fromStoA[-1,j]
    qA[j]    = q_fromStoA[-1,j]
    Point_A  = "PointA"
    cCl_A[:,j], cCd_A[:,j], AoA_tot_deg_A[j], AoA_tot_rad_A[j], cCz_A[:,j],\
    cCa_A[:,j], Norm_force_A[:,j], Axial_force_A[:,j], Shear_A[:,j],\
    Bending_A[:,j], m_distr_A[:,j], Torsion_A[:,j],\
    fig_A = obj2.CalcShearBendingTorsion_func(chord[:], y_halfspan[:],\
                                              cl_A[:,j], cd_A[:,j],\
                                              cm_A[:,j], alfaA[j], twist_angle,\
                                              qA[j], Point_A, n_fig, obj2)
    #######################
    ####### POINT C #######
    #######################
    if ( flight_env_case_pos == "Case1" ): 
        alfaC[j] = alfa_fromAtoGust1[-1,j]
        qC[j]    = q_fromAtoGust1[-1,j] 
    elif ( flight_env_case_pos == "Case2" ):
        alfaC[j] = alfa_fromAtoC[-1,j]
        qC[j]    = q_fromAtoC[-1,j]
    elif ( flight_env_case_pos == "Case3" ): 
        alfaC[j] = alfa_fromAtoGust1[-1,j]
        qC[j]    = q_fromAtoGust1[-1,j]  

    Point_C  = "PointC" 
    cCl_C[:,j], cCd_C[:,j], AoA_tot_deg_C[j], AoA_tot_rad_C[j], cCz_C[:,j],\
    cCa_C[:,j], Norm_force_C[:,j], Axial_force_C[:,j], Shear_C[:,j],\
    Bending_C[:,j], m_distr_C[:,j], Torsion_C[:,j],\
    fig_C = obj2.CalcShearBendingTorsion_func(chord[:], y_halfspan[:],\
                                              cl_C[:,j], cd_C[:,j],\
                                              cm_C[:,j], alfaC[j], twist_angle,\
                                              qC[j], Point_C, n_fig, obj2)     
    #######################
    ####### POINT D #######
    #######################
    alfaD[j] = alfa_fromDto0[0, j]
    qD[j]    = q_fromDto0[0, j]
    Point_D  = "PointD"   
    cCl_D[:,j], cCd_D[:,j], AoA_tot_deg_D[j], AoA_tot_rad_D[j], cCz_D[:,j],\
    cCa_D[:,j], Norm_force_D[:,j], Axial_force_D[:,j], Shear_D[:,j],\
    Bending_D[:,j], m_distr_D[:,j], Torsion_D[:,j],\
    fig_D = obj2.CalcShearBendingTorsion_func(chord[:], y_halfspan[:],\
                                              cl_D[:,j], cd_D[:,j],\
                                              cm_D[:,j], alfaD[j], twist_angle,\
                                              qD[j], Point_D, n_fig, obj2)      
    ################################
    ####### POINT S INVERTED #######
    ################################
    alfaSinv[j] = alfa_from0toSinv[0, j]
    qSinv[j]    = q_from0toSinv[0, j]
    Point_Sinv  = "PointSinverted"   
    cCl_Sinv[:,j], cCd_Sinv[:,j], AoA_tot_deg_Sinv[j], AoA_tot_rad_Sinv[j],\
    cCz_Sinv[:,j], cCa_Sinv[:,j], Norm_force_Sinv[:,j], Axial_force_Sinv[:,j],\
    Shear_Sinv[:,j], Bending_Sinv[:,j], m_distr_Sinv[:,j], Torsion_Sinv[:,j],\
    fig_Sinv = obj2.CalcShearBendingTorsion_func(chord[:], y_halfspan[:],\
                                                 cl_Sinv[:,j], cd_Sinv[:,j],\
                                                 cm_Sinv[:,j], alfaSinv[j],\
                                                 twist_angle, qSinv[j],\
                                                 Point_Sinv, n_fig, obj2)    
    #######################
    ####### POINT G #######
    #######################
    alfaG[j]    = alfa_fromSinvtoG[-1,j]
    qG[j]       = q_fromSinvtoG[-1,j]     
    Point_G  = "PointG"   
    cCl_G[:,j], cCd_G[:,j], AoA_tot_deg_G[j], AoA_tot_rad_G[j],\
    cCz_G[:,j], cCa_G[:,j], Norm_force_G[:,j], Axial_force_G[:,j],\
    Shear_G[:,j], Bending_G[:,j], m_distr_G[:,j], Torsion_G[:,j],\
    fig_G = obj2.CalcShearBendingTorsion_func(chord[:], y_halfspan[:],\
                                              cl_G[:,j], cd_G[:,j],\
                                              cm_G[:,j], alfaG[j], twist_angle,\
                                              qG[j], Point_G, n_fig, obj2)   
    #######################
    ####### POINT F #######
    #######################
    if ( neg_case_flag == "Case1_inverted" ): 
        alfaF[j] = alfa_fromGust1toF[-1,j]
        qF[j]    = q_fromGust1toF[-1,j] 
    elif ( neg_case_flag == "Case2_inverted" ):
        alfaF[j] = alfa_fromGtoF[-1,j]
        qF[j]    = q_fromGtoF[-1,j]
    elif ( neg_case_flag == "Case3_inverted" ): 
        alfaF[j] = alfa_fromGtoF[-1,j]
        qF[j]    = q_fromGtoF[-1,j]  
         
    Point_F  = "PointF"   
    cCl_F[:,j], cCd_F[:,j], AoA_tot_deg_F[j], AoA_tot_rad_F[j],\
    cCz_F[:,j], cCa_F[:,j], Norm_force_F[:,j], Axial_force_F[:,j],\
    Shear_F[:,j], Bending_F[:,j], m_distr_F[:,j], Torsion_F[:,j],\
    fig_F = obj2.CalcShearBendingTorsion_func(chord[:], y_halfspan[:],\
                                              cl_F[:,j], cd_F[:,j],\
                                              cm_F[:,j], alfaF[j], twist_angle,\
                                              qF[j], Point_F, n_fig, obj2)    
    #######################
    ####### POINT E #######
    #######################
    alfaE[j] = alfa_fromEto0[0, j]
    qE[j]    = q_fromEto0[0, j]
    Point_E  = "PointE"   
    cCl_E[:,j], cCd_E[:,j], AoA_tot_deg_E[j], AoA_tot_rad_E[j], cCz_E[:,j],\
    cCa_E[:,j], Norm_force_E[:,j], Axial_force_E[:,j], Shear_E[:,j],\
    Bending_E[:,j], m_distr_E[:,j], Torsion_E[:,j],\
    fig_E = obj2.CalcShearBendingTorsion_func(chord[:], y_halfspan[:],\
                                              cl_E[:,j], cd_E[:,j],\
                                              cm_E[:,j], alfaE[j], twist_angle,\
                                              qE[j], Point_E, n_fig, obj2) 
      
    # UPDATE FIG NUMBER 
    print("Figure number: " + str(n_fig))
    n_fig = n_fig + 1   
    
###############################################################################
############## STORE FIGURES INSIDE THE CORRECT FOLDER - S, M, T ##############
###############################################################################
fig_name = "\ShearBendingTorsion"
Point    = ["PointS", "PointA", "PointC", "PointD",\
            "PointSinverted", "PointG", "PointF", "PointE"]
fig_form = ".pdf"
main_dir = r'I:\PythonTesiConversion\csvla'
out_dir  = r'\Output\ShearBendingTorsion'
n        = 1

# Iteration - Moving figures inside Output folder
for i in range(0, len(Point)):
    # Print point
    print("\n")
    print("Point: " + Point[i])
    
    for j in range(0, n_Mass):
        figure_to_move = fig_name + Point[i] + str(n) + fig_form
        original       = main_dir + figure_to_move
        target         = main_dir + out_dir + figure_to_move
        shutil.move(original, target)
        
        # Figure number update
        print("Figure number: " + str(n))
        n = n + 1
        
    # REINITIALIZATION     
    print("\n")
    n = 1
#################################
####### CLOSE ALL FIGURES #######
#################################
plt.close("all")