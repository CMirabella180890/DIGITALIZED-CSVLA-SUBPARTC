# -*- coding: utf-8 -*-
"""
Created on Sun May  1 12:06:20 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  FINAL ENVELOPE                                                           %
%  Now we must develop a systematic approach to the final envelope, due to  %
%  the fact that the gust envelope and flight envelope must be combined to  %
%  obtain the final envelope.                                               %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@author: claum
"""
import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval
import csvla 
obj = csvla.csvla
from sympy import * 

#################################################
# ++++ POSITIVE SIDE OF THE FINAL ENVELOPE ++++ #
#################################################
N                 = 1000
VGust1            = np.ones((n_Mass,1))
VGust2            = np.ones((n_Mass,1))
nGust1            = np.ones((n_Mass,1))
nGust2            = np.ones((n_Mass,1))
new_VA            = np.ones((n_Mass,1))
new_nA            = np.ones((n_Mass,1))
V_fe_from0toS     = np.ones((N, n_Mass))
n_fe_from0toS     = np.ones((N, n_Mass))
V_fe_fromStoA     = np.ones((N, n_Mass))
n_fe_fromStoA     = np.ones((N, n_Mass))
n_fe_fromAtoGust1 = np.ones((N, n_Mass))
V_fe_fromAtoGust1 = np.ones((N, n_Mass))
n_fe_fromGust1toC = np.ones((N, n_Mass))
n_fe_fromAtoC     = np.ones((N, n_Mass))
V_fe_fromAtoC     = np.ones((N, n_Mass))
V_fe_fromGust1toC = np.ones((N, n_Mass))
n_fe_fromCtoGust2 = np.ones((N, n_Mass))
V_fe_fromCtoGust2 = np.ones((N, n_Mass))
n_fe_fromGust2toD = np.ones((N, n_Mass))
V_fe_fromGust2toD = np.ones((N, n_Mass))
n_fe_fromCtoD     = np.ones((N, n_Mass))
V_fe_fromCtoD     = np.ones((N, n_Mass))
n_fe_fromDto0     = np.ones((N, n_Mass))
V_fe_fromDto0     = np.ones((N, n_Mass))

# INVERTED FLIGHT
VGust1_inv        = np.ones((n_Mass,1))
VGust2_inv        = np.ones((n_Mass,1))
new_VG            = np.ones((n_Mass,1))
new_nG            = np.ones((n_Mass,1))
n_fe_from0toSinv  = np.ones((N, n_Mass))
V_fe_from0toSinv  = np.ones((N, n_Mass))
n_fe_fromSinvtoG  = np.ones((N, n_Mass))
V_fe_fromSinvtoG  = np.ones((N, n_Mass))
n_fe_fromGtoGust1 = np.ones((N, n_Mass))
V_fe_fromGtoGust1 = np.ones((N, n_Mass))
n_fe_fromGust1toF = np.ones((N, n_Mass))
V_fe_fromGust1toF = np.ones((N, n_Mass))
n_fe_fromGtoF     = np.ones((N, n_Mass))
V_fe_fromGtoF     = np.ones((N, n_Mass))
n_fe_fromFtoE     = np.ones((N, n_Mass))
V_fe_fromFtoE     = np.ones((N, n_Mass))
n_fe_fromEto0     = np.ones((N, n_Mass))
V_fe_fromEto0     = np.ones((N, n_Mass))

for i in range(0, n_Mass):
    
    # SYMBOLIC VARIABLES INITIALIZATION 
    V = symbols("V")
    
    # ASSIGNING CORRESPONDING VALUES 
    a  = float(rho0 * CL_max_clean)
    b  = float(rho0 * CLalfa_rad * kg_op[i] * Ude_cruise)
    c  = float(2 * WS[i])
    eq = a * V**2 - b * V - c
    
    # SOLVING THE EQUATION
    solution_pos = solveset( Eq(eq, 0), V )
    VA1          = abs(float(solution_pos.args[0]))
    VA2          = abs(float(solution_pos.args[1]))
    
    # TESTING OUTPUTS
    print(" -- Checking: \n", " -- Case" + str(i), ":\n")
    print(" -- New VA1: ", VA1)
    print(" -- New VA2: ", VA2, "\n")
    
    for k in range(len(solution_pos.args)): 
        if (VA1 > 0) and (VA1 > VA2):
            new_VA[i] = VA1 
        elif (VA2 > 0) and (VA2 > VA1):
            new_VA[i] = VA2
            
    # EVALUATING THE COMBINED FLIGHT ENVELOPE - FINAL ENVELOPE
    if (new_VA[i] > VA[i]):
        
        # FINAL ENVELOPE - FROM 0 TO S 
        V_fe_from0toS[:,i] = np.full(N, VS[i])
        n_fe_from0toS[:,i] = np.reshape(np.linspace(0.0, nS[i], N), (N,))
        
        # FINAL ENVELOPE - FROM S TO A 
        n_fe_fromStoA[:,i] = np.reshape(np.linspace(nS[i], nA[i], N), (N,))
        V_fe_fromStoA[:,i] = obj.calcvs(rho0, WS[i], CL_max_clean, n_fe_fromStoA[:,i])
        
        # ASSESSING LOAD FACTOR CASE 
        if ( np.max(ng_pos_cruise_op[:,i]) > nmax ):
            
            # IDENTIFICATION FLAG - POSITIVE SIDE
            flight_env_case_pos = "Case1"
            
            ################# NEW N_GUST1 VALUE #################
            #                   0.5 * rho0 * V * a * Kg * Ude   #
            #    n_gust = 1.0 + -----------------------------   #
            #                             ( W / S )             #
            #####################################################
            n_gust1          = symbols("n_gust1")
            a                = float(0.5 * rho0 * V_fe_fromStoA[-1,i] * CLalfa_rad * kg_op[i] * Ude_cruise)
            b                = float(WS[i])
            eq               = n_gust1 - (a / b) - 1.0
            solution_n_gust1 = solveset( Eq(eq, 0), n_gust1 )
            nGust1[i]        = float(solution_n_gust1.args[0])
            
            # FINAL ENVELOPE FROM A TO GUST1 
            n_fe_fromAtoGust1[:,i] = np.reshape(np.linspace(nA[i], nGust1[i], N), (N,))
            V_fe_fromAtoGust1[:,i] = obj.calcvs(rho0, WS[i], CL_max_clean, n_fe_fromAtoGust1[:,i])
            
            # FINAL ENVELOPE FROM GUST1 TO C 
            V_fe_fromGust1toC[:,i]  = np.reshape(np.linspace(V_fe_fromAtoGust1[-1,i], VC[i], N), (N,))
            n_function_fromGust1toC = lambdify(n_gust1, eq, "numpy")
            n_fe_fromGust1toC[:,i]  = n_function_fromGust1toC(V_fe_fromGust1toC[:,i])
            
            # FINAL ENVELOPE FROM C TO GUST2
            V_gust2          = symbols("V_gust2")
            y1_pos           = n_fe_fromGust1toC[-1,i]
            y2_pos           = nmax
            y_pos            = np.reshape(np.array([ y1_pos, y2_pos ]), (2,))
            x_pos            = np.reshape(np.array([ VC[i], VD[i]] ), (2,))
            p_fromCtoGust2   = chebfit(x_pos, y_pos, deg=1)
            a                = p_fromCtoGust2[1]
            b                = p_fromCtoGust2[0]
            eq               = a * V_gust2 + b - nmax
            solution_V_gust2 = solveset( Eq(eq,0), V_gust2 )
            VGust2[i]        = float(solution_V_gust2.args[0])
            # ACTUAL CALCULATION 
            V_fe_fromCtoGust2[:,i]  = np.reshape(np.linspace(VC[i], VGust2[i], N), (N,))
            n_function_fromCtoGust2 = lambdify(V_gust2, eq, "numpy")
            n_fe_fromCtoGust2[:,i]  = n_function_fromCtoGust(V_fe_fromCtoGust2[:,i])
            
            # FINAL ENVELOPE FROM GUST2 TO D 
            n_fe_fromGust2toD[:,i] = np.full(N, n_fe_fromCtoGust2[-1,i])
            V_fe_fromGust2toD[:,i] = np.reshape(np.linspace(V_fe_fromCtoGust2[-1,i], VD[i], N),(N,))
            
            # FINAL ENVELOPE FROM D TO 0 
            V_fe_fromDto0[:,i] = np.full(N, VD[i])
            n_fe_fromDto0[:,i] = np.reshape(np.linspace(n_fe_fromGust2toD[-1,i], 0.0, N), (N,))
            
        elif ( np.max(ng_pos_cruise_op[:,i]) < nmax ):

            # IDENTIFICATION FLAG - POSITIVE SIDE
            flight_env_case_pos = "Case2"
            
            # FINAL ENVELOPE FROM A TO C
            V_fe_fromAtoC[:,i] = np.reshape(np.linspace(VA[i], VC[i], N),(N,))
            n_fe_fromAtoC[:,i] = np.full(N, nmax)
            
            # FINAL ENVELOPE FROM C TO D 
            V_fe_fromCtoD[:,i] = np.reshape(np.linspace(VC[i], VD[i], N), (N,) )
            n_fe_fromCtoD[:,i] = np.full(N, nmax)
            
            # FINAL ENVELOPE FROM D TO 0 
            V_fe_fromDto0[:,i] = np.full(N, VD[i])
            n_fe_fromDto0[:,i] = np.reshape(np.linspace(nmax, 0.0, N), (N,))
        
    elif (new_VA[i] < VA[i]):
        
        # IDENTIFICATION FLAG - POSITIVE SIDE
        flight_env_case_pos = "Case3" 
            
        ################### NEW N_A VALUE ###################
        #                   0.5 * rho0 * V * a * Kg * Ude   #
        #    n_gust = 1.0 + -----------------------------   #
        #                             ( W / S )             #
        #####################################################
        nA_sym = symbols("nA")
        a      = float(0.5 * rho0 * new_VA[i] * CLalfa_rad * kg_op[i] * Ude_cruise)
        b      = float( WS[i] )
        eq     = b * nA_sym - b - a
        
        # SOLVING FOR N_A
        solution_nA = solveset( Eq(eq, 0), nA_sym )
        new_nA[i]   = float(solution_nA.args[0])
        
        # FINAL ENVELOPE - FROM 0 TO S
        n_fe_from0toS[:,i] = np.reshape(np.linspace(0.0, nS[i], N), (N,))
        V_fe_from0toS[:,i] = np.full(N, VS[i])
        
        # FINAL ENVELOPE FROM S TO A
        n_fe_fromStoA[:,i] = np.reshape(np.linspace(nS[i], nA[i], N), (N,))
        V_fe_fromStoA[:,i] = obj.calcvs(rho0, WS[i], CL_max_clean, n_fe_fromStoA[:,i])
        
        # FINAL ENVELOPE FROM A TO GUST1 
        V_gust1 = symbols("V_gust1")
        a       = float(0.5 * rho0 * CLalfa_rad * kg_op[i] * Ude_cruise)
        b       = float(WS[i] - nmax * WS[i]) 
        eq      = a * V_gust1 + b
        
        # SOLVING FOR GUST1 AIRSPEED 
        solution_V_gust1       = solveset( Eq(eq, 0), V_gust1) 
        VGust1[i]              = float(solution_V_gust1.args[0])
        n_fe_fromAtoGust1[:,i] = np.full(N, nA[i])
        V_fe_fromAtoGust1[:,i] = np.reshape(np.linspace(VA[i], VGust1[i], N), (N,))
        
        # FINAL ENVELOPE FROM GUST1 TO C
        V_fe_fromGust1toC[:,i] = np.reshape(np.linspace(VGust1[i], VC[i], N), (N,))
        n_fe_fromGust1toC[:,i] = obj.calc_n_gust(rho0, V_fe_fromGust1toC[:,i],\
                                                V_fe_fromGust1toC[:,i], CLalfa_rad,\
                                                kg_op[i], Ude_cruise, WS[i],\
                                                        pos_case_flag)
            
        # FINAL ENVELOPE FROM C TO GUST2 
        x_pos           = np.reshape(np.array([ VC[i], VD[i]] ), (2,))
        y1_pos          = n_fe_fromGust1toC[-1,i]
        y2_pos          = ng_pos_dive_op[-1,i]
        y_pos           = np.reshape(np.array([ y1_pos, y2_pos ]), (2,))
        p_gust_positive = chebfit(x_pos, y_pos, deg = 1)
        V_gust2         = symbols("V_gust2")
        a               = p_gust_positive[1]
        b               = p_gust_positive[0]
        eq              = a * V_gust2 + b 
        
        # SOLVING FOR GUST2 AIRSPEED
        solution_V_gust2       = solveset( Eq(eq, nmax), V_gust2) 
        VGust2[i]              = float(solution_V_gust2.args[0])
        V_fe_fromCtoGust2[:,i] = np.reshape(np.linspace(VC[i], VGust2[i], N), (N,))
        n_fe_fromCtoGust2[:,i] = chebval(V_fe_fromCtoGust2[:,i], p_gust_positive)
        
        # FINAL ENVELOPE FROM GUST2 TO D 
        V_fe_fromGust2toD[:,i] = np.reshape(np.linspace(VGust2[i], VD[i], N), (N,))
        n_fe_fromGust2toD[:,i] = np.full(N, nmax)
        
        # FINAL ENVELOPE FROM D TO 0 
        V_fe_fromDto0[:,i] = np.full(N, VD[i])
        n_fe_fromDto0[:,i] = np.reshape(np.linspace(nmax, 0.0, N), (N,))

#################################################
# ++++ INVERTED SIDE OF THE FINAL ENVELOPE ++++ #
#################################################        
for i in range(0, n_Mass):
    
    # SYMBOLIC VARIABLES INITIALIZATION 
    V = symbols("V")
    
    # ASSIGNING CORRESPONDING VALUES 
    a  = float(rho0 * abs(CL_inverted))
    b  = float(rho0 * CLalfa_rad * kg_op[i] * Ude_cruise)
    c  = float(2 * WS[i])
    eq = a * V**2 + b * V - c
    
    # SOLVING THE EQUATION
    solution_pos = solveset( Eq(eq, 0), V )
    VG1          = abs(float(solution_pos.args[0]))
    VG2          = abs(float(solution_pos.args[1]))
    
    # TESTING OUTPUTS
    print(" -- Checking: \n", " -- Case" + str(i), ":\n")
    print(" -- New VG1: ", VG1)
    print(" -- New VG2: ", VG2, "\n")
    
    for k in range(len(solution_pos.args)): 
        if (VG1 > 0) and (VG1 > VG2):
            new_VG[i] = VG1 
        elif (VG2 > 0) and (VG2 > VG1):
            new_VG[i] = VG2        
            
    # INVERTED SIDE DECISION LOOP 
    if (new_VG[i] > VG[i]):
        
        # IDENTIFICATION FLAG - POSITIVE SIDE
        flight_env_case_neg = "Case1_inverted" 
        
        # FINAL ENVELOPE FROM 0 TO S INVERTED 
        n_fe_from0toSinv[:,i] = np.reshape(np.linspace(0.0, nSinv[i], N), (N,))
        V_fe_from0toSinv[:,i] = np.full(N, VSinv[i])
        
        # FINAL ENVELOPE FROM S INVERTED TO G 
        n_fe_fromSinvtoG[:,i] = np.reshape(np.linspace(nSinv[i], nG[i], N), (N,))
        V_fe_fromSinvtoG[:,i] = obj.calcvs(rho0, WS[i], CL_inverted, n_fe_fromSinvtoG[:,i])
        
        # FINAL ENVELOPE FROM G TO GUST1 
        ################### NEW N_GUST1 VALUE ###################
        #                   0.5 * rho0 * V * a * Kg * Ude       #
        #    n_gust = 1.0 + -----------------------------       #
        #                             ( W / S )                 #
        #########################################################
        V_gust1          = symbols("V_gust1")
        a                = float( 0.5 * rho0 * CLalfa_rad * kg_op[i] * Ude_cruise )
        b                = float( WS[i] * ( 1 - nmin ))
        eq               = a * V_gust1 + b
        solution_V_gust1 = solveset( Eq(eq, 0), V_gust1 )
        VGust1_inv[i]    = abs(float(solution_V_gust1.args[0]))
        
        # FROM G TO GUST1 
        V_fe_fromGtoGust1[:,i] = np.reshape(np.linspace( V_fe_fromSinvtoG[-1,i], VGust1_inv[i], N ), (N,))
        n_fe_fromGtoGust1[:,i] = np.full(N, nmin)
        
        # FINAL ENVELOPE FROM GUST1 TO F
        V_fe_fromGust1toF[:,i] = np.reshape(np.linspace(VGust1_inv[i], VF[i], N), (N,))
        n_fe_fromGust1toF[:,i] = obj.calc_n_gust(rho0, V_fe_fromGust1toF[:,i],\
                                                V_fe_fromGust1toF[:,i], CLalfa_rad,\
                                                kg_op[i], Ude_cruise, WS[i],\
                                                        neg_case_flag)
        
        # FINAL ENVELOPE FROM F TO E 
        x_neg           = np.reshape(np.array([ VF[i], VE[i]] ), (2,))
        y1_neg          = n_fe_fromGust1toF[-1,i]
        y2_neg          = ng_neg_dive_op[-1,i]
        y_neg           = np.reshape(np.array([ y1_neg, y2_neg ]), (2,))
        p_gust_negative = chebfit(x_neg, y_neg, deg = 1)
        V_gust2         = symbols("V_gust2")
        a               = p_gust_negative[1]
        b               = p_gust_negative[0]
        eq              = a * V_gust2 + b 
        
        # SOLVING FOR GUST2 AIRSPEED
        solution_V_gust2   = solveset( Eq(eq, nmin), V_gust2) 
        VGust2_inv[i]      = float(solution_V_gust2.args[0])
        V_fe_fromFtoE[:,i] = np.reshape(np.linspace(VF[i], VE[i], N), (N,))
        n_fe_fromFtoE[:,i] = chebval(V_fe_fromFtoE[:,i], p_gust_negative) 
        
        # FINALE ENVELOPE FROM E TO 0
        V_fe_fromEto0[:,i] = np.full(N, VE[i])
        n_fe_fromEto0[:,i] = np.reshape(np.linspace(n_fe_fromFtoE[-1,i], 0.0, N), (N,))
        
        
    elif (new_VG[i] < VG[i]):
        # ASSESSING LOAD FACTOR CASE 
        if ( np.abs(np.max(ng_neg_cruise_op[:,i])) > abs(nmin) ):
            
            # IDENTIFICATION FLAG - POSITIVE SIDE
            flight_env_case_neg = "Case2_inverted" 
            
            # FINAL ENVELOPE FROM 0 TO S INVERTED 
            n_fe_from0toSinv[:,i] = np.reshape(np.linspace(0.0, nSinv[i], N), (N,))
            V_fe_from0toSinv[:,i] = np.full(N, VSinv[i])
            
            # FINAL ENVELOPE FROM S INVERTED TO G 
            new_nG[i] = obj.calc_n_gust(rho0, new_VG[i],\
                                        new_VG[i], CLalfa_rad,\
                                        kg_op[i], Ude_cruise, WS[i],\
                                                         neg_case_flag)
            n_fe_fromSinvtoG[:,i] = np.reshape(np.linspace(nSinv[i], new_nG[i], N), (N,))
            V_fe_fromSinvtoG[:,i] = obj.calcvs(rho0, WS[i], CL_inverted, n_fe_fromSinvtoG[:,i])
            
            # FINAL ENVELOPE FROM G TO F 
            V_fe_fromGtoF[:,i] = np.reshape(np.linspace(new_VG[i], VF[i], N), (N,))
            n_fe_fromGtoF[:,i] = obj.calc_n_gust(rho0, V_fe_fromGtoF[:,i],\
                                        V_fe_fromGtoF[:,i], CLalfa_rad,\
                                        kg_op[i], Ude_cruise, WS[i],\
                                                         neg_case_flag)
                
            # FINAL ENVELOPE FROM F TO E
            x_neg           = np.reshape(np.array([ V_fe_fromGtoF[-1,i], VE[i]] ), (2,))
            y1_neg          = n_fe_fromGust1toF[-1,i]
            y2_neg          = ng_neg_dive_op[-1,i]
            y_neg           = np.reshape(np.array([ y1_neg, y2_neg ]), (2,))
            p_gust_negative = chebfit(x_neg, y_neg, deg = 1)
            V_gust2         = symbols("V_gust2")
            a               = p_gust_negative[1]
            b               = p_gust_negative[0]
            eq              = a * V_gust2 + b 
            
            # SOLVING FOR GUST2 AIRSPEED
            solution_V_gust2   = solveset( Eq(eq, nmin), V_gust2) 
            VGust2_inv[i]      = float(solution_V_gust2.args[0])
            V_fe_fromFtoE[:,i] = np.reshape(np.linspace(V_fe_fromGtoF[-1,i], VE[i], N), (N,))
            n_fe_fromFtoE[:,i] = chebval(V_fe_fromFtoE[:,i], p_gust_negative) 
            
            # FINALE ENVELOPE FROM E TO 0
            V_fe_fromEto0[:,i] = np.full(N, VE[i])
            n_fe_fromEto0[:,i] = np.reshape(np.linspace(n_fe_fromFtoE[-1,i], 0.0, N), (N,))
            
            

        elif ( np.abs(np.max(ng_neg_cruise_op[:,i])) < abs(nmin) ):   
            
            # IDENTIFICATION FLAG - POSITIVE SIDE
            flight_env_case_neg = "Case3_inverted"  
            
            # FINAL ENVELOPE FROM 0 TO S INVERTED 
            n_fe_from0toSinv[:,i] = np.reshape(np.linspace(0.0, nSinv[i], N), (N,))
            V_fe_from0toSinv[:,i] = np.full(N, VSinv[i])
            
            # FINAL ENVELOPE FROM S INVERTED TO G 
            n_fe_fromSinvtoG[:,i] = np.reshape(np.linspace(nSinv[i], nG[i], N), (N,))
            V_fe_fromSinvtoG[:,i] = obj.calcvs(rho0, WS[i], CL_inverted, n_fe_fromSinvtoG[:,i])  
            
            # FINAL ENVELOPE FROM G TO F 
            n_fe_fromGtoF[:,i] = np.reshape(np.linspace(nG[i], nF[i], N), (N,))
            V_fe_fromGtoF[:,i] = np.reshape(np.linspace(V_fe_fromSinvtoG[-1,i], VF[i], N), (N,))
            
            # FINAL ENVELOPE FROM F TO E 
            n_fe_fromFtoE[:,i] = np.reshape(np.linspace(nF[i], nE[i], N), (N,))
            V_fe_fromFtoE[:,i] = np.reshape(np.linspace(V_fe_fromGtoF[-1,i], VE[i], N), (N,))
            
            # FINAL ENVELOPE FROM E TO 0
            n_fe_fromEto0[:,i] = np.reshape(np.linspace(nE[i], 0.0, N), (N,))
            V_fe_fromEto0[:,i] = np.full(N, VE[i])
n = 1            
for i in range(0,n_Mass):       
        figure_name = "fig" + str(n)          
        figure_name = obj.final_envelope(V_fe_from0toS[:,i], n_fe_from0toS[:,i], V_fe_fromStoA[:,i],\
                           n_fe_fromStoA[:,i], n_fe_fromAtoGust1[:,i], V_fe_fromAtoGust1[:,i],\
                           n_fe_fromGust1toC[:,i], V_fe_fromGust1toC[:,i], n_fe_fromCtoGust2[:,i],\
                           V_fe_fromCtoGust2[:,i], n_fe_fromGust2toD[:,i], V_fe_fromGust2toD[:,i],\
                           n_fe_fromDto0[:,i], V_fe_fromDto0[:,i], n_fe_from0toSinv[:,i],\
                           V_fe_from0toSinv[:,i], n_fe_fromSinvtoG[:,i], V_fe_fromSinvtoG[:,i],\
                           n_fe_fromGtoGust1[:,i], V_fe_fromGtoGust1[:,i], n_fe_fromGust1toF[:,i],\
                           V_fe_fromGust1toF[:,i], n_fe_fromFtoE[:,i], V_fe_fromFtoE[:,i],\
                           n_fe_fromEto0[:,i], V_fe_fromEto0[:,i], n_fe_fromGtoF[:,i],\
                           V_fe_fromGtoF[:,i], n_fe_fromAtoC[:,i], V_fe_fromAtoC[:,i],\
                           n_fe_fromCtoD[:,i], V_fe_fromCtoD[:,i], Reg, Aircraft_name, n,\
                           flight_env_case_pos, flight_env_case_neg)
        n = n + 1
# =============================================================================    
# MOVING FIGURES INSIDE OUTPUT FOLDER
# FinalEnvelope
n = 1
for i in range(0, n_Mass):
    figure_name = r'\FinalEnvelope' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\FinalEnvelope'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)
    n = n + 1
# ============================================================================= 
###############################   
###### CLOSE ALL FIGURES ######
###############################
plt.close("all")           
####################################################
##### STORE INSIDE THE SIMPLE NAMESPACE OBJECT #####
####################################################
Final_envelope = obj.final_envelope_store_points(V_fe_from0toS, n_fe_from0toS, V_fe_fromStoA,\
                   n_fe_fromStoA, n_fe_fromAtoGust1, V_fe_fromAtoGust1,\
                   n_fe_fromGust1toC, V_fe_fromGust1toC, n_fe_fromCtoGust2,\
                   V_fe_fromCtoGust2, n_fe_fromGust2toD, V_fe_fromGust2toD,\
                   n_fe_fromDto0, V_fe_fromDto0, n_fe_from0toSinv,\
                   V_fe_from0toSinv, n_fe_fromSinvtoG, V_fe_fromSinvtoG,\
                   n_fe_fromGtoGust1, V_fe_fromGtoGust1, n_fe_fromGust1toF,\
                   V_fe_fromGust1toF, n_fe_fromFtoE, V_fe_fromFtoE,\
                   n_fe_fromEto0, V_fe_fromEto0, n_fe_fromGtoF,\
                   V_fe_fromGtoF, n_fe_fromAtoC, V_fe_fromAtoC,\
                   n_fe_fromCtoD, V_fe_fromCtoD, flight_env_case_pos,\
                   flight_env_case_neg)
# UPDATING DICTIONARY
aircraft.update(Final_envelope)
 # UPDATING THE SIMPLE NAMESPACE OBJECT
aircraft_data = SimpleNamespace(**aircraft)           