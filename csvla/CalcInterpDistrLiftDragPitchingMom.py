# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:40:30 2022

@author: claum
"""
def plot_distributions(y_halfspan, unit_CL, unit_CD, unit_CM,\
                       cl_S, cl_A, cl_C, cl_D, cl_Sinv, cl_G, cl_F, cl_E,\
                       cd_S, cd_A, cd_C, cd_D, cd_Sinv, cd_G, cd_F, cd_E,\
                       cm_S, cm_A, cm_C, cm_D, cm_Sinv, cm_G, cm_F, cm_E, n): 
    
    # =======================================
    # ======= COEFFICIENTS =======
    # =======================================
    cl_comparison = plt.figure()
    plt.plot(y_halfspan, cl_S, color="#000000", linewidth=1.0, linestyle = "solid", label = r"Point S")
    plt.plot(y_halfspan, cl_A, color="#0000FF", linewidth=1.0, linestyle = "solid", label = r"Point A")
    plt.plot(y_halfspan, cl_C, color="#FF0000", linewidth=1.0, linestyle = "solid", label = r"Point C")
    plt.plot(y_halfspan, cl_D, color="#FF4500", linewidth=1.0, linestyle = "solid", label = r"Point D")
    
    plt.plot(y_halfspan, cl_Sinv, color="#9A0EEA", linewidth=1.0, linestyle = "solid", label = r"Point S inv")
    plt.plot(y_halfspan, cl_G,    color="#800000", linewidth=1.0, linestyle = "solid", label = r"Point G")
    plt.plot(y_halfspan, cl_F,    color="#008000", linewidth=1.0, linestyle = "solid", label = r"Point F")
    plt.plot(y_halfspan, cl_E,    color="#008080", linewidth=1.0, linestyle = "solid", label = r"Point E")
    
    plt.plot(y_halfspan, unit_CL, color="red", linewidth=1.0, linestyle = "dashed", label = r"Unitary load - $c_{l}(y)$")
    
    plt.xlabel(r'Station along wing semispan - $y$ [m]')                # x-label to the axes.
    plt.ylabel(r'Spanwise lift coefficient distribution - $C_{l}(y)$')  # y-label to the axes.
    plt.title(r'Lift coefficient curve distributions - Comparison') # Title to the axes.
    plt.minorticks_on()
    plt.grid(True, linestyle='-.', which="both")
    plt.legend(loc="best", prop={'size': 6})
    figure_name = 'cl_comparison' + str(n) + '.pdf'
    plt.savefig(figure_name, bbox_inches='tight')
    plt.show()
    # =======================================
    
    # =======================================
    cd_comparison = plt.figure()
    plt.plot(y_halfspan, cd_S, color="#000000", linewidth=1.0, linestyle = "solid", label = r"Point S")
    plt.plot(y_halfspan, cd_A, color="#0000FF", linewidth=1.0, linestyle = "solid", label = r"Point A")
    plt.plot(y_halfspan, cd_C, color="#FF0000", linewidth=1.0, linestyle = "solid", label = r"Point C")
    plt.plot(y_halfspan, cd_D, color="#FF4500", linewidth=1.0, linestyle = "solid", label = r"Point D")
    
    plt.plot(y_halfspan, cd_Sinv, color="#9A0EEA", linewidth=1.0, linestyle = "solid", label = r"Point S inv")
    plt.plot(y_halfspan, cd_G,    color="#800000", linewidth=1.0, linestyle = "solid", label = r"Point G")
    plt.plot(y_halfspan, cd_F,    color="#008000", linewidth=1.0, linestyle = "solid", label = r"Point F")
    plt.plot(y_halfspan, cd_E,    color="#008080", linewidth=1.0, linestyle = "solid", label = r"Point E")
    
    # plt.plot(y_halfspan, unit_CD, color="red", linewidth=1.0, linestyle = "dashed", label = r"Unitary load - $c_{d}(y)$")
    
    plt.xlabel(r'Station along wing semispan - $y$ [m]')                # x-label to the axes.
    plt.ylabel(r'Spanwise lift coefficient distribution - $C_{d}(y)$')  # y-label to the axes.
    plt.title(r'Drag coefficient curve distributions - Comparison') # Title to the axes.
    plt.minorticks_on()
    plt.grid(True, linestyle='-.', which="both")
    plt.legend(loc="best", prop={'size': 6})
    figure_name = 'cd_comparison' + str(n) + '.pdf'
    plt.savefig(figure_name, bbox_inches='tight')
    plt.show()
    # =======================================
    
    # =======================================
    cm_comparison = plt.figure()
    plt.plot(y_halfspan, cm_S, color="#000000", linewidth=1.0, linestyle = "solid", label = r"Point S")
    plt.plot(y_halfspan, cm_A, color="#0000FF", linewidth=1.0, linestyle = "solid", label = r"Point A")
    plt.plot(y_halfspan, cm_C, color="#FF0000", linewidth=1.0, linestyle = "solid", label = r"Point C")
    plt.plot(y_halfspan, cm_D, color="#FF4500", linewidth=1.0, linestyle = "solid", label = r"Point D")
    
    plt.plot(y_halfspan, cm_Sinv, color="#9A0EEA", linewidth=1.0, linestyle = "solid", label = r"Point S inv")
    plt.plot(y_halfspan, cm_G,    color="#800000", linewidth=1.0, linestyle = "solid", label = r"Point G")
    plt.plot(y_halfspan, cm_F,    color="#008000", linewidth=1.0, linestyle = "solid", label = r"Point F")
    plt.plot(y_halfspan, cm_E,    color="#008080", linewidth=1.0, linestyle = "solid", label = r"Point E")
    
    # plt.plot(y_halfspan, unit_CM, color="red", linewidth=1.0, linestyle = "dashed", label = r"Unitary load - $c_{m}(y)$")
    
    plt.xlabel(r'Station along wing semispan - $y$ [m]')                # x-label to the axes.
    plt.ylabel(r'Spanwise pitch. mom. coefficient distribution - $C_{m}(y)$')  # y-label to the axes.
    plt.title(r'Pitch. mom coefficient curve distributions - Comparison') # Title to the axes.
    plt.minorticks_on()
    plt.grid(True, linestyle='-.', which="both")
    plt.legend(loc="best", prop={'size': 6})
    figure_name = 'cm_comparison' + str(n) + '.pdf'
    plt.savefig(figure_name, bbox_inches='tight')
    plt.show()
    # =======================================
    
    return cl_comparison, cd_comparison, cm_comparison

# STRAIGHT FLIGHT INITIALIZATION

###############################
###### LIFT COEFFICIENTS ######
###############################
CL_S    = np.ones((n_Mass,1))
cl_S    = np.ones((N, n_Mass))
CL_A    = np.ones((n_Mass,1))
cl_A    = np.ones((N, n_Mass))
CL_C    = np.ones((n_Mass,1))
cl_C    = np.ones((N, n_Mass))
CL_D    = np.ones((n_Mass,1))
cl_D    = np.ones((N, n_Mass))
CL_Sinv = np.ones((n_Mass,1))
cl_Sinv = np.ones((N, n_Mass))
CL_G    = np.ones((n_Mass,1))
cl_G    = np.ones((N, n_Mass))
CL_F    = np.ones((n_Mass,1))
cl_F    = np.ones((N, n_Mass))
CL_E    = np.ones((n_Mass,1))
cl_E    = np.ones((N, n_Mass))

###############################
###### DRAG COEFFICIENTS ######
###############################
CD_S    = np.ones((n_Mass,1))
cd_S    = np.ones((N, n_Mass))
CD_A    = np.ones((n_Mass,1))
cd_A    = np.ones((N, n_Mass))
CD_C    = np.ones((n_Mass,1))
cd_C    = np.ones((N, n_Mass))
CD_D    = np.ones((n_Mass,1))
cd_D    = np.ones((N, n_Mass))
CD_Sinv = np.ones((n_Mass,1))
cd_Sinv = np.ones((N, n_Mass))
CD_G    = np.ones((n_Mass,1))
cd_G    = np.ones((N, n_Mass))
CD_F    = np.ones((n_Mass,1))
cd_F    = np.ones((N, n_Mass))
CD_E    = np.ones((n_Mass,1))
cd_E    = np.ones((N, n_Mass))

##########################################
###### PITCHING MOMENT COEFFICIENTS ######
##########################################
CM_S    = np.ones((n_Mass,1))
cm_S    = np.ones((N, n_Mass))
CM_A    = np.ones((n_Mass,1))
cm_A    = np.ones((N, n_Mass))
CM_C    = np.ones((n_Mass,1))
cm_C    = np.ones((N, n_Mass))
CM_D    = np.ones((n_Mass,1))
cm_D    = np.ones((N, n_Mass))
CM_Sinv = np.ones((n_Mass,1))
cm_Sinv = np.ones((N, n_Mass))
CM_G    = np.ones((n_Mass,1))
cm_G    = np.ones((N, n_Mass))
CM_F    = np.ones((n_Mass,1))
cm_F    = np.ones((N, n_Mass))
CM_E    = np.ones((n_Mass,1))
cm_E    = np.ones((N, n_Mass))

#############################
###### STRAIGHT FLIGHT ######
#############################
if ( flight_env_case_pos == "Case1" ): 
    
    for j in range(0,n_Mass):
        #################################
        ############ POINT S ############
        #################################
        CL_S[j] = CL_wb_from0toS[-1,j]
        CD_S[j] = CD_from0toS[-1,j]
        CM_S[j] = CM_CG_from0toS[-1,j]
        
        cl_S[:,j] = float(CL_S[j]) * unit_CL
        cd_S[:,j] = float(CD_S[j]) * unit_CD
        cm_S[:,j] = abs(float(CM_S[j])) * unit_CM
    
        #################################
        ############ POINT A ############
        #################################
        CL_A[j] = CL_wb_fromStoA[-1,j]
        CD_A[j] = CD_fromStoA[-1,j]
        CM_A[j] = CM_CG_fromStoA[-1,j]
        
        cl_A[:,j] = float(CL_A[j]) * unit_CL
        cd_A[:,j] = float(CD_A[j]) * unit_CD
        cm_A[:,j] = abs(float(CM_A[j])) * unit_CM
    
        #################################
        ############ POINT C ############
        #################################
        CL_C[j] = CL_wb_fromGust1toC[-1,j]
        CD_C[j] = CD_fromGust1toC[-1,j]
        CM_C[j] = CM_CG_fromGust1toC[-1,j]
        
        cl_C[:,j] = float(CL_C[j]) * unit_CL
        cd_C[:,j] = float(CD_C[j]) * unit_CD
        cm_C[:,j] = abs(float(CM_C[j])) * unit_CM
    
        #################################
        ############ POINT D ############
        #################################
        CL_D[j] = CL_wb_fromGust2toD[-1,j]
        CD_D[j] = CD_fromGust2toD[-1,j]
        CM_D[j] = CM_CG_fromGust2toD[-1,j]
        
        cl_D[:,j] = float(CL_D[j]) * unit_CL
        cd_D[:,j] = float(CD_D[j]) * unit_CD
        cm_D[:,j] = abs(float(CM_D[j])) * unit_CM
    
elif ( flight_env_case_pos == "Case2" ):
    
    for j in range(0,n_Mass):
        #################################
        ############ POINT S ############
        #################################
        CL_S[j] = CL_wb_from0toS[-1,j]
        CD_S[j] = CD_from0toS[-1,j]
        CM_S[j] = CM_CG_from0toS[-1,j]
        
        cl_S[:,j] = float(CL_S[j]) * unit_CL
        cd_S[:,j] = float(CD_S[j]) * unit_CD
        cm_S[:,j] = abs(float(CM_S[j])) * unit_CM
        
        #################################
        ############ POINT A ############
        #################################
        CL_A[j] = CL_wb_fromStoA[-1,j]
        CD_A[j] = CD_fromStoA[-1,j]
        CM_A[j] = CM_CG_fromStoA[-1,j]
        
        cl_A[:,j] = float(CL_A[j]) * unit_CL
        cd_A[:,j] = float(CD_A[j]) * unit_CD
        cm_A[:,j] = abs(float(CM_A[j])) * unit_CM
        
        #################################
        ############ POINT C ############
        #################################
        CL_C[j] = CL_wb_fromAtoC[-1,j]
        CD_C[j] = CD_fromAtoC[-1,j]
        CM_C[j] = CM_CG_fromAtoC[-1,j]
        
        cl_C[:,j] = float(CL_C[j]) * unit_CL
        cd_C[:,j] = float(CD_C[j]) * unit_CD
        cm_C[:,j] = abs(float(CM_C[j])) * unit_CM
        
        #################################
        ############ POINT D ############
        #################################
        CL_D[j] = CL_wb_fromCtoD[-1,j]
        CD_D[j] = CD_fromCtoD[-1,j]
        CM_D[j] = CM_CG_fromCtoD[-1,j]
        
        cl_D[:,j] = float(CL_D[j]) * unit_CL
        cd_D[:,j] = float(CD_D[j]) * unit_CD
        cm_D[:,j] = abs(float(CM_D[j])) * unit_CM
        
elif ( flight_env_case_pos == "Case3" ): 
    
    for j in range(0,n_Mass):
    
        #################################
        ############ POINT S ############
        #################################
        CL_S[j] = CL_wb_from0toS[-1,j]
        CD_S[j] = CD_from0toS[-1,j]
        CM_S[j] = CM_CG_from0toS[-1,j]
        
        cl_S[:,j] = float(CL_S[j]) * unit_CL
        cd_S[:,j] = float(CD_S[j]) * unit_CD
        cm_S[:,j] = abs(float(CM_S[j])) * unit_CM
        
        #################################
        ############ POINT A ############
        #################################
        CL_A[j] = CL_wb_fromStoA[-1,j]
        CD_A[j] = CD_fromStoA[-1,j]
        CM_A[j] = CM_CG_fromStoA[-1,j]
        
        cl_A[:,j] = float(CL_A[j]) * unit_CL
        cd_A[:,j] = float(CD_A[j]) * unit_CD
        cm_A[:,j] = abs(float(CM_A[j])) * unit_CM
        
        #################################
        ############ POINT C ############
        #################################
        CL_C[j] = CL_wb_fromGust1toC[-1,j]
        CD_C[j] = CD_fromGust1toC[-1,j]
        CM_C[j] = CM_CG_fromGust1toC[-1,j]
        
        cl_C[:,j] = float(CL_C[j]) * unit_CL
        cd_C[:,j] = float(CD_C[j]) * unit_CD
        cm_C[:,j] = abs(float(CM_C[j])) * unit_CM
        
        #################################
        ############ POINT D ############
        #################################
        CL_D[j] = CL_wb_fromGust2toD[-1,j]
        CD_D[j] = CD_fromGust2toD[-1,j]
        CM_D[j] = CM_CG_fromGust2toD[-1,j]
        
        cl_D[:,j] = float(CL_D[j]) * unit_CL
        cd_D[:,j] = float(CD_D[j]) * unit_CD
        cm_D[:,j] = abs(float(CM_D[j])) * unit_CM

# INVERTED FLIGHT INITIALIZATION

#############################
###### INVERTED FLIGHT ######
#############################
if ( flight_env_case_neg == "Case1_inverted" ): 
    
    for j in range(0,n_Mass):
        ##########################################
        ############ POINT S INVERTED ############
        ##########################################
        CL_Sinv[j] = CL_wb_from0toSinv[-1,j]
        CD_Sinv[j] = CD_from0toSinv[-1,j]
        CM_Sinv[j] = CM_CG_from0toSinv[-1,j]
        
        cl_Sinv[:,j] = abs(float(CL_Sinv[j])) * unit_CL
        cd_Sinv[:,j] = abs(float(CD_Sinv[j])) * unit_CD
        cm_Sinv[:,j] = abs(float(CM_Sinv[j])) * unit_CM
    
        #################################
        ############ POINT G ############
        #################################
        CL_G[j] = CL_wb_fromSinvtoG[-1,j]
        CD_G[j] = CD_fromSinvtoG[-1,j]
        CM_G[j] = CM_CG_fromSinvtoG[-1,j]
        
        cl_G[:,j] = abs(float(CL_G[j])) * unit_CL
        cd_G[:,j] = abs(float(CD_G[j])) * unit_CD
        cm_G[:,j] = abs(float(CM_G[j])) * unit_CM
    
        #################################
        ############ POINT F ############
        #################################
        CL_F[j] = CL_wb_fromGust1toF[-1,j]
        CD_F[j] = CD_fromGust1toF[-1,j]
        CM_F[j] = CM_CG_fromGust1toF[-1,j]
        
        cl_F[:,j] = abs(float(CL_F[j])) * unit_CL
        cd_F[:,j] = abs(float(CD_F[j])) * unit_CD
        cm_F[:,j] = abs(float(CM_F[j])) * unit_CM
    
        #################################
        ############ POINT E ############
        #################################
        CL_E[j] = CL_wb_fromFtoE[-1,j]
        CD_E[j] = CD_fromFtoE[-1,j]
        CM_E[j] = CM_CG_fromFtoE[-1,j]
        
        cl_E[:,j] = abs(float(CL_E[j])) * unit_CL
        cd_E[:,j] = abs(float(CD_E[j])) * unit_CD
        cm_E[:,j] = abs(float(CM_E[j])) * unit_CM
    
elif ( flight_env_case_neg == "Case2_inverted" ):
    
    for j in range(0,n_Mass):
        ##########################################
        ############ POINT S INVERTED ############
        ##########################################
        CL_Sinv[j] = CL_wb_from0toSinv[-1,j]
        CD_Sinv[j] = CD_from0toSinv[-1,j]
        CM_Sinv[j] = CM_CG_from0toSinv[-1,j]
        
        cl_Sinv[:,j] = abs(float(CL_Sinv[j])) * unit_CL
        cd_Sinv[:,j] = abs(float(CD_Sinv[j])) * unit_CD
        cm_Sinv[:,j] = abs(float(CM_Sinv[j])) * unit_CM
        
        #################################
        ############ POINT G ############
        #################################
        CL_G[j] = CL_wb_fromSinvtoG[-1,j]
        CD_G[j] = CD_fromSinvtoG[-1,j]
        CM_G[j] = CM_CG_fromSinvtoG[-1,j]
        
        cl_G[:,j] = abs(float(CL_G[j])) * unit_CL
        cd_G[:,j] = abs(float(CD_G[j])) * unit_CD
        cm_G[:,j] = abs(float(CM_G[j])) * unit_CM
        
        #################################
        ############ POINT F ############
        #################################
        CL_F[j] = CL_wb_fromGtoF[-1,j]
        CD_F[j] = CD_fromGtoF[-1,j]
        CM_F[j] = CM_CG_fromGtoF[-1,j]
        
        cl_F[:,j] = abs(float(CL_F[j])) * unit_CL
        cd_F[:,j] = abs(float(CD_F[j])) * unit_CD
        cm_F[:,j] = abs(float(CM_F[j])) * unit_CM
        
        #################################
        ############ POINT E ############
        #################################
        CL_E[j] = CL_wb_fromFtoE[-1,j]
        CD_E[j] = CD_fromFtoE[-1,j]
        CM_E[j] = CM_CG_fromFtoE[-1,j]
        
        cl_E[:,j] = abs(float(CL_E[j])) * unit_CL
        cd_E[:,j] = abs(float(CD_E[j])) * unit_CD
        cm_E[:,j] = abs(float(CM_E[j])) * unit_CM
        
elif ( flight_env_case_neg == "Case3_inverted" ): 
    
    for j in range(0,n_Mass):
    
        ##########################################
        ############ POINT S INVERTED ############
        ##########################################
        CL_Sinv[j] = CL_wb_from0toSinv[-1,j]
        CD_Sinv[j] = CD_from0toSinv[-1,j]
        CM_Sinv[j] = CM_CG_from0toSinv[-1,j]
        
        cl_Sinv[:,j] = abs(float(CL_Sinv[j])) * unit_CL
        cd_Sinv[:,j] = abs(float(CD_Sinv[j])) * unit_CD
        cm_Sinv[:,j] = abs(float(CM_Sinv[j])) * unit_CM
        
        #################################
        ############ POINT G ############
        #################################
        CL_G[j] = CL_wb_fromSinvtoG[-1,j]
        CD_G[j] = CD_fromSinvtoG[-1,j]
        CM_G[j] = CM_CG_fromSinvtoG[-1,j]
        
        cl_G[:,j] = abs(float(CL_G[j])) * unit_CL
        cd_G[:,j] = abs(float(CD_G[j])) * unit_CD
        cm_G[:,j] = abs(float(CM_G[j])) * unit_CM
        
        #################################
        ############ POINT F ############
        #################################
        CL_F[j] = CL_wb_fromGtoF[-1,j]
        CD_F[j] = CD_fromGtoF[-1,j]
        CM_F[j] = CM_CG_fromGtoF[-1,j]
        
        cl_F[:,j] = abs(float(CL_F[j])) * unit_CL
        cd_F[:,j] = abs(float(CD_F[j])) * unit_CD
        cm_F[:,j] = abs(float(CM_F[j])) * unit_CM
        
        #################################
        ############ POINT E ############
        #################################
        CL_E[j] = CL_wb_fromFtoE[-1,j]
        CD_E[j] = CD_fromFtoE[-1,j]
        CM_E[j] = CM_CG_fromFtoE[-1,j]
        
        cl_E[:,j] = abs(float(CL_E[j])) * unit_CL
        cd_E[:,j] = abs(float(CD_E[j])) * unit_CD
        cm_E[:,j] = abs(float(CM_E[j])) * unit_CM

n = 1
for j in range(0, n_Mass):
    cl_comparison, cd_comparison, cm_comparison =  plot_distributions(y_halfspan[:],\
                           unit_CL[:], unit_CD[:], unit_CM[:],\
                           cl_S[:,j], cl_A[:,j], cl_C[:,j], cl_D[:,j], cl_Sinv[:,j], cl_G[:,j], cl_F[:,j], cl_E[:,j],\
                           cd_S[:,j], cd_A[:,j], cd_C[:,j], cd_D[:,j], cd_Sinv[:,j], cd_G[:,j], cd_F[:,j], cd_E[:,j],\
                           cm_S[:,j], cm_A[:,j], cm_C[:,j], cm_D[:,j], cm_Sinv[:,j], cm_G[:,j], cm_F[:,j], cm_E[:,j], n) 
    n = n + 1
    
# =============================================================================    
# MOVING FIGURES INSIDE OUTPUT FOLDER
# FinalEnvelope
n = 1
for i in range(0, n_Mass):
    
    figure_name = r'\cl_comparison' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\AeroCoefficientDistributions'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)
    
    figure_name = r'\cd_comparison' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\AeroCoefficientDistributions'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)
    
    figure_name = r'\cm_comparison' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\AeroCoefficientDistributions'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)

    n = n + 1
# ============================================================================= 
plt.close("all")