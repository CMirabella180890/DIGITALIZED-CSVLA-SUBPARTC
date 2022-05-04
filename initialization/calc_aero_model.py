# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 20:33:22 2022

Script file to interpolate wing-body aerodynamics. This is done through simple 
interpolation techniques.  Care must be  exercised when the boundaries of  the 
interpolation process  are defined by the user.  Check also other  aerodynamic 
parameters.

@author: claum
"""
import sys
import numpy as np
from numpy.polynomial.chebyshev import chebfit
from numpy.polynomial.chebyshev import chebval
import matplotlib.pyplot as plt 
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# ============================================================================

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++ FUNCTION TO SAVE DATA VALUABLE FOR THE INTERPOLATION PROCESS +++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
def initialization(aircraft_data):
    
    """
    This function initialize data for the interpolation process, extracting 
    them from the conveniently defined SimpleNamespace object 'aircraft_data'
    """

    interpolation_flag = aircraft_data.Wing_body_aerodynamic["interpolation_flag"]["Value"]
    
    # =============================================================================
    # SWITCH CASE TO INTERPOLATION 
    # =============================================================================
    if interpolation_flag == "User":
        # ++++++++++++++++++++++++++++++++++++++++++++++++
        # ++++++++++++++++ Initialization ++++++++++++++++
        # ++++++++++++++++++++++++++++++++++++++++++++++++
        CL_max_clean   = aircraft_data.Wing_body_aerodynamic["MaX_CL_clean"]["Value"]    # Maximum lift coefficient in clean config.
        CL_max_takeoff = aircraft_data.Wing_body_aerodynamic["Max_CL_takeoff"]["Value"]  # Maximum lift coefficient in takeoff config.
        CL_max_landing = aircraft_data.Wing_body_aerodynamic["Max_CL_landing"]["Value"]  # Maximum lift coefficient in landing config. 
        CL_inverted    = aircraft_data.Wing_body_aerodynamic["Inverted_CL"]["Value"]     # Maximum lift coefficient during inverted flight
        CL_star        = aircraft_data.Wing_body_aerodynamic["CL_star"]["Value"]         # End of the linear part of the lift coefficient curve

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ++++++++++++++++ Interpolation data ++++++++++++++++
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        wb_aoa         = aircraft_data.Wing_body_aerodynamic["wb_angle_of_attack"]["Value"]
        wb_CL          = aircraft_data.Wing_body_aerodynamic["wb_CL"]["Value"]
        wb_CD          = aircraft_data.Wing_body_aerodynamic["wb_CD"]["Value"]
        wb_CM          = aircraft_data.Wing_body_aerodynamic["wb_CM"]["Value"]
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Numpy arrays to store  conveniently aerodynamic data
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        n_Interp       = len(wb_CL)
        CL_wb          = np.ones((n_Interp,1))
        CD_wb          = np.ones((n_Interp,1))
        CM_wb          = np.ones((n_Interp,1))
        AOA_wb         = np.ones((n_Interp,1))
        for i in range(0,n_Interp): 
            CL_wb[i]  = wb_CL[i]
            CD_wb[i]  = wb_CD[i]
            CM_wb[i]  = wb_CM[i]
            AOA_wb[i] = wb_aoa[i]
            
    elif interpolation_flag == "Open VSP":
        # ++++++++++++++++++++++++++++++++++++++++++++++++
        # ++++++++++++++++ Initialization ++++++++++++++++
        # ++++++++++++++++++++++++++++++++++++++++++++++++
        CL_max_clean   = aircraft_data.Wing_body_aerodynamic["MaX_CL_clean"]["Value"]    # Maximum lift coefficient in clean config.
        CL_max_takeoff = aircraft_data.Wing_body_aerodynamic["Max_CL_takeoff"]["Value"]  # Maximum lift coefficient in takeoff config.
        CL_max_landing = aircraft_data.Wing_body_aerodynamic["Max_CL_landing"]["Value"]  # Maximum lift coefficient in landing config. 
        CL_inverted    = aircraft_data.Wing_body_aerodynamic["Inverted_CL"]["Value"]     # Maximum lift coefficient during inverted flight
        CL_star        = aircraft_data.Wing_body_aerodynamic["CL_star"]["Value"]         # End of the linear part of the lift coefficient curve

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ++++++++++++++++ Interpolation data ++++++++++++++++
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        wb_aoa         = aircraft_data.Wing_body_aerodynamic["wb_angle_of_attack"]["Value"]
        wb_CL          = aircraft_data.Wing_body_aerodynamic["wb_CL"]["Value"]
        wb_CD          = aircraft_data.Wing_body_aerodynamic["wb_CD"]["Value"]
        wb_CM          = aircraft_data.Wing_body_aerodynamic["wb_CM"]["Value"]
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Numpy arrays to store  conveniently aerodynamic data
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        n_Interp       = len(wb_CL)
        CL_wb          = np.ones((n_Interp,1))
        CD_wb          = np.ones((n_Interp,1))
        CM_wb          = np.ones((n_Interp,1))
        AOA_wb         = np.ones((n_Interp,1))
        for i in wb_aoa: 
            CL_wb[i]  = wb_CL[i]
            CD_wb[i]  = wb_CD[i]
            CM_wb[i]  = wb_CM[i]
            AOA_wb[i] = wb_aoa[i]
    elif interpolation_flag == "Other":
        # ++++++++++++++++++++++++++++++++++++++++++++++++
        # ++++++++++++++++ Initialization ++++++++++++++++
        # ++++++++++++++++++++++++++++++++++++++++++++++++
        CL_max_clean   = aircraft_data.Wing_body_aerodynamic["MaX_CL_clean"]["Value"]    # Maximum lift coefficient in clean config.
        CL_max_takeoff = aircraft_data.Wing_body_aerodynamic["Max_CL_takeoff"]["Value"]  # Maximum lift coefficient in takeoff config.
        CL_max_landing = aircraft_data.Wing_body_aerodynamic["Max_CL_landing"]["Value"]  # Maximum lift coefficient in landing config. 
        CL_inverted    = aircraft_data.Wing_body_aerodynamic["Inverted_CL"]["Value"]     # Maximum lift coefficient during inverted flight
        CL_star        = aircraft_data.Wing_body_aerodynamic["CL_star"]["Value"]         # End of the linear part of the lift coefficient curve

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ++++++++++++++++ Interpolation data ++++++++++++++++
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        wb_aoa         = aircraft_data.Wing_body_aerodynamic["wb_angle_of_attack"]["Value"]
        wb_CL          = aircraft_data.Wing_body_aerodynamic["wb_CL"]["Value"]
        wb_CD          = aircraft_data.Wing_body_aerodynamic["wb_CD"]["Value"]
        wb_CM          = aircraft_data.Wing_body_aerodynamic["wb_CM"]["Value"]
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Numpy arrays to store  conveniently aerodynamic data
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++
        n_Interp       = len(wb_CL)
        CL_wb          = np.ones((n_Interp,1))
        CD_wb          = np.ones((n_Interp,1))
        CM_wb          = np.ones((n_Interp,1))
        AOA_wb         = np.ones((n_Interp,1))
        for i in wb_CL: 
            CL_wb[i]  = wb_CL[i]
            CD_wb[i]  = wb_CD[i]
            CM_wb[i]  = wb_CM[i]
            AOA_wb[i] = wb_aoa[i]
        
    return CL_max_clean, CL_max_takeoff, CL_max_landing, CL_inverted, CL_star, CL_wb, CD_wb, CM_wb, AOA_wb

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++ INTERPOLATION OF THE LIFT COEFFICIENT CURVE ++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
def cl_wb_interpolation(n_Points, n_Poly1, n_Poly2, CL_wb, CL_star, AOA_wb):
    
    tol1   = 1E-1
    tol2   = 1E-2
    check1 = np.abs(CL_wb - CL_star)<tol1
    check2 = np.abs(CL_wb - CL_star)<tol2
    if (any(check2)==True):
        index         = np.argwhere(np.abs(CL_wb - CL_star)<tol2)
        index_cl_star = index[0,0]
    elif (any(check1)==True) and (all(check2)==False):
        index         = np.argwhere(np.abs(CL_wb - CL_star)<tol1)
        index_cl_star = index[0,0]
        
    alpha_interp = np.linspace(AOA_wb[0], AOA_wb[-1], n_Points)
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # +++++++++++++++ ALPHA STAR INTERPOLATION +++++++++++++++ 
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    cl_alpha_star  = CL_wb[index_cl_star-1:index_cl_star+1]
    aoa_alpha_star = AOA_wb[index_cl_star-1:index_cl_star+1]
    p_alpha_star   = chebfit(cl_alpha_star[:, 0], aoa_alpha_star[:, 0], deg = n_Poly1)
    ALPHA_star     = chebval(CL_star, p_alpha_star)
    
    # LINEAR PART 
    AOA_wb1        = np.append(AOA_wb[0:index_cl_star-1], ALPHA_star)
    CL_wb1         = np.append(CL_wb[0:index_cl_star-1], CL_star)
    p_cl_wb1       = chebfit(AOA_wb1, CL_wb1, n_Poly1)
    CL_model1      = chebval(alpha_interp[alpha_interp<ALPHA_star], p_cl_wb1)
    
    # NON LINEAR, POST STALL PARABOLIC FIT
    AOA_wb2        = np.append(ALPHA_star, AOA_wb[index_cl_star+1:-1])
    CL_wb2         = np.append(CL_star, CL_wb[index_cl_star+1:-1])
    p_cl_wb2       = chebfit(AOA_wb2, CL_wb2, n_Poly2)
    CL_model2      = chebval(alpha_interp[alpha_interp>ALPHA_star], p_cl_wb2)
    
    # ASSEMBLING CURVES 
    CL_fullmodel   = np.append(CL_model1, CL_model2)   
    
    return index_cl_star, p_alpha_star, ALPHA_star, alpha_interp, CL_fullmodel, p_cl_wb1, p_cl_wb2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++ INTERPOLATION OF THE DRAG COEFFICIENT CURVE ++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
def cd_wb_interpolation(alpha_interp, alfa0L, n_Poly2, CD_wb, AOA_wb):
    
    p_cd_wb      = chebfit(AOA_wb[:,0], CD_wb[:,0],  n_Poly2)
    CD_fullmodel = chebval(alpha_interp, p_cd_wb)
    CD0          = chebval(alfa0L, p_cd_wb)
    
    return p_cd_wb, CD_fullmodel, CD0
    
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++ INTERPOLATION OF THE PITCHING MOMENT COEFFICIENT CURVE ++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def cm_wb_interpolation(n_Points, n_Poly1, CL_wb, CM_wb, CL_star, CL_fullmodel, index_cl_star):
    
    # INITIALIZATION
    tol1      = 1E-1
    tol2      = 1E-2
    check1    = np.abs(CL_wb - CL_star)<tol1
    check2    = np.abs(CL_wb - CL_star)<tol2
    CL_for_CM = CL_wb[CL_wb < CL_star]
    length_CL = len(CL_for_CM)
    
    # ESTABLISHING THE CORRECT INDEX - END OF THE LINEAR PART
    if (any(check2)==True):
        index         = np.argwhere(np.abs(CL_wb - CL_star)<tol2)
        index_cl_star = index[0,0]
    elif (any(check1)==True) and (all(check2)==False):
        index         = np.argwhere(np.abs(CL_wb - CL_star)<tol1)
        index_cl_star = index[0,0]
        
    if (length_CL != index_cl_star):
        sys.exit("ERROR: Uncorrect pitching moment interpolation! Review data format.")
        
    # ACTUAL INTERPOLATION 
    CM_for_CM    = CM_wb[0:index_cl_star,0]
    p_CM_wb      = chebfit(CL_for_CM, CM_for_CM, deg = n_Poly1)
    CM_fullmodel = chebval(CL_fullmodel, p_CM_wb)
    
    return p_CM_wb, CM_fullmodel    

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++ PLOTTING THE LIFT COEFFICIENT CURVE FITTING ++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
def plot_interpolation_model(AOA_wb, CL_wb, CD_wb, CM_wb, AOA_fullmodel, CL_fullmodel, CD_fullmodel, CM_fullmodel):
    
    # =======================================
    # ======= CL INTERPOLATION FIGURE =======
    # =======================================
    fig1 = plt.figure()
    plt.scatter(AOA_wb, CL_wb, marker=".", s=5, alpha=0.85, color='black')
    plt.plot(AOA_fullmodel, CL_fullmodel)
    plt.xlabel(r'Angle of attack - $\alpha_{wb}$ [deg]')    # x-label to the axes.
    plt.ylabel(r'Wing body Lift coefficient - $C_{L, wb}$') # y-label to the axes.
    plt.title(r'Lift coefficient curve interpolation')      # Title to the axes.
    plt.minorticks_on()
    plt.grid(True, linestyle='-.', which="both")
    plt.savefig('figure1.pdf', bbox_inches='tight')
    plt.show()
    # =======================================
    
    # =======================================
    # ======= CD INTERPOLATION FIGURE =======
    # =======================================
    fig2 = plt.figure()
    plt.scatter(AOA_wb, CD_wb, marker=".", s=5, alpha=0.85, color='black')
    plt.plot(AOA_fullmodel, CD_fullmodel)
    plt.xlabel(r'Angle of attack - $\alpha_{wb}$ [deg]')    # x-label to the axes.
    plt.ylabel(r'Wing body Drag coefficient - $C_{D, wb}$') # y-label to the axes.
    plt.title(r'Drag coefficient curve interpolation')      # Title to the axes.
    plt.minorticks_on()
    plt.grid(True, linestyle='-.', which="both")
    plt.savefig('figure2.pdf', bbox_inches='tight')
    plt.show()
    # =======================================
    
    # =======================================
    # ======= CM INTERPOLATION FIGURE =======
    # =======================================
    fig3 = plt.figure()
    plt.scatter(CL_wb, CM_wb, marker=".", s=5, alpha=0.85, color='black')
    plt.plot(CL_fullmodel, CM_fullmodel)
    plt.xlabel(r'Wing body Lift coefficient - $C_{L, wb}$')            # x-label to the axes.
    plt.ylabel(r'Wing body Pitching moment coefficient - $C_{M, wb}$') # y-label to the axes.
    plt.title(r'Pitching mom. coefficient curve interpolation')        # Title to the axes.
    plt.minorticks_on()
    plt.grid(True, linestyle='-.', which="both")
    plt.savefig('figure3.pdf', bbox_inches='tight')
    plt.show()
    # =======================================
    
    return fig1, fig2, fig3 
    