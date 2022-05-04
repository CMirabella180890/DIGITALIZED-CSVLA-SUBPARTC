# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 09:08:33 2022

@author: claum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CS-VLA 333 Flight envelope                                                 %
%  (a) GENERAL. Compliance with the strength requirements of this subpart     %
%      must be shown at any combination of airspeed and load factor on and    %
%      within the boundaries of a flight envelope (similar to the one in      %
%      sub-paragraph(d) of this paragraph) that represents the envelope of    %
%      the flight loading conditions specified by the manoeuvring and gust    %
%      criteria of sub-paragraphs(b) and (c) of this paragraph                %
%      respectively.                                                          %
%                                                                             %
%  (b) MANOEUVRING ENVELOPE. Except where limited by maximum static lift      %
%      coefficients, the aeroplane is assumed to be subjected to              %
%      symmetrical manoeuvres resulting in the following limit load           % 
%      factors:                                                               %
%      (1) the positive manoeuvring load factor specified in CS-VLA 337 at    %
%          speeds up to VD;                                                   %
%      (2) the negative manoeuvring load factor specified in CS-VLA 337 at    %
%          VC;                                                                % 
%      (3) factors varying linearly with speed from the specified value at    %
%          VC to 0.0 at VD.                                                   %
%                                                                             %
%  (c) GUST ENVELOPE.                                                         %
%      (1) The aeroplane is assumed to be subjected to symmetrical vertical   %
%          gusts in level flight. The resulting limit load factors must       %
%          correspond to the conditions determined as follows:                %
%          (i) positive (up) and negative (down) gusts of 15.24 m/s at VC     %
%              must be considered;                                            %
%         (ii) positive and negative gusts of 7.62 m/s at VD must be          %
%              considered.                                                    %
%      (2) The following assumptions must be made:                            %
%          (i) the shape of the gust is                                       %
%                  Ude   /          2*pi*S   \                                %
%              U = --- * | 1 - cos( ------ ) |                                %
%                   2    \          25*MGC   /                                %
%              where                                                          %
%              S   = distance penetrated into gust (m);                       %
%              MGC = mean geometric chord of wing (m);                        %
%              Ude = derived gust velocity referred to in                     %
%                    sub-paragraph(c)(1) (m/s)                                %
%         (ii) gust load factors vary linearly with speed between VC and      %  
%              VD.                                                            %
%                                                                             %
%  (d) FLIGHT ENVELOPE. Figure inside the CS-VLA airworthiness                %
%      requirements. NOTE: point G need not be investigated when the          %
%      supplementary condition specified in CS-VLA 369 is investigated.       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

"""
import shutil
import matplotlib.pyplot as plt
from types import SimpleNamespace
import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval
import csvla
obj = csvla.csvla

# AIRCRAFT INPUTS 
N            = 1000 
x_ac         = aircraft_data.Aerodynamic_centre["wing_x_ac"]["Value"]
max_fwd_ac   = aircraft_data.Aerodynamic_centre["max_fwd_ac"]["Value"]
max_aft_ac   = aircraft_data.Aerodynamic_centre["max_aft_ac"]["Value"]
max_up_bcg   = aircraft_data.Aerodynamic_centre["max_up_bcg"]["Value"]
max_down_bcg = aircraft_data.Aerodynamic_centre["max_down_bcg"]["Value"]
wing_xle     = aircraft_data.Wing["wing_xle"]["Value"]
MAC          = aircraft_data.Wing["MAC"]["Value"]
S_wing       = aircraft_data.Wing["S"]["Value"]
nmax         = aircraft_data.Aircraft_inputs["nmax"]["Value"]
nmin         = aircraft_data.Aircraft_inputs["nmin"]["Value"]
rho0         = aircraft_data.Standard_atmosphere["sea_level"]["Density"]["Value"]
rho_op       = aircraft_data.Standard_atmosphere["operative_ceiling"]["Density"]["Value"]
rho_th       = aircraft_data.Standard_atmosphere["theoretical_ceiling"]["Density"]["Value"]
g            = aircraft_data.Constants["g"]["Value"]
MTOM         = aircraft_data.Aircraft_mass["Max_takeoff_mass"]["Value"]
EM           = aircraft_data.Aircraft_mass["Empty_mass"]["Value"]
UM           = aircraft_data.Aircraft_mass["Useful_mass"]["Value"]
CM           = aircraft_data.Aircraft_mass["Crew_mass"]["Value"]
FM           = aircraft_data.Aircraft_mass["Fuel_mass"]["Value"]
OM           = aircraft_data.Aircraft_mass["Oil_mass"]["Value"]
Max_fwd_cg   = aircraft_data.Aircraft_mass["Max_fwd_cg"]["Value"]
Max_aft_cg   = aircraft_data.Aircraft_mass["Max_aft_cg"]["Value"] 
l_ht         = aircraft_data.Horizontal_Tailplane["ht_moment_arm"]["Value"]
x_vt         = aircraft_data.Vertical_Tailplane["vt_xloc"]["Value"]
CL_max_clean = aircraft_data.Wing_body_aerodynamic["MaX_CL_clean"]["Value"]
CL_inverted  = aircraft_data.Wing_body_aerodynamic["Inverted_CL"]["Value"]

# MASS ENVELOPE 
Minimum_takeoff_mass = EM + UM - (OM + 0.8 * FM)
n_Mass               = 4
Mass_envelope        = np.linspace(Minimum_takeoff_mass, MTOM, n_Mass)

# CENTRE OF GRAVITY ENVELOPE
xcg_envelope         = np.linspace(Max_fwd_cg, Max_aft_cg, n_Mass)

# AERODYNAMIC CENTRE ENVELOPE 
x_ac_envelope        = np.linspace(max_fwd_ac, max_aft_ac, n_Mass)

# HORIZONTAL TAILPLANE MOMENT ARM ENVELOPE
l_ht_envelope        = np.ones((n_Mass,1))
for i in range(0, n_Mass): 
    if (xcg_envelope[i] < x_ac_envelope[i]):
        l_ht_envelope[i] = l_ht + xcg_envelope[i] * MAC
    elif (xcg_envelope[i] > x_ac_envelope[i]):
        l_ht_envelope[i] = l_ht - xcg_envelope[i] * MAC
        
# VERTICAL TAILPLANE MOMENT ARM ENVELOPE
l_vt                 = x_vt - wing_xle - x_ac_envelope * MAC
l_vt_envelope        = np.ones((n_Mass,1))
for i in range(0, n_Mass): 
    if (xcg_envelope[i] < x_ac_envelope[i]):
        l_vt_envelope[i] = l_vt[i] + xcg_envelope[i] * MAC
    elif (xcg_envelope[i] > x_ac_envelope[i]):
        l_vt_envelope[i] = l_vt[i] - xcg_envelope[i] * MAC

# THRUST LINE DATA
max_up_thrust_line   = aircraft_data.Thrust_characteristics["max_up_thrust_line"]["Value"]
max_down_thrust_line = aircraft_data.Thrust_characteristics["max_down_thrust_line"]["Value"]
thrust_line_envelope = np.linspace(max_down_thrust_line, max_up_thrust_line, n_Mass) 

# BCG ENVELOPE 
bcg_envelope = np.linspace(max_down_bcg, max_up_bcg, n_Mass)

# AIRCRAFT WING LOADING IN [Pa] AND IN [kg/m^2]
WS = np.ones((n_Mass,1))
MS = np.ones((n_Mass,1))
for i in range(0, n_Mass):
    WS[i] = ( Mass_envelope[i] * g ) / ( S_wing )
    MS[i] = ( Mass_envelope[i] ) / ( S_wing )
    
# ------ ################################################### ------------------
# ------ ##### UPDATINGE THE SIMPLE NAMESPACE VARIABLE ##### ------------------
# ------ ################################################### ------------------
New_values = {
	"Aircraft_mass": {
	    "Max_takeoff_mass": {"Value": 450,  "Unit": "kg"},
	    "Empty_mass"      : {"Value": 281,  "Unit": "kg"},
	    "Useful_mass"     : {"Value": 169,  "Unit": "kg"},
	    "Fuel_mass"       : {"Value": 50.5, "Unit": "kg"},
	    "Crew_mass"       : {"Value": 180,  "Unit": "kg"},
		"Oil_mass"        : {"Value": 1.85, "Unit": "kg"},
		"Mass_values_num" : {"Value": 4,    "Unit": "Pure number"},
		"Max_fwd_cg"      : {"Value": 0.21, "Unit": "MAC percentage"},
		"Max_aft_cg"      : {"Value": 0.24, "Unit": "MAC percentage"},
        "xcg_envelope"    : {"Value": xcg_envelope, "Unit": "MAC percentage"},
		"IY_inertia_mom"  : {"Value": 513,  "Unit": "kg * m^2"},
        "Min_takeoff_mass": {"Value": Minimum_takeoff_mass, "Unit": "kg"},
        "n_Mass"          : {"Value": n_Mass, "Unit": "Pure number"},
        "Mass_envelope"   : {"Value": Mass_envelope, "Unit": "kg"}
	    },
	"Aerodynamic_centre": {
		"wing_x_ac"                    : {"Value": 0.25, "Unit": "MAC percentage"},
		"max_up_bcg"                   : {"Value": 0.00, "Unit": "m"},
		"max_down_bcg"                 : {"Value": 0.00, "Unit": "m"},
        "bcg_envelope"                 : {"Value": bcg_envelope, "Unit":"m"},
		"lift_horizontal_component_arm": {"Value": 0.00, "Unit": "m"},
		"max_fwd_ac"                   : {"Value": 0.25, "Unit": "MAC percentage"},
		"max_aft_ac"                   : {"Value": 0.25, "Unit": "MAC percentage"},
        "x_ac_envelope"                : {"Value": x_ac_envelope, "Unit": "MAC percentage"}
	    },
	"Horizontal_tail_loads" : {
		"time_step"      : {"Value": 0.006,       "Unit": "1/s"},
		"time_interval"  : {"Value": 50.0,        "Unit": "s"},
		"damping_factor" : {"Value": 0.30,        "Unit": "Non dimensional"},
		"type_flag"      : {"Value": "Aerobatic", "Unit": "flag"},
		"command_flag"   : {"Value": "Stick",     "Unit": "flag"},
        "l_ht_envelope"  : {"Value": l_ht_envelope, "Unit": "m"}
		},
	"Vertical_tail_loads" : {
		"CY0"                              : {"Value": 0.000036156860, "Unit": "Non dimensional"},
		"CYdr"                             : {"Value": 0.000644,       "Unit": "1/deg"},
		"CY_case_a2"                       : {"Value": 0.0245,         "Unit": "Non dimensional"},
		"CY_case_a3"                       : {"Value": 0.0233,         "Unit": "Non dimensional"},
		"radius_of_gyration_vertical_tail" : {"Value": 0.30,           "Unit": "m"},
		"max_yaw_angle"                    : {"Value": 15.0,           "Unit": "deg"},
        "l_vt_envelope"                    : {"Value": l_vt_envelope, "Unit": "m"}
	    },
	"Thrust_characteristics": {
		"thrust_line"          : {"Value": 0.0, "Unit": "m"},
		"max_up_thrust_line"   : {"Value": 0.0, "Unit": "m"},
		"max_down_thrust_line" : {"Value": 0.0, "Unit": "m"},
        "thrust_line_envelope" : {"Value": thrust_line_envelope, "Unit": "m"}
		},
    "Wing_loading": {
        "Wing_loading_Pa" : {"Value": WS, "Unit":"Pa"},
        "Wing_loading_kg" : {"Value": MS, "Unit": "kg/m^2"}
        }
    }
# UPDATING DICTIONARY
aircraft.update(New_values)
 # UPDATING THE SIMPLE NAMESPACE OBJECT
aircraft_data = SimpleNamespace(**aircraft)  
    
################################################
####### LOAD FACTORS VECTOR CALCULATIONS #######
################################################
npos = obj.calcn(nmax, N)
nneg = obj.calcn(nmin, N)

# MINIMUM CRUISE DESIGN AIRSPEED 
min_cruise_design_airspeed = aircraft_data.Design_airspeed["Minimum_design_cruise_airspeed"]["Value"]
    
################################################
########### STALL SPEED CALCULATIONS ###########
################################################
# x = calcvs(obj, rho, WingLoading, MaxLiftCoeff, PositiveLoadFactors)
# This function defines a vector with stall airspeed for the chosen
# aircraft, within the precribed range of load factors. Check the
# class file csvla.m to have a complete documentation.
# Positive stall speed 
VSpos = np.ones((N, n_Mass))
VSneg = np.ones((N, n_Mass))
for j in range(0, n_Mass):
    for i in range(0, N):
        VSpos[i,j] = obj.calcvs(rho0, WS[j], CL_max_clean, npos[i])
        VSneg[i,j] = obj.calcvs(rho0, WS[j], CL_inverted, nneg[i])
 
#######################################################
########### CALCULATION OF THE CRUISE SPEED ###########
#######################################################
# x = calcvc(obj, WingLoading, MaxContinuousPowerSpeedVH)
#  This function identifies (following CS-VLA airworthiness reg.)
#  maximum cruise speed (Point C) for flight envelope calculations. 
#  To have a complete documentation check the class file csvla.m
#  VH design speed for max continous power: this airspeed is not available
#  but must be known. From CS - VLA Airworthiness rules
nC = np.full(n_Mass, nmax)
cruise_speed = np.ones((n_Mass,1))
for i in range(0,n_Mass):
    cruise_speed[i] = obj.calcvc(WS[i], 0.0)
    
VF = cruise_speed    
#######################################################
############ CALCULATION OF THE DIVE SPEED ############
#######################################################
#  x = calcvd(obj, MinDesignCruiseSpeed, CruiseSpeedVC)
#  This function identifies (following CS-VLA airworthiness reg.)
#  the maximum dive speed (Point D) for flight envelope
#  calculations. To have a complete documentation check the class
#  file csvla.m
VD = np.ones((n_Mass,1))
VE = np.ones((n_Mass,1))
for i in range(0,n_Mass):
    VD[i] = obj.calcvd(min_cruise_design_airspeed, cruise_speed[i])
    
VE = VD
# CORRESPONDING LOAD FACTORS
nD = np.full(n_Mass, nmax) 
nE = np.full(n_Mass, nmin)    

######################################################################
############ FLIGHT ENVELOPE POINTS - POINT 0 AND POINT S ############
######################################################################
n0 = 0.0
nS = np.full(n_Mass, 1.0) 
VS = np.ones((n_Mass, 1))
for i in range(0,n_Mass):
    VS[i] = obj.calcvs(rho0, WS[i], CL_max_clean, nS[i])
    
V0 = VS

# FLIGHT ENVELOPE - FROM 0 TO S 
n_from0toS = np.zeros((N,n_Mass))
V_from0toS = np.zeros((N,n_Mass))
for i in range(0, n_Mass):
    n_from0toS[:,i] = np.linspace(0.0, nS[i], N)
    V_from0toS[:,i] = np.full(N, VS[i])
    
###############################################################################
############ FLIGHT ENVELOPE POINTS - POINT 0 AND POINT S INVERTED ############
###############################################################################
n0inv = 0.0 
nSinv = np.full(n_Mass, -1.0) 
VSinv  = np.ones((n_Mass,1))
for i in range(0,n_Mass):
    VSinv[i] = obj.calcvs(rho0, WS[i], CL_inverted, nSinv[i])
    
V0inv = VSinv 

# FLIGHT ENVELOPE - FROM 0 TO S 
n_from0toSinv = np.zeros((N,n_Mass))
V_from0toSinv = np.zeros((N,n_Mass))
for i in range(0, n_Mass):
    n_from0toSinv[:,i] = np.linspace(0.0, nSinv[i], N)
    V_from0toSinv[:,i] = np.full(N, VSinv[i])
        
##########################################################
############ FLIGHT ENVELOPE POINTS - POINT A ############
##########################################################
nA = np.full(n_Mass, nmax) 
VA = np.ones((n_Mass,1))
for i in range(0,n_Mass):
    VA[i] = obj.calcvs(rho0, WS[i], CL_max_clean, nA[i])
    
# FLIGHT ENVELOPE - FROM S TO A 
n_fromStoA = np.zeros((N,n_Mass))
V_fromStoA = np.zeros((N,n_Mass))
for i in range(0, n_Mass):
    n_fromStoA[:,i] = np.linspace(nS[i], nA[i], N)
    V_fromStoA[:,i] = obj.calcvs(rho0, WS[i], CL_max_clean, n_fromStoA[:,i])
    
##########################################################
############ FLIGHT ENVELOPE POINTS - POINT G ############
##########################################################
nG = np.full(n_Mass, nmin) 
VG = np.ones((n_Mass,1))
for i in range(0,n_Mass):
    VG[i] = obj.calcvs(rho0, WS[i], CL_inverted, nG[i])
    
# FLIGHT ENVELOPE - FROM S INVERTED TO G 
n_fromSinvtoG = np.zeros((N,n_Mass))
V_fromSinvtoG = np.zeros((N,n_Mass))
for i in range(0, n_Mass):
    n_fromSinvtoG[:,i] = np.linspace(nSinv[i], nG[i], N)
    V_fromSinvtoG[:,i] = obj.calcvs(rho0, WS[i], CL_inverted, n_fromSinvtoG[:,i])
    
##########################################################
############ FLIGHT ENVELOPE POINTS - POINT C ############
##########################################################
# FLIGHT ENVELOPE - FROM A TO C 
n_fromAtoC = np.zeros((N,n_Mass))
V_fromAtoC = np.zeros((N,n_Mass))
for i in range(0, n_Mass):
    n_fromAtoC[:,i] = np.full(N, nmax)
    V_fromAtoC[:,i] = np.reshape(np.linspace(VA[i], cruise_speed[i], N), (N,))
    
##########################################################
############ FLIGHT ENVELOPE POINTS - POINT C ############
##########################################################
# FLIGHT ENVELOPE - FROM C TO D 
n_fromCtoD = np.zeros((N,n_Mass))
V_fromCtoD = np.zeros((N,n_Mass))
for i in range(0, n_Mass):
    n_fromCtoD[:,i] = np.full(N, nmax)
    V_fromCtoD[:,i] = np.reshape(np.linspace(cruise_speed[i], VD[i], N), (N,))
      
##########################################################
############ FLIGHT ENVELOPE POINTS - POINT D ############
##########################################################
# FLIGHT ENVELOPE - FROM D TO 0 
n_fromDto0 = np.zeros((N,n_Mass))
V_fromDto0 = np.zeros((N,n_Mass))
for i in range(0, n_Mass):
    n_fromDto0[:,i] = np.linspace(nmax, 0.0, N) 
    V_fromDto0[:,i] = np.full(N, VD[i])

##########################################################
############ FLIGHT ENVELOPE POINTS - POINT F ############
##########################################################
# FLIGHT ENVELOPE - FROM G TO F 
n_fromGtoF = np.zeros((N,n_Mass))
V_fromGtoF = np.zeros((N,n_Mass))
for i in range(0, n_Mass):
    n_fromGtoF[:,i] = np.full(N, nmin)
    V_fromGtoF[:,i] = np.reshape(np.linspace(VG[i], cruise_speed[i], N), (N,))    

####################################################################
############ FLIGHT ENVELOPE POINTS - FROM POINT F TO 0 ############
####################################################################
# FLIGHT ENVELOPE - FROM G TO F 
n_fromFto0 = np.zeros((N,n_Mass))
V_fromFto0 = np.zeros((N,n_Mass))
nF         = np.full(n_Mass, nmin)
for i in range(0, n_Mass):
    V_fromFto0[:,i] = np.reshape(np.linspace(cruise_speed[i], VE[i], N), (N,)) 
    x               = np.reshape(np.array([cruise_speed[i], VE[i]]), (2,))
    y               = np.reshape(np.array([nF[i], 0.0]), (2,))
    p_fromFto0      = chebfit(x, y, deg = 1)
    for j in range(0, N):
        n_fromFto0[j,i] = chebval(V_fromFto0[j,i], p_fromFto0) 
        
# -- ###################################################################### --
# -- ### STORE FLIGHT ENVELOPE DATA INSIDE THE SIMPLE NAMESPACE OBJECT #### --
# -- ###################################################################### --
Flight_envelope = {"Flight_envelope": { 
    "n_from0toS": {"Value": n_from0toS, "Unit": "g"},
    "V_from0toS": {"Value": V_from0toS, "Unit": "m/s"}, 
    "n_fromStoA": {"Value": n_fromStoA, "Unit": "g"}, 
    "V_fromStoA": {"Value": V_fromStoA, "Unit": "m/s"},
    "n_fromAtoC": {"Value": n_fromAtoC, "Unit": "g"},
    "V_fromAtoC": {"Value": V_fromAtoC, "Unit": "m/s"}, 
    "n_fromCtoD": {"Value": n_fromCtoD, "Unit": "g"}, 
    "V_fromCtoD": {"Value": V_fromCtoD, "Unit": "m/s"},
    "n_fromDto0": {"Value": n_fromDto0, "Unit": "g"},
    "V_fromDto0": {"Value": V_fromDto0, "Unit": "m/s"},
    "n_from0toSinv": {"Value": n_from0toSinv, "Unit": "g"},
    "V_from0toSinv": {"Value": V_from0toSinv, "Unit": "m/s"}, 
    "n_fromSinvtoG": {"Value": n_fromSinvtoG, "Unit": "g"},
    "V_fromSinvtoG": {"Value": V_fromSinvtoG, "Unit": "m/s"}, 
    "n_fromGtoF"   : {"Value": n_fromGtoF, "Unit": "g"},
    "V_fromGtoF"   : {"Value": V_fromGtoF, "Unit": "m/s"}, 
    "n_fromFto0"   : {"Value": n_fromFto0, "Unit": "g"},
    "V_fromFto0"   : {"Value": V_fromFto0, "Unit": "m/s"},
               "Point S": {"nS": {"Value": nS, "Unit": "g"}, 
                           "VS": {"Value": VS, "Unit": "m/s"}
                                                  },
               "Point A": {"nA": {"Value": nA, "Unit": "g"}, 
                           "VA": {"Value": VA, "Unit": "m/s"}
                                                  }, 
               "Point C": {"nC": {"Value": nC, "Unit": "g"}, 
                           "VC": {"Value": cruise_speed, "Unit": "m/s"}
                                                  },   
               "Point D": {"nD": {"Value": nD, "Unit": "g"}, 
                           "VD": {"Value": VD, "Unit": "m/s"}
                                                  },  
               "Point S inverted": {"nSinv": {"Value": nSinv, "Unit": "g"}, 
                                    "VSinv": {"Value": VSinv, "Unit": "m/s"}
                                                  },  
               "Point G": {"nG": {"Value": nG, "Unit": "g"}, 
                           "VG": {"Value": VG, "Unit": "m/s"}
                                                  }, 
               "Point F": {"nF": {"Value": nF, "Unit": "g"}, 
                           "VF": {"Value": VF, "Unit": "m/s"}
                                                  },  
               "Point E": {"nE": {"Value": nE, "Unit": "g"}, 
                           "VE": {"Value": VE, "Unit": "m/s"}
                                                  },                             
                           } }    
# UPDATING DICTIONARY
aircraft.update(Flight_envelope)
 # UPDATING THE SIMPLE NAMESPACE OBJECT
aircraft_data = SimpleNamespace(**aircraft)      

# PLOTTING RESULTS 
Reg           = aircraft_data.Aircraft_inputs["Regulation"]["Value"]
Aircraft_name = aircraft_data.Aircraft_inputs["Aircraft_name"]["Value"]
n             = 1
for i in range(0,n_Mass):
    figure_name = 'fig' + str(n)
    figure_name = obj.V_n_diagram(n_from0toS[:,i], n_fromStoA[:,i], n_fromAtoC[:,i], n_fromCtoD[:,i], n_fromDto0[:,i],\
                V_from0toS[:,i], V_fromStoA[:,i], V_fromAtoC[:,i], V_fromCtoD[:,i], V_fromDto0[:,i],\
                n_from0toSinv[:,i], n_fromSinvtoG[:,i], n_fromGtoF[:,i], n_fromFto0[:,i],\
                V_from0toSinv[:,i], V_fromSinvtoG[:,i], V_fromGtoF[:,i], V_fromFto0[:,i],\
                Reg, Aircraft_name, n,\
                nS[i], VS[i], nA[i], VA[i], nC[i], cruise_speed[i], nD[i], VD[i],\
                nSinv[i], VSinv[i], nG[i], VG[i], nF[i], VF[i]) 
    n = n + 1
# =============================================================================    
# MOVING FIGURES INSIDE OUTPUT FOLDER
# FlightEnvelope
n = 1
for i in range(0, n_Mass):
    figure_name = r'\FlightEnvelope' + str(n) + '.pdf'
    main_dir    = r'I:\PythonTesiConversion\csvla'
    original    = main_dir + figure_name
    out_dir     = r'\Output\FlightEnvelope'
    target      = main_dir + out_dir + figure_name
    shutil.move(original, target)
    n = n + 1
# ============================================================================= 
###############################   
###### CLOSE ALL FIGURES ######
###############################
plt.close("all")