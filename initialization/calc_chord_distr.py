# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 07:38:54 2022

@author: claum
"""
import math
import numpy as np 

def initialization_calc_chord(aircraft_data):
    
    # Main wing parameters 
    b_wing         = aircraft_data.Wing["b"]["Value"]
    croot          = aircraft_data.Wing["c_root"]["Value"]
    ctip           = aircraft_data.Wing["c_tip"]["Value"]
    S_wing         = 0.5 * ( croot + ctip ) * b_wing
    AR_wing        = ( b_wing**2 ) / (S_wing) 
    MGC_wing       = S_wing / b_wing 
    chord_kink_one = aircraft_data.Wing["c_kink_one"]["Value"]
    chord_kink_two = aircraft_data.Wing["c_kink_two"]["Value"]
    wing_type_flag = aircraft_data.Wing["wing_type"]["Value"]
    taper_ratio1   = chord_kink_one / croot  
    taper_ratio2   = chord_kink_two / chord_kink_one  
    taper_ratio3   = ctip / chord_kink_two     
                                                                             
    # Horizontal tail parameters 
    ht_span        = aircraft_data.Horizontal_Tailplane["ht_span"]["Value"]
    ht_croot       = aircraft_data.Horizontal_Tailplane["ht_croot"]["Value"]
    ht_ctip        = aircraft_data.Horizontal_Tailplane["ht_ctip"]["Value"]
    ht_taper_ratio = ht_croot / ht_ctip
    ht_S           = 0.5 * ( ht_croot + ht_ctip ) * ht_span
    ht_AR          = ( ht_span**2 ) / ( ht_S )
    ht_MGC         = ht_S / ht_span  
                                                                             
    # Vertical tail parameters 
    vt_span        = aircraft_data.Vertical_Tailplane["vt_span"]["Value"]
    vt_croot       = aircraft_data.Vertical_Tailplane["vt_croot"]["Value"]
    vt_ctip        = aircraft_data.Vertical_Tailplane["vt_ctip"]["Value"]
    vt_taper_ratio = vt_croot / vt_ctip
    vt_S           = 0.5 * ( vt_croot + vt_ctip ) * vt_span
    vt_AR          = ( vt_span**2 ) / ( vt_S )
    vt_MGC         = vt_S / vt_span 
    
    return b_wing, croot, ctip, S_wing, AR_wing, MGC_wing, chord_kink_one,\
        chord_kink_two, wing_type_flag, taper_ratio1, taper_ratio2, taper_ratio3,\
            ht_span, ht_croot, ht_ctip, ht_taper_ratio, ht_S, ht_AR, ht_MGC,\
                vt_span, vt_croot, vt_ctip, vt_taper_ratio, vt_S, vt_AR, vt_MGC
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++ FUNCTION THAT CALCULATES CHORD DISTRIBUTIONS ++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
def chord_distribution(b_wing, S_wing, N, croot, ctip, chord_kink_one,\
                           chord_kink_two, wing_type_flag, taper_ratio1,\
                               taper_ratio2, taper_ratio3): 
    """
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    +++++++ GENERAL FUNCTION TO EVALUATE CHORD DISTRIBUTION OF A WING +++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Parameters
    ----------
    b_wing : float
        Wing span.
    S_wing : float
        Wing surface.
    N : int
        Number of stations along the wing semi-span.
    croot : float
        Root chord.
    ctip : float
        Tip chord.
    chord_kink_one : float
        Chord measured at the first kink.
    chord_kink_two : float
        Chord measured at the second kink.
    wing_type_flag : float
        Flag that identifies wing planform.
    taper_ratio1 : float
        First taper ratio: lambda1 = chord_kink1 / croot.
    taper_ratio2 : float
        Second taper ratio: lambda2 = chord_kink2 / chord_kink1.
    taper_ratio3 : float
        Third taper ratio: lambda3 = ctip / chord_kink2.

    Returns
    -------
    chord : float
        Chord distribution along the wing semi-span.
    y_halfspan : float
        Stations along the semi-span.

    """
    # Half-span vector 
    b_half     = 0.5 * b_wing 
    y_halfspan = np.linspace(0.0, b_half, N)
    # CHORD CALCULATOR FUNCTION 
    chord_calculator = lambda S_wing, b_wing, taper_ratio, y_halfspan: ( (2 * S_wing) / ( (1 + taper_ratio) * b_wing)) * (1 - (((1 - taper_ratio))/b_wing) * np.abs(2 * y_halfspan))
    
    # Switch to assess wing type and chord calculator
    if (wing_type_flag == "Rectangular"):
        
        taper_ratio = ctip / croot 
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # + CALCULATION OF A CHORD DISTRIBUTION WITH A CONVENIENT, SIMPLE FUNCTION +
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        chord = chord_calculator(S_wing, b_wing, taper_ratio, np.flip(y_halfspan))    
    elif (wing_type_flag == "With_kins"):
        if (chord_kink_one != chord_kink_two): 
            
            # SETTING UP Y ALONG THE SPAN 
            y_halfspan1 = y_halfspan[1:math.ceil(N/3)]
            y_halfspan2 = y_halfspan[math.ceil(N/3)+1:math.ceil(2*N/3)]
            y_halfspan3 = y_halfspan[math.ceil(2*N/3)+1:-1]
            
            b_wing1     = y_halfspan1[-1]
            b_wing2     = y_halfspan2[-1]
            b_wing3     = y_halfspan3[-1]
            
            S_wing1     = 0.5 * (croot + chord_kink_one) * b_wing1
            S_wing2     = 0.5 * (chord_kink_two + chord_kink_one) * b_wing2
            S_wing3     = 0.5 * (ctip + chord_kink_two) * b_wing3
            
            # CHORD DISTRIBUTION CALCULATIONS
            chord1 = chord_calculator(S_wing1, b_wing1, taper_ratio1, np.flip(y_halfspan1))
            chord2 = chord_calculator(S_wing2, b_wing2, taper_ratio2, np.flip(y_halfspan2))
            chord3 = chord_calculator(S_wing3, b_wing3, taper_ratio3, np.flip(y_halfspan3))
            
            # ASSEMBLING THE VECTOR 
            chord = np.concatenate((chord1, chord2, chord3))          
            
        elif (chord_kink_one == chord_kink_two):
            
            # SETTING UP Y ALONG THE SPAN 
            y_halfspan1 = y_halfspan[1:math.ceil(N/2)]
            y_halfspan2 = y_halfspan[math.ceil(N/2)+1:-1]
            
            b_wing1     = y_halfspan1[-1]
            b_wing2     = y_halfspan2[-1]
            
            S_wing1     = 0.5 * (croot + chord_kink_one) * b_wing1
            S_wing2     = 0.5 * (chord_kink_two + chord_kink_one) * b_wing2
            
            # CHORD DISTRIBUTION CALCULATIONS
            chord1 = chord_calculator(S_wing1, b_wing1, taper_ratio1, np.flip(y_halfspan1))
            chord2 = chord_calculator(S_wing2, b_wing2, taper_ratio2, np.flip(y_halfspan2))
            
            # ASSEMBLING THE VECTOR 
            chord = np.concatenate((chord1, chord2))        
    
    return chord, y_halfspan

def mean_aerodynamic_chord(S_wing, b_wing, taper_ratio, chord, y_halfspan):
    
    """
        +++++++++++++++++++++++++++++++++++++++
        + FORMULA:                            +
        +           2       /               \ +
        + M.A.C. = --- * Int| ( c(y) )^2 dy | +
        +           S       \               / +
        +++++++++++++++++++++++++++++++++++++++

    Parameters
    ----------
    S_wing : float
        Wing surface.
    chord : float
        Chord distributions along the wing half span.
    y_halfspan : float
        Stations along the wing half span.

    Returns
    -------
    y_MAC : float
        MAC position along the wing half span.
    MAC : float
        Mean aerodynamic chord.

    """
    chord_squared = np.power(chord, 2)
    MAC           = (2 / S_wing) * np.trapz(chord_squared, y_halfspan)
    y_MAC         = ( b_wing / 6 ) * ( (1 + 2 * taper_ratio) / (1 + taper_ratio) )
    return y_MAC, MAC
    