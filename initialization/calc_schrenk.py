# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 14:52:02 2022

@author: claum
"""

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# ========================================

def calc_schrenk_distribution(b_wing, S_wing, chord, y_halfspan):
    """
    A function that produce Schrenk's loads distribution (unitary CL 
    distribution along wing half span).

    Parameters
    ----------
    b_wing : float
        Wing span [m].
    S_wing : float
        Wing surface [squared meters].
    chord : float
        Chord distribution along the win semi span [m].
    y_halfspan : float
        Stations distribution along the wing semi span [m].

    Returns
    -------
    eta : float
        Non dimensional stations distribution along the wing semi span [m].
    elliptical_load : flot
        Elliptical load distribution.
    schrenk_cCL : float
        Schrenk load times chord distribution [m].
    unit_CL : float
        Unitary lift coefficient distribution.

    """
    
    # NON DIMENSIONAL STATIONS ALONG THE WING SEMI SPAN 
    eta = (y_halfspan) / (0.5 * b_wing)
    
    # ELLIPTICAL LOAD 
    elliptical_load = ((4 * S_wing) / (np.pi * b_wing)) * np.sqrt(1 - eta**2)
    
    # SCHRENK CHORD * CL 
    N = len(chord)
    schrenk_cCL = np.ones(N)
    for i in range(N-1):
        schrenk_cCL[i] = 0.5 * (chord[i] + elliptical_load[i])
    
    schrenk_cCL[-1] = 0
    
    # CHECKING RESULTS 
    global_CL = (2 / S_wing) * np.trapz(schrenk_cCL, y_halfspan)
    if (global_CL < 1.0) or (global_CL > 1.0):
        global_CL   = 1.0 / global_CL 
        unit_CL     = (global_CL * schrenk_cCL) / (chord)
        schrenk_cCL = unit_CL * chord 
        
    print(" Global unit lift coefficient: \n", " ", (2 / S_wing) * np.trapz(schrenk_cCL, y_halfspan))
    
    return eta, elliptical_load, schrenk_cCL, unit_CL

def plot_schrenk(y_halfspan, chord, elliptical_load, schrenk_cCL, unit_CL): 
    
    # =======================================
    # ======= CL INTERPOLATION FIGURE =======
    # =======================================
    fig4 = plt.figure()
    plt.plot(y_halfspan, chord, color="black", linewidth=1.0, linestyle = "dashed", label = r"Chord $c(y)$")
    plt.plot(y_halfspan, elliptical_load, color="blue", linewidth=1.0,\
             linestyle = "dashed", label = r"Elliptical load - $C_{l}(y)$")
    plt.plot(y_halfspan, unit_CL, color="red", linewidth=1.0, label = r"Unitary load - $C_{l}(y)$")
    plt.xlabel(r'Station along wing semispan - $y$ [m]')                # x-label to the axes.
    plt.ylabel(r'Spanwise lift coefficient distribution - $C_{l}(y)$')  # y-label to the axes.
    plt.title(r'Lift coefficient curve distributions - Schrenk method') # Title to the axes.
    plt.minorticks_on()
    plt.grid(True, linestyle='-.', which="both")
    plt.legend(loc="best")
    plt.savefig('figure4.pdf', bbox_inches='tight')
    plt.show()
    # =======================================
    
    # =======================================
    # ======= CL INTERPOLATION FIGURE =======
    # =======================================
    fig5 = plt.figure()
    plt.plot(y_halfspan, schrenk_cCL)
    plt.xlabel(r'Station along wing semispan - $y$ [m]')                  # x-label to the axes.
    plt.ylabel(r'Lift coeff. distr. times chord - $c(y) \cdot C_{l}(y)$') # y-label to the axes.
    plt.title(r'Lift coeff. distr. times chord - Schrenk method')         # Title to the axes.
    plt.minorticks_on()
    plt.grid(True, linestyle='-.', which="both")
    # plt.legend(loc="best")
    plt.savefig('figure5.pdf', bbox_inches='tight')
    plt.show()
    # =======================================
    
    return fig4, fig5