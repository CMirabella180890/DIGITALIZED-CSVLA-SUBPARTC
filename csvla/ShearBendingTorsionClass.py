# -*- coding: utf-8 -*-
"""
Created on Fri May 13 07:19:47 2022

 SHEAR - BENDING - TORSION Class with simple functions.
   A generic collection of useful function to evaluate Shear and
   Bending moment distributions along the main wing span. 
   In this class, one can find: 
    1. A function to calculate a chord distribution, knowing the wing
       surface Sw, the taper ratio lambda, the span b and a vector
       which contains stations along the span, called y inside this
       file.
    2. Two functions to evaluate axial and normal forces acting on the
       wing per unit length along the wing span. To be used the
       function requires a knowledge of the product c(y)*Cl(y) and 
       c(y)*Cd(y) plus the actual total angle of attack of the wing
       with respect to the flow direction, defined as 
       
                  Alpha_TOT = Alpha(y) + iw(y)
       
       with iw the twist angle of the wing.
    3. Two functions to evaluate shear and bending moment distribution
       along the span of the wing. Shear force distribution is based
       on the normal force distribution along the wing, while bending
       moment is calculated from the shear distribution and a proper
       arm from the wing root/aircraft longitudinal simmetry axis. 
    4. A function that automatically plot the shear and bending moment
       diagram along the wing span. 
 It is important to remember the following flight airworthiness rules:
 FLIGHT LOADS
  CS - VLA 321 
  (a) Flight load factors represent the ratio of the aerodynamic force
      component (acting normal to the assumed longitudinal axis of the
      aeroplane) to the weight of the aeroplane. A positive flight load
      factor is one in which the aerodynamic force acts upward, with
      respect to the aeroplane. 
  (b) Compliance with the flight load requirements of this subpart must
      shown - 
      (1) At each critical altitude within the range in which the
          aeroplane may be expected to operate;
      (2) At each practicable combination of weight and disposable load
          within the operating limitations specified in the Flight Manual.
  
  CS - VLA 331
  (a) The appropriate balancing horizontal tail load must be accounted for
      in a rational or conservative manner when determining the wing loads
      and linear inertia loads corresponding to any of the symmetrical
      flight conditions specified in CS - VLA 331 to 345.
  (b) The incremental horizontal tail loads due to manoeuvring and gusts
      must be reacted by the angular inertia of the aeroplane in a
      rational or conservative manner.
      
@author: claum
"""

import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# plt.rcParams['axes.grid'] = True
import numpy as np

class ShearBendingTorsionClass: 
    #################################################
    #### SHEAR, BENDING AND TORSION CALCULATIONS ####
    #################################################
    
    #################################################
    ########### NORMAL FORCE CALCULATIONS ###########
    #################################################
    def calc_normal_force(AoA_tot, c_Cl, c_Cd, flag):
        """
        A function that calculates normal force coefficients in body axes for 
        any total angle of attack, expressed in radians. The total angle of at-
        tack is given by the sum of the current vehicle angle of attack and the 
        wing geometrical twist (i_w).

        Parameters
        ----------
        AoA_tot : float
            Total angle of attack given by the following:
                AoA_tot = alpha + i_w.
        c_Cl : float
            Distribution of lift coefficients times the chord measured at some 
            station along the wing semi-span.
        c_Cd : float
            Distribution of drag coefficients times the chord measured at some 
            station along the wing semi-span.
        flag : string
            A flag that permits to assess if the AoA_tot is given in degrees or
            in radians. To properly use the np.cos and np.sin functions it must
            be expressed in [rad].

        Returns
        -------
        Normal_force : float
            Normal force acting on the wing.

        """
        
        if ( flag == "rad" ): 
            Normal_force = c_Cl * np.cos(AoA_tot) + c_Cd * np.sin(AoA_tot)
        elif ( flag == "deg" ): 
            AoA_tot_rad  = np.deg2rad(AoA_tot)
            Normal_force = c_Cl * np.cos(AoA_tot_rad) + c_Cd * np.sin(AoA_tot_rad)
        
        return Normal_force
        
    ################################################
    ########### AXIAL FORCE CALCULATIONS ###########
    ################################################
    def calc_axial_force(AoA_tot, c_Cl, c_Cd, flag):
        """
        A function that calculates axial force coefficients in body axes for 
        any total angle of attack, expressed in radians. The total angle of at-
        tack is given by the sum of the current vehicle angle of attack and the 
        wing geometrical twist (i_w).

        Parameters
        ----------
        AoA_tot : float
            Total angle of attack given by the following:
                AoA_tot = alpha + i_w.
        c_Cl : float
            Distribution of lift coefficients times the chord measured at some 
            station along the wing semi-span.
        c_Cd : float
            Distribution of drag coefficients times the chord measured at some 
            station along the wing semi-span.
        flag : string
            A flag that permits to assess if the AoA_tot is given in degrees or
            in radians. To properly use the np.cos and np.sin functions it must
            be expressed in [rad].

        Returns
        -------
        Axial_force : float
            Axial force acting on the wing.

        """
        
        if ( flag == "rad" ): 
            Axial_force = c_Cd * np.cos(AoA_tot) - c_Cl * np.sin(AoA_tot)
        elif ( flag == "deg" ): 
            AoA_tot_rad  = np.deg2rad(AoA_tot)
            Axial_force = c_Cd * np.cos(AoA_tot_rad) - c_Cl * np.sin(AoA_tot_rad)
        
        return Axial_force
            
    ################################################
    ########### SHEAR FORCE CALCULATIONS ###########
    ################################################
    def calc_shear_force(y_halfspan, Normal_force):
        """
        Function that evaluates Shear forces acting on the main wing of 
        an aircraft; the calculations are done as follow: 
        
        Shear[i] = Shear[i-1] + a
        
        where vectors inside this expression must be flipped; this pre-
        liminary operation is done to ensure that the calculation start 
        from wing tip and goes on until wing root is reached. The a is 
        calculated as follow:                
                   
        a = 0.5 * [ Fz[i-1] + Fz[i] ] * [ y[i-1] - y[i] ]
        
        SHEAR CALCULATIONS 
        
        Wing root                            Wing tip
        #
        # ----- / ----- / ----- / ... / ----- /  
        #          y[i]-station          Starting point for calc.
        
        Parameters
        ----------
        y_halfspan : float
            Wing semi-span stations; before these data can be used to
            evaluate shear the array must be flipped.
        Normal_force : float
            Wing normal force values; before these data can be used to
            evaluate shear the array must be flipped.

        Returns
        -------
        Shear : float
            Shear acting on the wing semi-span.

        """
        
        # INITIALIZATION OF THE VECTOR CONTAINING SHEAR FORCE VALUES
        Shear = np.zeros( ( len(y_halfspan), ) )
        
        # WE FIRST FLIP THE Y_HALFSPAN VECTOR and NORMAL FORCE VECTOR 
        # WE SELECT FZ AS A PROPER SYMBOL TO INDICATE NORMAL FORCE
        y_flip  = np.flip(y_halfspan) 
        FZ_flip = np.flip(Normal_force) 
        
        # SHEAR FORCE CALCULATIONS LOOP
        for i in range(1, len(y_halfspan)):
            a  = 0.5 * (FZ_flip[i-1] + FZ_flip[i]) * (y_flip[i-1] - y_flip[i])
            Shear[i] = Shear[i-1] + a
            
        return Shear 
                
    #########################################################
    ########### BENDING MOMENT FORCE CALCULATIONS ###########
    #########################################################
    def calc_bending_moment_force(y_halfspan, Shear):
        """
        Function that evaluates bending moments  acting on the main
        wing of an aircraft; the calculations are done as follow: 
        
        Bend_mom[i] = Bend_mom[i-1] + a
        
        where vectors inside this expression must be flipped; this pre-
        liminary operation is done to ensure that the calculation start 
        from wing tip and goes on until wing root is reached. The a is 
        calculated as follow:                
                   
        a = 0.5 * [ Shear[i-1] + Shear[i] ] * [ y[i-1] - y[i] ]
        
        BENDING MOMENT CALCULATIONS 
        
        Wing root                            Wing tip
        #
        # ----- / ----- / ----- / ... / ----- /  
        #          y[i]-station          Starting point for calc.
        
        Parameters
        ----------
        y_halfspan : float
            Wing semi-span stations; before these data can be used to
            evaluate shear the array must be flipped.
        Shear : float
            Wing shear values; before these data can be used to
            evaluate bending moment the array must be flipped.

        Returns
        -------
        Bend_mom : float
            Bending moment acting on the wing semi-span.

        """
        
        # INITIALIZATION OF THE VECTOR CONTAINING BENDING MOMENTS VALUES
        Bend_mom = np.zeros( ( len(y_halfspan), ) )
        
        # WE FIRST FLIP THE Y_HALFSPAN VECTOR and SHEAR FORCE VECTOR 
        # WE SELECT S AS A PROPER SYMBOL TO INDICATE SHEAR FORCE
        y_flip  = np.flip(y_halfspan) 
        
        # SHEAR FORCE CALCULATIONS LOOP
        for i in range(1, len(y_halfspan)):
            a  = 0.5 * (Shear[i-1] + Shear[i]) * (y_flip[i-1] - y_flip[i])
            Bend_mom[i] = Bend_mom[i-1] + a
            
        return Bend_mom   
                
    ####################################################
    ########### PITCHING MOMENT CALCULATIONS ###########
    #################################################### 
    def calc_m_distr(Cm, Chord, q): 
        """
        Function that calculates the pitching moment distribution along
        the wing semi-span. The calculation is simply performed through
        the following formula: 
            
            m(y) = q * Cm(y) * C(y)^2
         
        This distribution is used to assess torsion along the wing semi-
        span in the next function.

        Parameters
        ----------
        Cm : float
            Pitching moment coefficient distribution along the wing 
            semi-span.
        Chord : float
            Chord distribution along the wing semi-span.
        q : float
            Selected dynamic pressure expressed in [Pa]; it is a single
            constant value.

        Returns
        -------
        m_distr : float
            Dimensional pitching moment distribution along the wing 
            semi-span.

        """
        
        # CALCULATION OF THE PITCHING MOMENT PER UNIT LENGTH
        m_distr = np.ones( (len(Cm), ) )
        
        # CALCULATION LOOP
        for i in range(0, len(Cm)):
            m_distr[i] = q * Cm[i] * ( Chord[i]**2 ) 
            
        return m_distr
                
    ###################################################
    ########### TORSION MOMENT CALCULATIONS ###########
    ################################################### 
    def calc_tors_mom(y_halfspan, m_distr):
        """
        Torsion moment distribution calculator function. This function 
        is relatively similar to the previous functions illustrated in 
        this  file.  The  calculation  technique  is  analogous to the 
        previously adopted one.

        Parameters
        ----------
        y_halfspan : float
            Stations along the wing semi-span; the vector must be 
            flipped before it can be used in the calculation pro-
            cess.
        m_distr : float
            Dimensional pitching moment distribution per unit length
            along the wing semi-span.

        Returns
        -------
        Tors_mom : float
            Torsion moment distribution along the wing semi-span.

        """
        
        # INITIALIZATION
        Tors_mom = np.zeros( (len(y_halfspan), ) )
        
        # FIRST FLIP THE Y_HALFSPAN VECTOR
        y_flip = np.flip(y_halfspan)
        
        # TORSION CALCULATIONS
        for i in range(1, len(y_halfspan)):
            a = 0.5 * (m_distr[i-1] + m_distr[i]) * (y_flip[i-1] - y_flip[i])
            Tors_mom[i] = Tors_mom[i-1] + a
            
        return Tors_mom
                
    ###################################################
    ########### PLOTTING RESULTS - DIAGRAMS ###########
    ################################################### 
    def shear_bending_torsion_diagram(Shear, Bend_mom, Tors_mom, Point,\
                                      y_halfspan, fig_num):
        """
        A function that shows Shear, Bending and Torsion moment distr.
        along the wing halfspan.

        Parameters
        ----------
        Shear : float
            Shear force distribution along the wing semispan.
        Bend_mom : float
            Bending moment distribution along the wing semispan.
        Tors_mom : float
            Torsion moment distribution along the wing semispan.
        Point : float
            Final envelope point at which calculations are performed.
        y_halfspan : float
            Stations along the wing semispan.
        fig_num : string
            Figure number.
        Returns
        -------
        fig : figure
            Subplots with Shear, Bending and Torsion moments diagrams.

        """
        
        # INITIALIZATION 
        y_flip = np.flip(y_halfspan)                
        # TEST - sharex = True - To be checked
        fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
        fig.suptitle(r'Shear, Bending and Torsion distribution' + Point)
        ax1.plot(y_flip, Shear)
        ax1.set_ylabel(r'Shear - $T(y)$ - $[daN]$', fontsize = 'xx-small')
        ax1.grid(visible=True, which="minor")
        ax2.plot(y_flip, Bend_mom)
        ax2.set_ylabel(r'Bending mom. - $M(y)$ - $[daN \cdot m]$', fontsize = 'xx-small')
        ax2.grid(visible=True, which="minor")
        ax3.plot(y_flip, Tors_mom)
        ax3.set_ylabel(r'Torsion mom. - $T(y)$ - $[daN \cdot m]$', fontsize = 'xx-small')
        ax3.set_xlabel(r'Halfspan stations - $y$ - $[m]$')
        ax3.grid(visible=True, which="minor")
        fig.tight_layout() 
        fig_name = 'ShearBendingTorsion' + Point + str(fig_num) + '.pdf'
        plt.savefig(fig_name, bbox_inches='tight')
        
        return fig
                
    def CalcShearBendingTorsion_func(Chord, y_halfspan, Cl, Cd, Cm, alfa,\
                                     twist_angle, q, Point, fig_num, obj):
        """
        A convenient function that applies all the methods defined in-
        side the class ShearBendingTorsionClass.py. Each function called
        has an extensive doc string to clearify its role in current 
        calculations.

        Parameters
        ----------
        Chord : float
            Chord distribution along the wing semispan.
        y_halfspan : float
            Stations along the wing semispan.
        Cl : float
            Lift coefficient distribution along the wing semispan.
        Cd : float
            Drag coefficient distribution along the wing semispan.
        Cm : float
            Pitching mom. distribution along the wing semispan.
        alfa : float
            Angle of attack of the wing.
        twist_angle : float
            Geometrical twist angle of the wing (it could be the wing 
            root twist).
        q : float
            Dynamic pressure at the flight envelope point.
        Point : string
            Final, combined flight envelope.
        numb_case : int
            Case studied (it depends from the number of mass values 
            examined).
        obj : object
            Call to the class ShearBendingTorsionClass.

        Returns
        -------
        c_Cl : float
            Lift coefficient distribution times chord distribution 
            along the wing semispan.
        c_Cd : float
            Drag coefficient distribution times chord distribution 
            along the wing semispan.
        AoA_tot_deg : float
            Total angle of attack expressed in [deg].
        AoA_tot_rad : float
            Total angle of attack expressed in [rad].
        c_Cz : float
            Normal force coefficient distribution times chord distribu-
            tion along the wing semispan.
        c_Ca : float
            Axial force coefficient distribution times chord distribu-
            tion along the wing semispan.
        Normal_force : float
            Normal force distribution along the wing semispan.
        Axial_force : float
            Axial force distribution along the wing semispan.
        Shear : float
            Shear distribution along the wing semispan.
        Bending : float
            Bending moment distribution along the wing semispan.
        m_distr : float
            Distribution of pitching moment along the wing semispan.
        Torsion : float
            Torsion moment distribution along the wing semispan.
        results_plot : figure
            Figure with shear force, bending and torsion moments dia-
            grams collected in a convenient subplots format.

        """
        
        # INITIALIZATION
        c_Cl        = np.multiply(Chord, Cl)
        c_Cd        = np.multiply(Chord, Cd)
        AoA_tot_deg = alfa + twist_angle 
        AoA_tot_rad = np.deg2rad(AoA_tot_deg)
        
        # Calculation of the normal force coefficient. c_Cz = c_Cz(y)
        # This function will be used to evaluate the normal force coefficients
        # distribution along the semi-span; it is possible to find a complete
        # documentation inside this file.
        flag         = "rad"
        c_Cz         = obj.calc_normal_force(AoA_tot_rad, c_Cl, c_Cd, flag)
        Normal_force = np.multiply(c_Cz, q)
        
        # Calculation of the axial force coefficient. c_Ca = c_Ca(y) 
        # This function will be used to evaluate the axial force coefficients
        # distribution along the span; it is possible to find a complete
        # documentation inside this file.
        flag        = "rad"
        c_Ca        = obj.calc_axial_force(AoA_tot_rad, c_Cl, c_Cd, flag)
        Axial_force = np.multiply(c_Ca, q)
        
        # Calculation of the shear force distribution along the wing 
        # semispan. A complete description of this function can be 
        # found inside this file.
        Shear = obj.calc_shear_force(y_halfspan, Normal_force) * 1E-1
        
        # Calculation of the bending mom. distribution along the wing 
        # semispan. A complete description of this function can be 
        # found inside this file.
        Bending = obj.calc_bending_moment_force(y_halfspan, Shear)
        
        # Calculation of the pitching moment distribution along the wing 
        # semispan. The complete description of this function can be 
        # found inside this file.
        m_distr = obj.calc_m_distr(Cm, Chord, q)
        
        # Calculation of the Torsion moment distribution along the wing 
        # semispan. The complete description of this function can be 
        # found inside this file.
        Torsion = obj.calc_tors_mom(y_halfspan, m_distr) * 1E-1 
        
        results_plot = obj.shear_bending_torsion_diagram(Shear, Bending,\
                                                         Torsion, Point,\
                                                         y_halfspan, fig_num)
        
        return c_Cl, c_Cd, AoA_tot_deg, AoA_tot_rad, c_Cz, c_Ca,\
               Normal_force, Axial_force, Shear, Bending, m_distr,\
               Torsion, results_plot
                
                
            
                       
                