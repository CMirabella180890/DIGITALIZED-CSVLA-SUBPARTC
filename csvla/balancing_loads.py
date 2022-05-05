# -*- coding: utf-8 -*-
"""
Created on Tue May  3 06:57:13 2022

@author: claum
"""
from numpy.polynomial.chebyshev import chebfit, chebval
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

class balancing_loads: 
    """

    + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    + BALANCING LOADS CLASS                                                     +
    + This class collects methods and functions involved during balancing loads +
    + calculations. Balancing loads are loads acting upon main wing and horiz.  + 
    + tailplane during straight, non-accelerated flight. These are the base     +
    + loads to size aero-structures; manoeuvring and gust loads must be added   +
    + to the balancing loads, as prescribed by airworthiness regulations. In    +
    + this case, onli CS - VLA will be considered.                              +
    + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    
    CS - VLA 331 Symmetrical flight conditions 
    
    (a) The appropriate balancing horizontal tail load must be accounted for in
        a rational or  conservative  manner when determining  the win loads and 
        linear  inertia loads  corresponding  to any of the  symmetrical flight 
        conditions specified in CS - VLA 331 to 345.
        
    (b) The incremental horizontal tail loads due to manoeuvring and gusts must
        be reacted  by the angular  inertia of the aeroplane  in a rational and 
        conservative manner.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    %   INPUT                                                                   %
    %   rho              --> Air density                                        %
    %   V                --> Airspeed from the flight envelope                  %
    %   WS               --> Wing loading                                       %
    %   n                --> Load factor from the flight envelope               %
    %   CLmax_aero_model --> Max lift coefficient from the aerodynamic model    %
    %   a2               --> Interpolation coefficient, quadratic term          %
    %   b2               --> Interpolation coefficient, linear term             %
    %   c2               --> Interpolation coefficient, constant term           %
    %   CL_star          --> Lift coefficient at the end of the linear part of  %
    %                        the lift curve                                     %
    %   CL0              --> Lift at zero angle of attack                       %
    %   CLalfa           --> Lift curve slope in 1/deg                          %
    %   p_CD             --> Interpolation coefficients for the drag polar      %
    %   XAC              --> Aerodynamic centre abscissa in MAC percentage      %
    %   XCG              --> Centre of gravity absissa in MAC percentage        %
    %   bCG              --> Dimensional, vertical arm of the lift force        %
    %                        horizontal component; this must be divided by the  %
    %                        Mean Aerodynamic chord before to apply the         %
    %                        equation for the pitching mom. equilibrium         %
    %   MAC              --> Mean aerodynamic chord                             %
    %   h                --> Moment arm measured in meters of the thrust action %
    %                        line with respect to the center of gravity; this   %
    %                        value must be divided by the Mean Aerodynamic      %
    %                        chord                                              %
    %   CM0              --> Pitching moment coefficient at zero angle of       %
    %                        attack                                             %
    %   CM_landing_gear  --> Pitching moment contribute to the landing gear     %
    %   CM_slope         --> Slope of the pitching moment curve (a straight     %
    %                        line)                                              %
    %   S                --> Wing surface                                       %
    %   S_ht             --> Horizontal tail surface                            % 
    %   l_ht             --> Horizontal tail moment arm                         %
    %                                                                           %
    %   OUTPUT                                                                  % 
    %   CL_wb  --> Wing body lift coefficient                                   %
    %   alfa   --> Angle of attack                                              %
    %   CD     --> Wing body drag coefficient                                   %
    %   q      --> Dynamic pressure                                             %
    %   CM_CL  --> Lift pitching moment coefficient contribute                  %
    %   CM_CD  --> Drag pitching moment coefficient contribute                  %
    %   CM_CT  --> Thrust coefficient pitching moment coefficient contribute    %
    %   CM_CG  --> Total pitching moment coefficient                            %
    %   CL_ht  --> Lift coefficient of the horizontal tail                      %
    %   CL_new --> Full-vehicle lift coefficient                                %
    %   L_wb   --> Wing-body lift in daN                                        %
    %   L_new  --> Full-vehicle lift in daN                                     %
    %   L_ht   --> Horizontal tail in daN                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    """
    def alpha_fullmodel(p_lift1, p_lift2, CL, CL_max, CL_star, alfa0l, flag):
        """
        alfa = alpha_fullmodel(p_lift1, p_lift2, CL, CL_max, CL_star, alfa0l, flag)
        
        Function that calculates the angle of attack corresponding to the flight 
        lift coefficient, as prescribed by the interpolated aerodynamic model.

        Parameters
        ----------
        p_lift1 : array of float
            Interpolation coefficients relative to the drag curve. The variable
            is an array  containing Chebyshev coefficients ordered  from low to 
            high. The array is structured as follow: 
                                   [c_0, c_1, ... , c_n]
            with n the degree(s) of the fitting polynomials.
        p_lift2 : array of float
            Interpolation coefficients relative to the drag curve. The variable
            is an array  containing Chebyshev coefficients ordered  from low to 
            high. The array is structured as follow: 
                                   [c_0, c_1, ... , c_n]
            with n the degree(s) of the fitting polynomials.
        CL : float
            Lift coefficient.
        CL_max : float
            Maximum lift coefficient value calculated from the aerodynamic model.
        CL_star : float
            Lift coefficient at the end of the linear part of the lift curve.
        alfa0l : float
            Zero lift angle of attack in [deg].
        flag : string
            A flag to assess straight or inverted flight.

        Returns
        -------
        alfa : float
            Angle of attack calculated from the linear or quadratic lift coeff.
            curve,  as  defined  by  the aerodynamic  model  interpolation 
            coefficients.

        """
        
        if (flag == "straight"):
            
            # INITIALIZATION
            CL0          = p_lift1[0]
            CLalfa_deg   = p_lift1[1] 
            c            = p_lift2[0]
            b            = p_lift2[1]
            a            = p_lift2[2]
            alfa_lin     = lambda CL, CL0, CLalfa_deg: ( CL - CL0 ) / CLalfa_deg
            alfa_non_lin = lambda CL, a, b, c: (-b + np.sqrt(b**2 - 4*a*(c - CL)))/(2*a)
            
            # CHECKING THE LIFT COEFFICIENT VALUE
            if ( CL > CL_max ):
                print("\n ERROR: There is no solution for alfa! \n")
                print(" CL will be considered equal to CL_max! \n ")
                CL = CL_max
                
            # BACK ENVELOPE ALONG THE WING LIFT CURVE 
            # If the lift coefficient is less than CL_star, the resulting
            # lift coefficient is calculated from the linear part of the
            # lift curve. If the lift coefficient is higher than CL_star,
            # the resulting lift coefficient is calculated from the
            # inerpolation data (parabolic curve-fitting).  
            if ( CL <= CL_star ):
                
                alfa = alfa_lin(CL, CL0, CLalfa_deg)
                
            elif ( CL > CL_star ):
                
                alfa = alfa_non_lin(CL, a, b, c)
                
                if ( np.isreal(alfa) == False ):
                    print("\n WARNING: Negative number under square root. \n" )
                    print(" CL is greater than CL_max wing-body \n")
                    exit()
                elif ( np.isreal(alfa) == True ):
                    print(" Alfa value is correct. \n")
            
        elif (flag == "inverted"):
            
            # INITIALIZATION
            CL0          = p_lift1[0]
            CLalfa_deg   = p_lift1[1] 
            CL_max_inv   = -CL_max
            alfa_lin_inv = lambda CL, CL0, CLalfa_deg: ( CL - CL0 ) / CLalfa_deg
            # alfa_lin_inv = lambda CL, CL0, CLalfa_deg, alfa0l: (( CL - CL0 ) / CLalfa_deg) - alfa0l
            
            # BACK ENVELOPE ALONG THE WING LIFT CURVE 
            # If the lift coefficient is less than CL_star, the resulting
            # lift coefficient is calculated from the linear part of the
            # lift curve. If the lift coefficient is higher than CL_star,
            # the resulting lift coefficient is calculated from the
            # inerpolation data (parabolic curve-fitting).  
            if ( abs(CL) <= CL_max_inv ):
                
                alfa = alfa_lin_inv(CL, CL0, CLalfa_deg) - alfa0l
                # alfa = alfa_lin_inv(CL, CL0, CLalfa_deg, alfa0l)
                
            elif ( abs(CL) > CL_max_inv ):
                
                CL   = -CL_max_inv
                alfa = alfa_lin_inv(CL, CL0, CLalfa_deg) - alfa0l
                # alfa = alfa_lin_inv(CL, CL0, CLalfa_deg, alfa0l) 
        
        return alfa

    #########################################################################
    ############### PITCHING MOMENT CONTRIBUTION DUE TO LIFT ################
    #########################################################################    
    def CL_max_function(rho, V, WS, n, CLmax_aero_model):
        """
        % CL_wb = CL_max_function(rho, V, WS, n)
        %   A function that calculates the maximum lift coefficient from airspeed
        %   and load factor values. The user must check the calculated values to
        %   avoid that the maximum equilibrium lift coefficient results higher than
        %   the maximum lift coefficient defined in the aerodynamic model.
        %
        %   INPUT 
        %   rho              --> Air density (typically sea level)
        %   V                --> Airspeed values from the combined flight envelope
        %   WS               --> Wing loading expressed in N/m^2 
        %   n                --> Load factor values from the combined flight envelope
        %   CLmax_aero_model --> Maximum lift coefficient relative to the aero
        %                        model.
        %
        %   OUTPUT
        %   CL_wb --> Maximum equilibrium lift coefficient values; the resulting 
        %             lift  coefficient values  are related  to the wing - body 
        %             configuration         

        Parameters
        ----------
        rho : TYPE
            DESCRIPTION.
        V : TYPE
            DESCRIPTION.
        WS : TYPE
            DESCRIPTION.
        n : TYPE
            DESCRIPTION.
        CLmax_aero_model : TYPE
            DESCRIPTION.

        Returns
        -------
        CL_wb : TYPE
            DESCRIPTION.

        """
        # WING BODY LIFT COEFFICIENT
        CL_wb = (2 / rho) * ( 1 / V**2 ) * WS * n
        
        # CHECKING LIFT COEFFICIENT VALUE AGAINST THE AERO MODEL 
        if ( abs(CL_wb) <= abs(CLmax_aero_model) ): 
            print(" \n Wing body lift coefficient under the CLmax \n")
            CL_wb = CL_wb
        elif ( abs(CL_wb) >= abs(CLmax_aero_model) ):  
            print (" \n CLwb = CLmax_aero_model \n")
            CL_wb = CLmax_aero_model
        
        return CL_wb
        
    #########################################################################
    ############### PITCHING MOMENT CONTRIBUTION DUE TO LIFT ################
    #########################################################################    
    def CLWB_pitch_contrib(CL, alfa, xAC, xCG, bCG, MAC, flag):
        """
        % Lift contribution to pitching mom. = CLWBpitch_contrib(CL, alfa, xAC,\
        %                                                        xCG, bCG, MAC)
        % 
        %  This function is able to calculate the lift contribution
        %  related to the pitching moment equilibrium during a
        %  particular flight condition, defined in terms of airspeed V
        %  and load factor n. The equilibrium is taken about the
        %  aircraft center of gravity. This function take into account 
        %  all the terms that can be derived from simplified
        %  representation of the equilibrium condition; projection of 
        %  lift in a two dimensional fashion is performed, and all the 
        %  contributions are summed. 
        %  To be more specific, two contribution to the pitching moment
        %  are relate to the lift produced by the main lifting surface: 
        %  - a vertical component, which produce pitching moment if its
        %    arm with respect to the aircraft center of gravity is
        %    different from zero; 
        %  - an horizontal component, which produce pitching moment if
        %    its arm with respect to the aircraft center of gravity is
        %    different from zero. It must be remembered that, for a
        %    proper designed aircraft, the lift force horizontal
        %    component is relatively small.
        %  
        %  INPUTS 
        %  CL    --> Input lift coefficient, related to a particular
        %            flight condition which is completely defined in
        %            terms of airspeed V and load factor n (aircraft
        %            weight is a known quantity). It must be 
        %            remembered that the CL used here is relative to
        %            a wing - body configuration (wing - fuselage CL
        %            without horizontal empennage)
        %  alfa  --> Input angle of attack, used in this function to
        %            correctly project along wind axes the lift force 
        %            components in its vertical and horizontal
        %            component
        %  xAC   --> Position along the X - Axis of the aircraft 
        %            aerodynamic center; usually xac = 0.25 for a
        %            typical, subsonic airfoil when Xac is divided by
        %            the Mean Aerodynamic Chord.
        %  xCG   --> Position along the X - Axis of the aircraft center
        %            of gravity; for typical overall geom. arrangement
        %            and loading condition, the xcg could be anywhere
        %            between 0.2 and 0.5 when Xcg is divided by the
        %            Mean Aerodynamic chord
        %  bCG   --> Dimensional, vertical arm of the lift force
        %            horizontal component; this must be divided by the
        %            Mean Aerodynamic chord before to apply the
        %            equation for the pitching mom. equilibrium
        %  MAC   --> Mean Aerodynamic chord
        %  flag  --> Flag to assess how to deal with angles (expressed in [deg]
        %            or in [rad]). Values:
        %             ----> deg
        %             ----> rad
        %
        %  OUTPUTS
        %  LIFT PITCHING MOMENT CONTRIBUTION
        %  It's the main lifting surface contribution to pitching
        %  moment and its lift.
        %
        %  The function calculates the lift coefficient contribution to 
        %  the pitching moment contribution is:
        %
        %                                                          bcg
        %  CLWB_contr = CL*cos(alpha)*(xac - xcg) - CL*sin(alpha)* ---
        %                                                          MAC  
        """       
        if ( flag == "deg"):
            
            alfa_rad = np.deg2rad(alfa)
        
            temp1 = bCG / MAC
            temp2 = xAC - xCG
            
            CM_due_to_CL = CL * temp2 * np.cos(alfa_rad) - CL * temp1 * np.sin(alfa_rad)
        elif ( flag == "rad" ): 
            
            alfa_rad = alfa
        
            temp1 = bCG / MAC
            temp2 = xAC - xCG
            
            CM_due_to_CL = CL * temp2 * np.cos(alfa_rad) - CL * temp1 * np.sin(alfa_rad)
        
        return CM_due_to_CL
    
    #########################################################################
    ############### PITCHING MOMENT CONTRIBUTION DUE TO DRAG ################
    ######################################################################### 
    def CDWB_pitch_contrib(CD, alfa, xAC, xCG, bCG, MAC, flag):
        """
        % Drag contribution to pitching mom. = CDWBpitch_contrib(CD, alfa, xAC,\
        %                                                        xCG, bCG, MAC)
        % 
        %  This function is able to calculate the drag contribution
        %  related to the pitching moment equilibrium during a
        %  particular flight condition, defined in terms of airspeed V
        %  and load factor n. The equilibrium is taken about the
        %  aircraft center of gravity. This function take into account 
        %  all the terms that can be derived from simplified
        %  representation of the equilibrium condition; projection of 
        %  drag in a two dimensional fashion is performed, and all the 
        %  contributions are summed. 
        %  To be more specific, two contribution to the pitching moment
        %  are relate to the drag produced by the main lifting surface: 
        %  - a vertical component, which produce pitching moment if its
        %    arm with respect to the aircraft center of gravity is
        %    different from zero. It must be remembered that, for a
        %    typical design, the drag force vertical component is 
        %    relatively small;
        %  - an horizontal component, which produce pitching moment if
        %    its arm with respect to the aircraft center of gravity is
        %    different from zero. 
        %  
        %  INPUTS 
        %  CD    --> Input drag coefficient, related to a particular
        %            flight condition which is completely defined in
        %            terms of airspeed V and load factor n (aircraft
        %            weight is a known quantity). It must be 
        %            remembered that the CL used here is relative to
        %            a wing - body configuration (wing - fuselage CL
        %            without horizontal empennage)
        %  alfa  --> Input angle of attack, used in this function to
        %            correctly project along wind axes the lift force 
        %            components in its vertical and horizontal
        %            component
        %  xAC   --> Position along the X - Axis of the aircraft 
        %            aerodynamic center; usually xac = 0.25 for a
        %            typical, subsonic airfoil when Xac is divided by
        %            the Mean Aerodynamic Chord.
        %  xCG   --> Position along the X - Axis of the aircraft center
        %            of gravity; for typical overall geom. arrangement
        %            and loading condition, the xcg could be anywhere
        %            between 0.2 and 0.5 when Xcg is divided by the
        %            Mean Aerodynamic chord
        %  bCG   --> Dimensional, vertical arm of the lift force
        %            horizontal component; this must be divided by the
        %            Mean Aerodynamic chord before to apply the
        %            equation for the pitching mom. equilibrium
        %  MAC   --> Mean Aerodynamic chord
        %  flag  --> Flag to assess how to deal with angles (expressed in [deg]
        %            or in [rad]). Values:
        %             ----> deg
        %             ----> rad
        %
        %  OUTPUTS
        %  DRAG PITCHING MOMENT CONTRIBUTION
        %  It's the main lifting surface contribution to pitching
        %  moment and its drag.
        %
        %  The function calculates the drag coefficient contribution to 
        %  the pitching moment contribution is:
        %
        %                              bcg
        %  CDWB_contr = CD*cos(alpha)* --- - CD*sin(alpha)*(xac - xcg)
        %                              MAC  
        
        """   
        if ( flag == "deg" ):
            
            alfa_rad = np.deg2rad(alfa)
            
            temp1 = bCG / MAC
            temp2 = xAC - xCG
            
            CM_due_to_CD = CD * temp1 * np.cos(alfa_rad) - CD * temp2 * np.sin(alfa_rad)
            
        elif ( flag == "rad" ):
            
            alfa_rad = alfa
            
            temp1 = bCG / MAC
            temp2 = xAC - xCG
            
            CM_due_to_CD = CD * temp1 * np.cos(alfa_rad) - CD * temp2 * np.sin(alfa_rad)
        
        return CM_due_to_CD

    ###########################################################################
    ############### PITCHING MOMENT CONTRIBUTION DUE TO THRUST ################
    ########################################################################### 
    def CTWB_pitch_contrib(CD, h, MAC):
        """
        % Thrust contribution to pitching mom. = CTpitch_contrib(obj, CD, h, MAC)
        % 
        %  This function is able to calculate the thrust contribution
        %  related to the pitching moment equilibrium during a
        %  particular flight condition, defined in terms of airspeed V
        %  and load factor n. The equilibrium is taken about the
        %  aircraft center of gravity. To model the thrust, we start 
        %  the obvious consideration that thrust must be at least equal
        %  to drag for sustained, levelled, steady flight; therefore,
        %  in this function the thrust coefficient is considered equal
        %  to the drag coefficient. The moment arm of thrust action 
        %  line with respect to the aircraft center of gravity is
        %  indicated with h.
        %  
        %  INPUTS 
        %  CD    --> Input drag coefficient, related to a particular
        %            flight condition which is completely defined in
        %            terms of airspeed V and load factor n (aircraft
        %            weight is a known quantity). It must be 
        %            remembered that the CD used here is relative to
        %            a wing - body configuration (wing - fuselage CD
        %            without horizontal empennage)
        %  h     --> Moment arm measured in meters of the thrust action
        %            line with respect to the center of gravity; this
        %            value must be divided by the Mean Aerodynamic
        %            chord 
        %  MAC   --> Mean Aerodynamic chord
        %
        %  OUTPUTS
        %  THRUST PITCHING MOMENT CONTRIBUTION
        %  It's the pitching moment contribution associated with thrust
        %  produced by the power train.
        %
        %  The function calculates the drag coefficient contribution to 
        %  the pitching moment contribution is:
        %
        %                  h
        %  CTcontr = CD * ---
        %                 MAC    
        """
        
        temp1 = h / MAC
        
        CM_due_to_CT = CD * temp1
        
        return CM_due_to_CT

    ###################################################################
    ############### GLOBAL PITCHING MOMENT COEFFICIENT ################
    ###################################################################  
    def CM_aboutcg(CM0, CM_gear, CMCL, CL):
        """
        % Pitching mom. about c.g. = CM_aboutcg(obj, CM0, CM_landgear, CMCL, CL)
        % 
        %  This function is able to calculate the thrust contribution
        %  related to the pitching moment equilibrium during a
        %  particular flight condition, defined in terms of airspeed V
        %  and load factor n. The equilibrium is taken about the
        %  aircraft center of gravity. To model the thrust, we start 
        %  the obvious consideration that thrust must be at least equal
        %  to drag for sustained, levelled, steady flight; therefore,
        %  in this function the thrust coefficient is considered equal
        %  to the drag coefficient. The moment arm of thrust action 
        %  line with respect to the aircraft center of gravity is
        %  indicated with h.
        %  
        %  INPUTS 
        %  CM0         --> Zero - lift pitching moment of the wing -
        %                  body configuration 
        %  CM_gear     --> Pitching moment associated with the
        %                  extendend landing gear flight condition; it
        %                  is a constant value
        %  CMCL        --> Pitching moment curve slope, defined in
        %                  terms of lift coefficient variation
        %                                              d Cm
        %                  Pitching mom. curve slope = ----
        %                                              d CL
        %  CL          --> Input lift coefficient, related to a
        %                  particular flight condition which is 
        %                  completely defined in terms of airspeed V 
        %                  and load factor n (aircraft weight is a 
        %                  known quantity). It must be remembered that
        %                  the CL used here is relative to a 
        %                  wing - body configuration (wing - fuselage 
        %                  CL without horizontal empennage)
        %
        %  OUTPUTS
        %  PURE PITCHING MOMENT CONTRIBUTION
        %  Pure pitching moment associated with a typical flight 
        %  vehicle.
        %
        %  The function calculates the CM wing - body coefficient:
        %
        %                 d Cm
        %  CMcontr = CL * ---- + CM0 + (CM)_Landing_Gear
        %                 d CL 
        """ 
        CM_CG = CL * CMCL + CM_gear + CM0
        
        return CM_CG

    ###################################################################
    ############### GLOBAL PITCHING MOMENT COEFFICIENT ################
    ###################################################################  
    def CL_Tail(CM_due_to_CL, CM_due_to_CD, CM_due_to_CT, CM_CG, l_ht, MAC,\
                xAC, xCG, alfa):
        """
        % Horizontal tail lift coefficient = CL_Tail(CM_due_to_CL, CM_due_to_CD,\
                                                     CM_due_to_CT, CM_CG, l_tail,\
                                                     MAC, xAC, xCG, alfa)
        % 
        %  This function is able to calculate the horizontal tail lift
        %  coefficient associated with the prescribed flight condition,
        %  defined in terms of airspeed V and load factor n. This lift
        %  coefficient is essential to correctly size the horizontal
        %  tail structure and to get a knowledge about the extra lift
        %  that the main lifting surface must produce to achieve
        %  steady, non - accelerated, levelled flight to compensate for
        %  the airloads acting on the horizontal tail. 
        %  
        %  INPUTS 
        %  CM_due_to_CL --> Pitching mom. contribution associated with
        %                   wing - body lift
        %  CM_due_to_CD --> Pitching mom. contribution associated with
        %                   wing - body drag
        %  CM_due_to_CT --> Pitching mom. associated with non baricentral
        %                   thrust 
        %  CM_CG        --> Pitching mom. of the wing body configuration 
        %  l_ht         --> Horizontal tail mom. arm 
        %  MAC          --> Mean Aerodynamic chord
        %  xAC          --> Aerodynamic center position in unit of mean
        %                   aerodynamic chord
        %  xCG          --> Center of gravity position in unit of mean
        %                   aerodynamic chord
        %  alfa         --> Angle of attack at the prescribed flight
        %                   condition
        %
        %  OUTPUTS
        %  CL_Tail   ---> Tail lift coefficient
        %
        %  The function calculates the CL_Tail with the following
        %  equilibrium relationship:
        %
        %             (CLwb_contr + CDwb_contr + CT_contr + CM)
        %  CL_Tail =  -----------------------------------------
        %                l_tail    
        %               -------- + (xac - xcg) * cos(alpha)
        %                  MAC
        %
        """
        
        temp1   = CM_CG + CM_due_to_CL + CM_due_to_CD + CM_due_to_CT
        temp2   = l_ht / MAC
        temp3   = xAC - xCG
        
        CL_tail = ( temp1 ) / ( temp2 + temp3 * np.cos(alfa) )
        
        return CL_tail
        
    #########################################################################
    ########### FLIGHT BALANCING PARAMETERS CALCULATIONS FUNCTION ###########
    #########################################################################
    def flight_balancing_parameters(rho, V, WS, n, CLmax_aero_model, alfa0l,\
                                    p_drag, CL_star, p_lift1, p_lift2, xAC,\
                                    xCG, bCG, MAC, h, CM0, CM_gear, S, l_ht,\
                                    flag1, flag2, CM_CL, obj1): 
        """
        Function that calculates parameters related to balancing loads. The 
        aircraft is considered in equilibrium during normal flight operation.

        Parameters
        ----------
        rho : float
            Air density.
        V : float
            Airspeed from the flight envelope.
        WS : float
            Wing loading.
        n : float
            Load factor from the flight envelope.
        CLmax_aero_model : TYPE
            Max lift coefficient from the aerodynamic model.
        alfa0l : float
            Zero lift angle of attack from the aerodynamic model.
        p_drag : array of float
            Interpolation coefficients relative to the drag curve. The variable
            is an array  containing Chebyshev coefficients ordered  from low to 
            high. The array is structured as follow: 
                                   [c_0, c_1, ... , c_n]
            with n the degree(s) of the fitting polynomials.
        CL_star : float
            Lift coefficient at the end of the linear part of the lift curve.
        CL0 : float
            Lift coefficient at zero angle of attack.
        CLalfa : float
            Lift curve slope expressed in [1/deg].
        p_lift1 : array of float
            Interpolation coefficients relative to the lift curve. The variable
            is an array  containing Chebyshev coefficients ordered  from low to 
            high. The array is structured as follow: 
                                   [c_0, c_1, ... , c_n]
            with n the degree(s) of the fitting polynomials. This vector cont-
            ains the values relative to the linear part of the lift curve.
        p_lift2 : array of float
            Interpolation coefficients relative to the lift curve. The variable
            is an array  containing Chebyshev coefficients ordered  from low to 
            high. The array is structured as follow: 
                                   [c_0, c_1, ... , c_n]
            with n the degree(s) of the fitting polynomials. This vector cont-
            ains the values relative to the non-linear part of the lift curve.
        CL_max_inv : float
            Maximum lift coefficient during inverted flight ("negative" lift).
        xAC : float
            Aerodynamic centre abscissa in [MAC percentage].
        xCG : float
            Centre of gravity abscissa in [MAC percentage].
        bCG : float
            Dimensional, vertical arm of the lift force horizontal component.
            This must be divided by the MAC before to apply to the equation for
            the pitching moment equilibrium. In this instance, bCG is expressed
            in [m].
        MAC : float
            Mean Aerodynamic Chord in [m].
        h : float
            Moment arm measured in meters of the thrust action line with 
            respect to the centre of gravity. This value must be divided by the
            Mean Aerodynamic Chord.
        CM0 : float
            Pitching moment coefficient at zero angle of attack.
        CM_gear : float
            Pitching moment contribute corresponding to the landing gear.
        S : float
            Main wing surface expresse in [m^2].
        l_ht : float
            Horizontal tail moment arm.
        flag1 : string 
            Flag to identify flight condition. Values: 
                --> "straight" 
                --> "inverted"
        flag2 : string 
            Flag to identify angles unit of measure. Values: 
                --> "deg" 
                --> "rad"
        CM_CL : float
            Pitching moment curve slope in terms of lift coefficient.
        obj1 : object - balancing_loads
            Call to a class of methods for balancing loads calculations.
            
        Returns
        -------
        CL_wb : float
            Wing body lift coefficient.
        alfa : float
            Angle of attack.
        CD : float
            Wing body drag coefficient.
        q : float
            Dynamic pressure.
        CM_CL : float
            Pitching moment coefficient contribute due to lift.
        CM_CD : float
            Pitching moment coefficient contribute due to drag.  
        CM_CT : float
            Pitching moment coefficient contribute due to thrust.
        CM_CG : float
            Total pitching moment coefficient; this coefficient is globally 
            evaluated with respect to the aircraft centre of gravity.
        CL_ht_new : float
            Horizontal tail lift coefficient.
        CL_new : float
            Full-vehicle (Wing + Body + Horiz. tail) lift coefficient.
        L_wb : float
            Wing body lift expressed in [daN].
        L_new : float
            Full-vehicle lift expressed in [daN].  
        L_ht : float
            Horizontal tail lift expressed in [daN].  
            
        """
        # INITIALIZATION
        global alfa_new
        alfa_new = 0.0
        
        CL_wb = obj1.CL_max_function(rho, V, WS, n, CLmax_aero_model)
        
        if ( abs(CL_wb) >= abs(CLmax_aero_model) ):
            
            CL_wb        = CLmax_aero_model
            alfa         = obj1.alpha_fullmodel(p_lift1, p_lift2, CL_wb,\
                                                CLmax_aero_model, CL_star, alfa0l,\
                                                flag1)
            CD           = chebval(alfa, p_drag)
            q            = 0.5 * rho * V**2
            CM_due_to_CL = obj1.CLWB_pitch_contrib(CL_wb, alfa, xAC, xCG, bCG, MAC, flag2)
            CM_due_to_CD = obj1.CDWB_pitch_contrib(CD, alfa, xAC, xCG, bCG, MAC, flag2)
            CM_due_to_CT = obj1.CTWB_pitch_contrib(CD, h, MAC)
            CM_CG        = obj1.CM_aboutcg(CM0, CM_gear, CM_CL, CL_wb)
            CL_ht        = obj1.CL_Tail(CM_due_to_CL, CM_due_to_CD, CM_due_to_CT,\
                                        CM_CG, l_ht, MAC, xAC, xCG, alfa)
            
            # FULL-VEHICLE LIFT COEFFICIENT
            # wb  == wing-body
            # new == full-vehicle
            # ht  == horizontal-tail 
            CL_new = CL_wb - CL_ht 
            tol    = 1E-3
            k      = 1
            while ( np.abs(CL_wb - CL_new) > tol ):
                alfa_new = obj1.alpha_fullmodel(p_lift1, p_lift2, CL_new, CLmax_aero_model,\
                                                    CL_star, alfa0l, flag1)
                CD_new       = chebval(alfa_new, p_drag)
                CM_due_to_CL = obj1.CLWB_pitch_contrib(CL_new, alfa, xAC, xCG, bCG, MAC, flag2)
                CM_due_to_CD = obj1.CDWB_pitch_contrib(CD_new, alfa, xAC, xCG, bCG, MAC, flag2)
                CM_due_to_CT = obj1.CTWB_pitch_contrib(CD_new, h, MAC)
                CM_CG        = obj1.CM_aboutcg(CM0, CM_gear, CM_CL, CL_new)
                CL_ht_new    = obj1.CL_Tail(CM_due_to_CL, CM_due_to_CD, CM_due_to_CT,\
                                            CM_CG, l_ht, MAC, xAC, xCG, alfa_new)
                # UPDATING LIFT COEFFICIENTS
                CL_ht  = CL_ht_new
                CL_new = CL_wb - CL_ht
                CL_wb  = CL_wb + ( CL_new - CL_wb ) * 1E-1
                
                # UPDATING COUNTER 
                k = k + 1 
                # CHECK THE COUNTER TO BREAK THE LOOP
                if ( k == 20 ): 
                    break
                
            # LIFT CALCULATIONS IN [daN]
            L_wb  = q * S * CL_wb * 1E-1
            L_new = q * S * CL_new * 1E-1
            L_ht  = q * S * CL_ht * 1E-1
        
        elif ( abs(CL_wb) <= abs(CLmax_aero_model) ):
            
            alfa         = obj1.alpha_fullmodel(p_lift1, p_lift2, CL_wb,\
                                                CLmax_aero_model, CL_star, alfa0l,\
                                                flag1)
            CD           = chebval(alfa, p_drag)
            q            = 0.5 * rho * V**2
            CM_due_to_CL = obj1.CLWB_pitch_contrib(CL_wb, alfa, xAC, xCG, bCG, MAC, flag2)
            CM_due_to_CD = obj1.CDWB_pitch_contrib(CD, alfa, xAC, xCG, bCG, MAC, flag2)
            CM_due_to_CT = obj1.CTWB_pitch_contrib(CD, h, MAC)
            CM_CG        = obj1.CM_aboutcg(CM0, CM_gear, CM_CL, CL_wb)
            CL_ht        = obj1.CL_Tail(CM_due_to_CL, CM_due_to_CD, CM_due_to_CT,\
                                        CM_CG, l_ht, MAC, xAC, xCG, alfa)
            
            # FULL-VEHICLE LIFT COEFFICIENT
            # wb  == wing-body
            # new == full-vehicle
            # ht  == horizontal-tail 
            CL_new = CL_wb - CL_ht 
            tol    = 1E-3
            k      = 1
            while ( np.abs(CL_wb - CL_new) > tol ):
                alfa_new     = obj1.alpha_fullmodel(p_lift1, p_lift2, CL_new, CLmax_aero_model,\
                                                    CL_star, alfa0l, flag1)
                CD_new       = chebval(alfa_new, p_drag)
                CM_due_to_CL = obj1.CLWB_pitch_contrib(CL_new, alfa, xAC, xCG, bCG, MAC, flag2)
                CM_due_to_CD = obj1.CDWB_pitch_contrib(CD_new, alfa, xAC, xCG, bCG, MAC, flag2)
                CM_due_to_CT = obj1.CTWB_pitch_contrib(CD_new, h, MAC)
                CM_CG        = obj1.CM_aboutcg(CM0, CM_gear, CM_CL, CL_new)
                CL_ht_new    = obj1.CL_Tail(CM_due_to_CL, CM_due_to_CD, CM_due_to_CT,\
                                            CM_CG, l_ht, MAC, xAC, xCG, alfa_new)
                # UPDATING LIFT COEFFICIENTS
                CL_ht  = CL_ht_new
                CL_new = CL_wb - CL_ht
                CL_wb  = CL_wb + ( CL_new - CL_wb ) * 1E-1
                
                # UPDATING COUNTER 
                k = k + 1 
                # CHECK THE COUNTER TO BREAK THE LOOP
                if ( k == 20 ): 
                    break
                
            # LIFT CALCULATIONS IN [daN]
            L_wb  = q * S * CL_wb * 1E-1
            L_new = q * S * CL_new * 1E-1
            L_ht  = q * S * CL_ht * 1E-1
            
        return CL_wb, alfa, alfa_new, CD, q, CM_due_to_CL, CM_due_to_CD, CM_due_to_CT,\
               CM_CG, CL_ht, CL_new, L_wb, L_new, L_ht       
               
###############################################################################
################## MAIN WING LIFT DIAGRAM - BALANCING LOADS ###################
###############################################################################
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # V - n DIAGRAM 
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def main_wing_lift(V_fe_from0toS, V_fe_fromStoA, V_fe_fromAtoGust1,\
                       V_fe_fromGust1toC, V_fe_fromCtoGust2, V_fe_fromGust2toD,\
                       V_fe_fromDto0, V_fe_from0toSinv, V_fe_fromSinvtoG,\
                       V_fe_fromGtoGust1, V_fe_fromGust1toF, V_fe_fromFtoE,\
                       V_fe_fromEto0, V_fe_fromGtoF, V_fe_fromAtoC,\
                       V_fe_fromCtoD, L_wb_from0toS, L_new_from0toS, L_wb_fromStoA,\
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
                       L_ht_fromGtoF, L_ht_fromAtoC, L_ht_fromCtoD, L_ht_unit_load_factor,\
                       V_unit_load_factor, Aircraft_name, n, pos_case_flag, neg_case_flag):
        """
        Function that plots main wing lift diagrams, both wing-body and full-
        vehicle. 
        wb == wing-body
        fv == full-vehicle

        Parameters
        ----------
        L_wb_from0toS : TYPE
            DESCRIPTION.
        L_new_from0toS : TYPE
            DESCRIPTION.
        L_wb_fromStoA : TYPE
            DESCRIPTION.
        L_new_fromStoA : TYPE
            DESCRIPTION.
        L_new_fromAtoGust1 : TYPE
            DESCRIPTION.
        L_wb_fromAtoGust1 : TYPE
            DESCRIPTION.
        L_new_fromGust1toC : TYPE
            DESCRIPTION.
        L_wb_fromGust1toC : TYPE
            DESCRIPTION.
        L_new_fromCtoGust2 : TYPE
            DESCRIPTION.
        L_wb_fromCtoGust2 : TYPE
            DESCRIPTION.
        L_new_fromGust2toD : TYPE
            DESCRIPTION.
        L_wb_fromGust2toD : TYPE
            DESCRIPTION.
        L_new_fromDto0 : TYPE
            DESCRIPTION.
        L_wb_fromDto0 : TYPE
            DESCRIPTION.
        L_new_from0toSinv : TYPE
            DESCRIPTION.
        L_wb_from0toSinv : TYPE
            DESCRIPTION.
        L_new_fromSinvtoG : TYPE
            DESCRIPTION.
        L_wb_fromSinvtoG : TYPE
            DESCRIPTION.
        L_new_fromGtoGust1 : TYPE
            DESCRIPTION.
        L_wb_fromGtoGust1 : TYPE
            DESCRIPTION.
        L_new_fromGust1toF : TYPE
            DESCRIPTION.
        L_wb_fromGust1toF : TYPE
            DESCRIPTION.
        L_new_fromFtoE : TYPE
            DESCRIPTION.
        L_wb_fromFtoE : TYPE
            DESCRIPTION.
        L_new_fromEto0 : TYPE
            DESCRIPTION.
        L_wb_fromEto0 : TYPE
            DESCRIPTION.
        L_new_fromGtoF : TYPE
            DESCRIPTION.
        L_wb_fromGtoF : TYPE
            DESCRIPTION.
        L_new_fromAtoC : TYPE
            DESCRIPTION.
        V_fe_fromAtoC : TYPE
            DESCRIPTION.
        L_new_fromCtoD : TYPE
            DESCRIPTION.
        V_fe_fromCtoD : TYPE
            DESCRIPTION.
        Reg : TYPE
            DESCRIPTION.
        Aircraft_name : TYPE
            DESCRIPTION.
        n : TYPE
            DESCRIPTION.
        pos_case_flag : TYPE
            DESCRIPTION.
        neg_case_flag : TYPE
            DESCRIPTION.

        Returns
        -------
        fig : TYPE
            DESCRIPTION.

        """
        
        # ===================================================================
        # ===================================================================
        rc('font',**{'family':'serif','serif':['Palatino']})
        rc('text', usetex=True)
        # ===================================================================
        if   (pos_case_flag == "Case1") and (neg_case_flag == "Case1_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_new_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_new_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_new_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_new_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_new_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_new_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_new_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_new_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_new_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_new_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()   
                        
            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_wb_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_wb_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_wb_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_wb_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_wb_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_wb_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_wb_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_wb_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_wb_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_wb_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()  
                        
            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_ht_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_ht_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_ht_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_ht_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_ht_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_ht_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_ht_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_ht_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_ht_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_ht_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()  
            
        elif (pos_case_flag == "Case1") and (neg_case_flag == "Case2_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_new_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_new_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_new_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_new_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_new_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_new_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_new_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()     
                
            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_wb_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_wb_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_wb_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_wb_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_wb_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_wb_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_wb_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show() 
                
            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_ht_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_ht_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_ht_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_ht_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_ht_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_ht_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_ht_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()                
            
        elif (pos_case_flag == "Case1") and (neg_case_flag == "Case3_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_new_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_new_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_new_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_new_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_new_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_new_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_new_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()  
              
            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_wb_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_wb_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_wb_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_wb_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_wb_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_wb_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_wb_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()  
              
            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_ht_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_ht_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_ht_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_ht_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_ht_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_ht_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_ht_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2) 

            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()               
            
        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case1_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_new_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_new_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_new_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_new_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_new_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_new_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_new_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_new_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_new_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_new_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()      

            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_wb_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_wb_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_wb_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_wb_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_wb_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_wb_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_wb_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_wb_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_wb_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_wb_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()   

            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_ht_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_ht_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_ht_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_ht_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_ht_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_ht_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_ht_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_ht_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_ht_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_ht_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2) 

            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()         
            
        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case2_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_new_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_new_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_new_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_new_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_new_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_new_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_new_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()       
            
            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_wb_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_wb_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_wb_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_wb_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_wb_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_wb_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_wb_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()   
            
            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_ht_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_ht_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_ht_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_ht_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_ht_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_ht_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_ht_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)    

            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()            
            
        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case3_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_new_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_new_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_new_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_new_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_new_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_new_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_new_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()             
            
            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_wb_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_wb_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_wb_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_wb_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_wb_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_wb_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_wb_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()  
            
            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     L_ht_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     L_ht_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_ht_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], L_ht_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], L_ht_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, L_ht_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], L_ht_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)    

            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()    
            
        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case1_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_new_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_new_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_new_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_new_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_new_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_new_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_new_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_new_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_new_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_new_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()     
            
            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_wb_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_wb_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_wb_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_wb_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_wb_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_wb_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_wb_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_wb_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_wb_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_wb_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()
            
            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_ht_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_ht_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_ht_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_ht_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, L_ht_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, L_ht_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_ht_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_ht_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], L_ht_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], L_ht_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)

            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()   

        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case2_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_new_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_new_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_new_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_new_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_new_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_new_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_new_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()                
            
            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_wb_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_wb_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_wb_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_wb_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_wb_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_wb_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_wb_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()           
            
            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_ht_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_ht_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_ht_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_ht_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_ht_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_ht_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_ht_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)             

            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()   
            
        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case3_inverted"):
            
            fig1  = plt.figure()
            plt.plot(V_fe_from0toS,     L_new_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_new_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_new_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_new_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_new_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_new_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_new_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_new_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_new_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_new_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_new_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_new_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_new_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_new_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_new_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_new_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_new_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_new_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_new_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_new_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_new_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_new_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_new_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_new_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{fv}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Full vehicle lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFullVehicle' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()       
                
            fig2  = plt.figure()
            plt.plot(V_fe_from0toS,     L_wb_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     L_wb_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_wb_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_wb_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_wb_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_wb_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_wb_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_wb_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_wb_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_wb_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_wb_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_wb_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_wb_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_wb_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_wb_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_wb_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_wb_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_wb_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_wb_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_wb_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_wb_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_wb_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_wb_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_wb_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Main wing lift} ~ $L_{wb}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Wing body lift on the main wing' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'MainwingLiftFWingBody' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()   
                
            fig3  = plt.figure()
            plt.plot(V_fe_from0toS,     L_ht_from0toS,     color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromStoA,     L_ht_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, L_ht_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, L_ht_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, L_ht_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, L_ht_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     L_ht_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            plt.plot(V_unit_load_factor,L_ht_unit_load_factor, color="black", linewidth=1.8, linestyle = "dashed", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  L_ht_from0toSinv,  color="red", linewidth=0.2, linestyle = "dashed", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  L_ht_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     L_ht_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     L_ht_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     L_ht_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], L_ht_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], L_ht_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], L_ht_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], L_ht_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], L_ht_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, L_ht_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], L_ht_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], L_ht_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], L_ht_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], L_ht_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, L_ht_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], L_ht_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)        

            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Horiz. tailplane lift} ~ $L_{ht}$ ~ $[daN]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Horizontal tailplane lift' + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'HorizontalTailplaneLift' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()             
        
        return fig1, fig2, fig3      

    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # BALANCING VALUES STORED INSIDE THE SIMPLE NAMESPACE OBJECT
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def balancing_values_storage(          
                       Balancing_loads,\
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
                       pos_case_flag, neg_case_flag):   
        
        if   (pos_case_flag == "Case1") and (neg_case_flag == "Case1_inverted"):
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS"    : {"Value": V_fe_from0toS,     "Unit": "m/s"},
                    "n_fe_from0toS"    : {"Value": n_fe_from0toS,     "Unit": "g"},
                    "V_fe_fromStoA"    : {"Value": V_fe_fromStoA,     "Unit": "m/s"},
                    "n_fe_fromStoA"    : {"Value": n_fe_fromStoA,     "Unit": "g"},
                    "V_fe_fromAtoGust1": {"Value": V_fe_fromAtoGust1, "Unit": "m/s"},
                    "n_fe_fromAtoGust1": {"Value": n_fe_fromAtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toC": {"Value": V_fe_fromGust1toC, "Unit": "m/s"},
                    "n_fe_fromGust1toC": {"Value": n_fe_fromGust1toC, "Unit": "g"},
                    "V_fe_fromCtoGust2": {"Value": V_fe_fromCtoGust2, "Unit": "m/s"},
                    "n_fe_fromCtoGust2": {"Value": n_fe_fromCtoGust2, "Unit": "g"},
                    "V_fe_fromGust2toD": {"Value": V_fe_fromGust2toD, "Unit": "m/s"},
                    "n_fe_fromGust2toD": {"Value": n_fe_fromGust2toD, "Unit": "g"},
                    "V_fe_fromDto0"    : {"Value": V_fe_fromDto0,     "Unit": "m/s"},
                    "n_fe_fromDto0"    : {"Value": n_fe_fromDto0,     "Unit": "g"},
                    "V_fe_from0toSinv" : {"Value": V_fe_from0toSinv,  "Unit": "m/s"},
                    "n_fe_from0toSinv" : {"Value": n_fe_from0toSinv,  "Unit": "g"},
                    "V_fe_fromSinvtoG" : {"Value": V_fe_fromSinvtoG,  "Unit": "m/s"},
                    "n_fe_fromSinvtoG" : {"Value": n_fe_fromSinvtoG,  "Unit": "g"},
                    "V_fe_fromGtoGust1": {"Value": V_fe_fromGtoGust1, "Unit": "m/s"},
                    "n_fe_fromGtoGust1": {"Value": n_fe_fromGtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toF": {"Value": V_fe_fromGust1toF, "Unit": "m/s"},
                    "n_fe_fromGust1toF": {"Value": n_fe_fromGust1toF, "Unit": "g"},
                    "V_fe_fromFtoE"    : {"Value": V_fe_fromFtoE,     "Unit": "m/s"},
                    "n_fe_fromFtoE"    : {"Value": n_fe_fromFtoE,     "Unit": "g"},
                    "V_fe_fromEto0"    : {"Value": V_fe_fromEto0,     "Unit": "m/s"},
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0,     "Unit": "g"},
                      
                    "q_from0toS"     : {"Value": q_from0toS,      "Unit": "Pa"},
                    "CD_from0toS"    : {"Value": CD_from0toS,     "Unit": "Non dimensional"},
                    "q_fromStoA"     : {"Value": q_fromStoA,      "Unit": "Pa"},
                    "CD_fromStoA"    : {"Value": CD_fromStoA,     "Unit": "Non dimensional"},
                    "q_fromAtoGust1" : {"Value": q_fromAtoGust1,  "Unit": "Pa"},
                    "CD_fromAtoGust1": {"Value": CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toC" : {"Value": q_fromGust1toC,  "Unit": "Pa"},
                    "CD_fromGust1toC": {"Value": CD_fromGust1toC, "Unit": "Non dimensional"},
                    "q_fromCtoGust2" : {"Value": q_fromCtoGust2,  "Unit": "Pa"},
                    "CD_fromCtoGust2": {"Value": CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "q_fromGust2toD" : {"Value": q_fromGust2toD,  "Unit": "Pa"},
                    "CD_fromGust2toD": {"Value": CD_fromGust2toD, "Unit": "Non dimensional"},
                    "q_fromDto0"     : {"Value": q_fromDto0,      "Unit": "Pa"},
                    "CD_fromDto0"    : {"Value": CD_fromDto0,     "Unit": "Non dimensional"},
                    "q_from0toSinv"  : {"Value": q_from0toSinv,   "Unit": "Pa"},
                    "CD_from0toSinv" : {"Value": CD_from0toSinv,  "Unit": "Non dimensional"},
                    "q_fromSinvtoG"  : {"Value": q_fromSinvtoG,   "Unit": "Pa"},
                    "CD_fromSinvtoG" : {"Value": CD_fromSinvtoG,  "Unit": "Non dimensional"},
                    "q_fromGtoGust1" : {"Value": q_fromGtoGust1,  "Unit": "Pa"},
                    "CD_fromGtoGust1": {"Value": CD_fromGtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toF" : {"Value": q_fromGust1toF,  "Unit": "Pa"},
                    "CD_fromGust1toF": {"Value": CD_fromGust1toF, "Unit": "Non dimensional"},
                    "q_fromFtoE"     : {"Value": q_fromFtoE,      "Unit": "Pa"},
                    "CD_fromFtoE"    : {"Value": CD_fromFtoE,     "Unit": "Non dimensional"},
                    "q_fromEto0"     : {"Value": q_fromEto0,      "Unit": "Pa"},
                    "CD_fromEto0"    : {"Value": CD_fromEto0,     "Unit": "Non dimensional"},
                      
                    "CL_wb_from0toS"     : {"Value": CL_wb_from0toS,      "Unit": "Non dimensional"},
                    "CL_ht_from0toS"     : {"Value": CL_ht_from0toS,      "Unit": "Non dimensional"},
                    "CL_new_from0toS"    : {"Value": CL_new_from0toS,     "Unit": "Non dimensional"},
                    "CL_wb_fromStoA"     : {"Value": CL_wb_fromStoA,      "Unit": "Non dimensional"},
                    "CL_ht_fromStoA"     : {"Value": CL_ht_fromStoA,      "Unit": "Non dimensional"},
                    "CL_new_fromStoA"    : {"Value": CL_new_fromStoA,     "Unit": "Non dimensional"},
                    "CL_wb_fromAtoGust1" : {"Value": CL_wb_fromAtoGust1,  "Unit": "Non dimensional"},
                    "CL_ht_fromAtoGust1" : {"Value": CL_ht_fromAtoGust1,  "Unit": "Non dimensional"},
                    "CL_new_fromAtoGust1": {"Value": CL_new_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_wb_fromGust1toC" : {"Value": CL_wb_fromGust1toC,  "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toC" : {"Value": CL_ht_fromGust1toC,  "Unit": "Non dimensional"},
                    "CL_new_fromGust1toC": {"Value": CL_new_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_wb_fromCtoGust2" : {"Value": CL_wb_fromCtoGust2,  "Unit": "Non dimensional"},
                    "CL_ht_fromCtoGust2" : {"Value": CL_ht_fromCtoGust2,  "Unit": "Non dimensional"},
                    "CL_new_fromCtoGust2": {"Value": CL_new_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_wb_fromGust2toD" : {"Value": CL_wb_fromGust2toD,  "Unit": "Non dimensional"},
                    "CL_ht_fromGust2toD" : {"Value": CL_ht_fromGust2toD,  "Unit": "Non dimensional"},
                    "CL_new_fromGust2toD": {"Value": CL_new_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_wb_fromDto0"     : {"Value": CL_wb_fromDto0,      "Unit": "Non dimensional"},
                    "CL_ht_fromDto0"     : {"Value": CL_ht_fromDto0,      "Unit": "Non dimensional"},
                    "CL_new_fromDto0"    : {"Value": CL_new_fromDto0,     "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv"  : {"Value": CL_wb_from0toSinv,   "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv"  : {"Value": CL_ht_from0toSinv,   "Unit": "Non dimensional"},
                    "CL_new_from0toSinv" : {"Value": CL_new_from0toSinv,  "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG"  : {"Value": CL_wb_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG"  : {"Value": CL_ht_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG" : {"Value": CL_new_fromSinvtoG,  "Unit": "Non dimensional"},
                    "CL_wb_fromGtoGust1" : {"Value": CL_wb_fromGtoGust1,  "Unit": "Non dimensional"},
                    "CL_ht_fromGtoGust1" : {"Value": CL_ht_fromGtoGust1,  "Unit": "Non dimensional"},
                    "CL_new_fromGtoGust1": {"Value": CL_new_fromGtoGust1, "Unit": "Non dimensional"},
                    "CL_wb_fromGust1toF" : {"Value": CL_wb_fromGust1toF,  "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toF" : {"Value": CL_ht_fromGust1toF,  "Unit": "Non dimensional"},
                    "CL_new_fromGust1toF": {"Value": CL_new_fromGust1toF, "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE"     : {"Value": CL_wb_fromFtoE,      "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE"     : {"Value": CL_ht_fromFtoE,      "Unit": "Non dimensional"},
                    "CL_new_fromFtoE"    : {"Value": CL_new_fromFtoE,     "Unit": "Non dimensional"},
                    "CL_wb_fromEto0"     : {"Value": CL_wb_fromEto0,      "Unit": "Non dimensional"},
                    "CL_ht_fromEto0"     : {"Value": CL_ht_fromEto0,      "Unit": "Non dimensional"},
                    "CL_new_fromEto0"    : {"Value": CL_new_fromEto0,     "Unit": "Non dimensional"},

                    "alfa_from0toS"        : {"Value": alfa_from0toS,         "Unit": "deg"},
                    "alfa_new_from0toS"    : {"Value": alfa_new_from0toS,     "Unit": "deg"},
                    "alfa_fromStoA"        : {"Value": alfa_fromStoA,         "Unit": "deg"},
                    "alfa_new_fromStoA"    : {"Value": alfa_new_fromStoA,     "Unit": "deg"},
                    "alfa_fromAtoGust1"    : {"Value": alfa_fromAtoGust1,     "Unit": "deg"},
                    "alfa_new_fromAtoGust1": {"Value": alfa_new_fromAtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toC"    : {"Value": alfa_fromGust1toC,     "Unit": "deg"},
                    "alfa_new_fromGust1toC": {"Value": alfa_new_fromGust1toC, "Unit": "deg"},
                    "alfa_fromCtoGust2"    : {"Value": alfa_fromCtoGust2,     "Unit": "deg"},
                    "alfa_new_fromCtoGust2": {"Value": alfa_new_fromCtoGust2, "Unit": "deg"},
                    "alfa_fromGust2toD"    : {"Value": alfa_fromGust2toD,     "Unit": "deg"},
                    "alfa_new_fromGust2toD": {"Value": alfa_new_fromGust2toD, "Unit": "deg"},
                    "alfa_fromDto0"        : {"Value": alfa_fromDto0,         "Unit": "deg"},
                    "alfa_new_fromDto0"    : {"Value": alfa_new_fromDto0,     "Unit": "deg"},
                    "alfa_from0toSinv"     : {"Value": alfa_from0toSinv,      "Unit": "deg"},
                    "alfa_new_from0toSinv" : {"Value": alfa_new_from0toSinv,  "Unit": "deg"},
                    "alfa_fromSinvtoG"     : {"Value": alfa_fromSinvtoG,      "Unit": "deg"},
                    "alfa_new_fromSinvtoG" : {"Value": alfa_new_fromSinvtoG,  "Unit": "deg"},
                    "alfa_fromGtoGust1"    : {"Value": alfa_fromGtoGust1,     "Unit": "deg"},
                    "alfa_new_fromGtoGust1": {"Value": alfa_new_fromGtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toF"    : {"Value": alfa_fromGust1toF,     "Unit": "deg"},
                    "alfa_new_fromGust1toF": {"Value": alfa_new_fromGust1toF, "Unit": "deg"},
                    "alfa_fromFtoE"        : {"Value": alfa_fromFtoE,         "Unit": "deg"},
                    "alfa_new_fromFtoE"    : {"Value": alfa_new_fromFtoE,     "Unit": "deg"},
                    "alfa_fromEto0"        : {"Value": alfa_fromEto0,         "Unit": "deg"},
                    "alfa_new_fromEto0"    : {"Value": alfa_new_fromEto0,     "Unit": "deg"},
                      
                    "L_wb_from0toS"     : {"Value": L_wb_from0toS,      "Unit": "daN"},
                    "L_ht_from0toS"     : {"Value": L_ht_from0toS,      "Unit": "daN"},
                    "L_new_from0toS"    : {"Value": L_new_from0toS,     "Unit": "daN"},
                    "L_wb_fromStoA"     : {"Value": L_wb_fromStoA,      "Unit": "daN"},
                    "L_ht_fromStoA"     : {"Value": L_ht_fromStoA,      "Unit": "daN"},
                    "L_new_fromStoA"    : {"Value": L_new_fromStoA,     "Unit": "daN"},
                    "L_wb_fromAtoGust1" : {"Value": L_wb_fromAtoGust1,  "Unit": "daN"},
                    "L_ht_fromAtoGust1" : {"Value": L_ht_fromAtoGust1,  "Unit": "daN"},
                    "L_new_fromAtoGust1": {"Value": L_new_fromAtoGust1, "Unit": "daN"},
                    "L_wb_fromGust1toC" : {"Value": L_wb_fromGust1toC,  "Unit": "daN"},
                    "L_ht_fromGust1toC" : {"Value": L_ht_fromGust1toC,  "Unit": "daN"},
                    "L_new_fromGust1toC": {"Value": L_new_fromGust1toC, "Unit": "daN"},
                    "L_wb_fromCtoGust2" : {"Value": L_wb_fromCtoGust2,  "Unit": "daN"},
                    "L_ht_fromCtoGust2" : {"Value": L_ht_fromCtoGust2,  "Unit": "daN"},
                    "L_new_fromCtoGust2": {"Value": L_new_fromCtoGust2, "Unit": "daN"},
                    "L_wb_fromGust2toD" : {"Value": L_wb_fromGust2toD,  "Unit": "daN"},
                    "L_ht_fromGust2toD" : {"Value": L_ht_fromGust2toD,  "Unit": "daN"},
                    "L_new_fromGust2toD": {"Value": L_new_fromGust2toD, "Unit": "daN"},
                    "L_wb_fromDto0"     : {"Value": L_wb_fromDto0,      "Unit": "daN"},
                    "L_ht_fromDto0"     : {"Value": L_ht_fromDto0,      "Unit": "daN"},
                    "L_new_fromDto0"    : {"Value": L_new_fromDto0,     "Unit": "daN"},
                    "L_wb_from0toSinv"  : {"Value": L_wb_from0toSinv,   "Unit": "daN"},
                    "L_ht_from0toSinv"  : {"Value": L_ht_from0toSinv,   "Unit": "daN"},
                    "L_new_from0toSinv" : {"Value": L_new_from0toSinv,  "Unit": "daN"},
                    "L_wb_fromSinvtoG"  : {"Value": L_wb_fromSinvtoG,   "Unit": "daN"},
                    "L_ht_fromSinvtoG"  : {"Value": L_ht_fromSinvtoG,   "Unit": "daN"},
                    "L_new_fromSinvtoG" : {"Value": L_new_fromSinvtoG,  "Unit": "daN"},
                    "L_wb_fromGtoGust1" : {"Value": L_wb_fromGtoGust1,  "Unit": "daN"},
                    "L_ht_fromGtoGust1" : {"Value": L_ht_fromGtoGust1,  "Unit": "daN"},
                    "L_new_fromGtoGust1": {"Value": L_new_fromGtoGust1, "Unit": "daN"},
                    "L_wb_fromGust1toF" : {"Value": L_wb_fromGust1toF,  "Unit": "daN"},
                    "L_ht_fromGust1toF" : {"Value": L_ht_fromGust1toF,  "Unit": "daN"},
                    "L_new_fromGust1toF": {"Value": L_new_fromGust1toF, "Unit": "daN"},
                    "L_wb_fromFtoE"     : {"Value": L_wb_fromFtoE,      "Unit": "daN"},
                    "L_ht_fromFtoE"     : {"Value": L_ht_fromFtoE,      "Unit": "daN"},
                    "L_new_fromFtoE"    : {"Value": L_new_fromFtoE,     "Unit": "daN"},
                    "L_wb_fromEto0"     : {"Value": L_wb_fromEto0,      "Unit": "daN"},
                    "L_ht_fromEto0"     : {"Value": L_ht_fromEto0,      "Unit": "daN"},
                    "L_new_fromEto0"    : {"Value": L_new_fromEto0,     "Unit": "daN"},
                      
                    "CM_due_to_CL_from0toS"     : {"Value": CM_due_to_CL_from0toS,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS"     : {"Value": CM_due_to_CD_from0toS,      "Unit": "Non dimensional"},
                    "CM_CG_from0toS"            : {"Value": CM_CG_from0toS,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA"     : {"Value": CM_due_to_CL_fromStoA,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA"     : {"Value": CM_due_to_CD_fromStoA,      "Unit": "Non dimensional"},
                    "CM_CG_fromStoA"            : {"Value": CM_CG_fromStoA,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoGust1" : {"Value": CM_due_to_CL_fromAtoGust1,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoGust1" : {"Value": CM_due_to_CD_fromAtoGust1,  "Unit": "Non dimensional"},
                    "CM_CG_fromAtoGust1"        : {"Value": CM_CG_fromAtoGust1,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toC" : {"Value": CM_due_to_CL_fromGust1toC,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toC" : {"Value": CM_due_to_CD_fromGust1toC,  "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toC"        : {"Value": CM_CG_fromGust1toC,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoGust2" : {"Value": CM_due_to_CL_fromCtoGust2,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoGust2" : {"Value": CM_due_to_CD_fromCtoGust2,  "Unit": "Non dimensional"},
                    "CM_CG_fromCtoGust2"        : {"Value": CM_CG_fromCtoGust2,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust2toD" : {"Value": CM_due_to_CL_fromGust2toD,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust2toD" : {"Value": CM_due_to_CD_fromGust2toD,  "Unit": "Non dimensional"},
                    "CM_CG_fromGust2toD"        : {"Value": CM_CG_fromGust2toD,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0"     : {"Value": CM_due_to_CL_fromDto0,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0"     : {"Value": CM_due_to_CD_fromDto0,      "Unit": "Non dimensional"},
                    "CM_CG_fromDto0"            : {"Value": CM_CG_fromDto0,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv"  : {"Value": CM_due_to_CL_from0toSinv,   "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv"  : {"Value": CM_due_to_CD_from0toSinv,   "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv"         : {"Value": CM_CG_from0toSinv,          "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG"  : {"Value": CM_due_to_CL_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG"  : {"Value": CM_due_to_CD_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG"         : {"Value": CM_CG_fromSinvtoG,          "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoGust1" : {"Value": CM_due_to_CL_fromGtoGust1,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoGust1" : {"Value": CM_due_to_CD_fromGtoGust1,  "Unit": "Non dimensional"},
                    "CM_CG_fromGtoGust1"        : {"Value": CM_CG_fromGtoGust1,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toF" : {"Value": CM_due_to_CL_fromGust1toF,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toF" : {"Value": CM_due_to_CD_fromGust1toF,  "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toF"        : {"Value": CM_CG_fromGust1toF,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE"     : {"Value": CM_due_to_CL_fromFtoE,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE"     : {"Value": CM_due_to_CD_fromFtoE,      "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE"            : {"Value": CM_CG_fromFtoE,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0"     : {"Value": CM_due_to_CL_fromEto0,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0"     : {"Value": CM_due_to_CD_fromEto0,      "Unit": "Non dimensional"},
                    "CM_CG_fromEto0"            : {"Value": CM_CG_fromEto0,             "Unit": "Non dimensional"}
                    }
                }

        elif (pos_case_flag == "Case1") and (neg_case_flag == "Case2_inverted"):
            
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS"    : {"Value": V_fe_from0toS,     "Unit": "m/s"},
                    "n_fe_from0toS"    : {"Value": n_fe_from0toS,     "Unit": "g"},
                    "V_fe_fromStoA"    : {"Value": V_fe_fromStoA,     "Unit": "m/s"},
                    "n_fe_fromStoA"    : {"Value": n_fe_fromStoA,     "Unit": "g"},
                    "V_fe_fromAtoGust1": {"Value": V_fe_fromAtoGust1, "Unit": "m/s"},
                    "n_fe_fromAtoGust1": {"Value": n_fe_fromAtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toC": {"Value": V_fe_fromGust1toC, "Unit": "m/s"},
                    "n_fe_fromGust1toC": {"Value": n_fe_fromGust1toC, "Unit": "g"},
                    "V_fe_fromCtoGust2": {"Value": V_fe_fromCtoGust2, "Unit": "m/s"},
                    "n_fe_fromCtoGust2": {"Value": n_fe_fromCtoGust2, "Unit": "g"},
                    "V_fe_fromGust2toD": {"Value": V_fe_fromGust2toD, "Unit": "m/s"},
                    "n_fe_fromGust2toD": {"Value": n_fe_fromGust2toD, "Unit": "g"},
                    "V_fe_fromDto0"    : {"Value": V_fe_fromDto0,     "Unit": "m/s"},
                    "n_fe_fromDto0"    : {"Value": n_fe_fromDto0,     "Unit": "g"},
                    "V_fe_from0toSinv" : {"Value": V_fe_from0toSinv,  "Unit": "m/s"},
                    "n_fe_from0toSinv" : {"Value": n_fe_from0toSinv,  "Unit": "g"},
                    "V_fe_fromSinvtoG" : {"Value": V_fe_fromSinvtoG,  "Unit": "m/s"},
                    "n_fe_fromSinvtoG" : {"Value": n_fe_fromSinvtoG,  "Unit": "g"},
                    "V_fe_fromGtoF"    : {"Value": V_fe_fromGtoF,     "Unit": "m/s"},
                    "n_fe_fromGtoF"    : {"Value": n_fe_fromGtoF,     "Unit": "g"},
                    "V_fe_fromFtoE"    : {"Value": V_fe_fromFtoE,     "Unit": "m/s"},
                    "n_fe_fromFtoE"    : {"Value": n_fe_fromFtoE,     "Unit": "g"},
                    "V_fe_fromEto0"    : {"Value": V_fe_fromEto0,     "Unit": "m/s"},
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0,     "Unit": "g"},

                    "q_from0toS"     : {"Value": q_from0toS,      "Unit": "Pa"},
                    "CD_from0toS"    : {"Value": CD_from0toS,     "Unit": "Non dimensional"},
                    "q_fromStoA"     : {"Value": q_fromStoA,      "Unit": "Pa"},
                    "CD_fromStoA"    : {"Value": CD_fromStoA,     "Unit": "Non dimensional"},
                    "q_fromAtoGust1" : {"Value": q_fromAtoGust1,  "Unit": "Pa"},
                    "CD_fromAtoGust1": {"Value": CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toC" : {"Value": q_fromGust1toC,  "Unit": "Pa"},
                    "CD_fromGust1toC": {"Value": CD_fromGust1toC, "Unit": "Non dimensional"},
                    "q_fromCtoGust2" : {"Value": q_fromCtoGust2,  "Unit": "Pa"},
                    "CD_fromCtoGust2": {"Value": CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "q_fromGust2toD" : {"Value": q_fromGust2toD,  "Unit": "Pa"},
                    "CD_fromGust2toD": {"Value": CD_fromGust2toD, "Unit": "Non dimensional"},
                    "q_fromDto0"     : {"Value": q_fromDto0,      "Unit": "Pa"},
                    "CD_fromDto0"    : {"Value": CD_fromDto0,     "Unit": "Non dimensional"},
                    "q_from0toSinv"  : {"Value": q_from0toSinv,   "Unit": "Pa"},
                    "CD_from0toSinv" : {"Value": CD_from0toSinv,  "Unit": "Non dimensional"},
                    "q_fromSinvtoG"  : {"Value": q_fromSinvtoG,   "Unit": "Pa"},
                    "CD_fromSinvtoG" : {"Value": CD_fromSinvtoG,  "Unit": "Non dimensional"},
                    "q_fromGtoF"     : {"Value": q_fromGtoF,      "Unit": "Pa"},
                    "CD_fromGtoF"    : {"Value": CD_fromGtoF,     "Unit": "Non dimensional"},
                    "q_fromFtoE"     : {"Value": q_fromFtoE,      "Unit": "Pa"},
                    "CD_fromFtoE"    : {"Value": CD_fromFtoE,     "Unit": "Non dimensional"},
                    "q_fromEto0"     : {"Value": q_fromEto0,      "Unit": "Pa"},
                    "CD_fromEto0"    : {"Value": CD_fromEto0,     "Unit": "Non dimensional"},

                    "CL_wb_from0toS"     : {"Value": CL_wb_from0toS,      "Unit": "Non dimensional"},
                    "CL_ht_from0toS"     : {"Value": CL_ht_from0toS,      "Unit": "Non dimensional"},
                    "CL_new_from0toS"    : {"Value": CL_new_from0toS,     "Unit": "Non dimensional"},
                    "CL_wb_fromStoA"     : {"Value": CL_wb_fromStoA,      "Unit": "Non dimensional"},
                    "CL_ht_fromStoA"     : {"Value": CL_ht_fromStoA,      "Unit": "Non dimensional"},
                    "CL_new_fromStoA"    : {"Value": CL_new_fromStoA,     "Unit": "Non dimensional"},
                    "CL_wb_fromAtoGust1" : {"Value": CL_wb_fromAtoGust1,  "Unit": "Non dimensional"},
                    "CL_ht_fromAtoGust1" : {"Value": CL_ht_fromAtoGust1,  "Unit": "Non dimensional"},
                    "CL_new_fromAtoGust1": {"Value": CL_new_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_wb_fromGust1toC" : {"Value": CL_wb_fromGust1toC,  "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toC" : {"Value": CL_ht_fromGust1toC,  "Unit": "Non dimensional"},
                    "CL_new_fromGust1toC": {"Value": CL_new_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_wb_fromCtoGust2" : {"Value": CL_wb_fromCtoGust2,  "Unit": "Non dimensional"},
                    "CL_ht_fromCtoGust2" : {"Value": CL_ht_fromCtoGust2,  "Unit": "Non dimensional"},
                    "CL_new_fromCtoGust2": {"Value": CL_new_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_wb_fromGust2toD" : {"Value": CL_wb_fromGust2toD,  "Unit": "Non dimensional"},
                    "CL_ht_fromGust2toD" : {"Value": CL_ht_fromGust2toD,  "Unit": "Non dimensional"},
                    "CL_new_fromGust2toD": {"Value": CL_new_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_wb_fromDto0"     : {"Value": CL_wb_fromDto0,      "Unit": "Non dimensional"},
                    "CL_ht_fromDto0"     : {"Value": CL_ht_fromDto0,      "Unit": "Non dimensional"},
                    "CL_new_fromDto0"    : {"Value": CL_new_fromDto0,     "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv"  : {"Value": CL_wb_from0toSinv,   "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv"  : {"Value": CL_ht_from0toSinv,   "Unit": "Non dimensional"},
                    "CL_new_from0toSinv" : {"Value": CL_new_from0toSinv,  "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG"  : {"Value": CL_wb_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG"  : {"Value": CL_ht_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG" : {"Value": CL_ht_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CL_wb_fromGtoF"     : {"Value": CL_wb_fromGtoF,      "Unit": "Non dimensional"},
                    "CL_ht_fromGtoF"     : {"Value": CL_ht_fromGtoF,      "Unit": "Non dimensional"},
                    "CL_new_fromGtoF"    : {"Value": CL_new_fromGtoF,     "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE"     : {"Value": CL_wb_fromFtoE,      "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE"     : {"Value": CL_ht_fromFtoE,      "Unit": "Non dimensional"},
                    "CL_new_fromFtoE"    : {"Value": CL_new_fromFtoE,     "Unit": "Non dimensional"},
                    "CL_wb_fromEto0"     : {"Value": CL_wb_fromEto0,      "Unit": "Non dimensional"},
                    "CL_ht_fromEto0"     : {"Value": CL_ht_fromEto0,      "Unit": "Non dimensional"},
                    "CL_new_fromEto0"    : {"Value": CL_new_fromEto0,     "Unit": "Non dimensional"},

                    "CM_due_to_CL_from0toS"        : {"Value": CM_due_to_CL_from0toS,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS"        : {"Value": CM_due_to_CD_from0toS,      "Unit": "Non dimensional"},
                    "CM_CG_from0toS"               : {"Value": CM_CG_from0toS,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA"        : {"Value": CM_due_to_CL_fromStoA,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA"        : {"Value": CM_due_to_CD_fromStoA,      "Unit": "Non dimensional"},
                    "CM_CG_fromStoA"               : {"Value" : CM_CG_fromStoA,            "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoGust1"    : {"Value": CM_due_to_CL_fromAtoGust1,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoGust1"    : {"Value": CM_due_to_CD_fromAtoGust1,  "Unit": "Non dimensional"},
                    "CM_CG_fromAtoGust1"           : {"Value" : CM_CG_fromAtoGust1,        "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toC"    : {"Value": CM_due_to_CL_fromGust1toC,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toC"    : {"Value": CM_due_to_CD_fromGust1toC,  "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toC"           : {"Value" : CM_CG_fromGust1toC,        "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoGust2"    : {"Value": CM_due_to_CL_fromCtoGust2,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoGust2"    : {"Value": CM_due_to_CD_fromCtoGust2,  "Unit": "Non dimensional"},
                    "CM_CG_fromCtoGust2"           : {"Value" : CM_CG_fromCtoGust2,        "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust2toD"    : {"Value": CM_due_to_CL_fromGust2toD,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust2toD"    : {"Value": CM_due_to_CD_fromGust2toD,  "Unit": "Non dimensional"},
                    "CM_CG_fromGust2toD"           : {"Value" : CM_CG_fromGust2toD,        "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0"        : {"Value": CM_due_to_CL_fromDto0,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0"        : {"Value": CM_due_to_CD_fromDto0,      "Unit": "Non dimensional"},
                    "CM_CG_fromDto0"               : {"Value" : CM_CG_fromDto0,            "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv"     : {"Value": CM_due_to_CL_from0toSinv,   "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv"     : {"Value": CM_due_to_CD_from0toSinv,   "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv"            : {"Value" : CM_CG_from0toSinv,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG"     : {"Value": CM_due_to_CL_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG"     : {"Value": CM_due_to_CD_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG"            : {"Value" : CM_due_to_CD_fromSinvtoG,  "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoF"        : {"Value": CM_due_to_CL_fromGtoF,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoF"        : {"Value": CM_due_to_CD_fromGtoF,      "Unit": "Non dimensional"},
                    "CM_CG_fromGtoF"               : {"Value": CM_CG_fromGtoF,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE"        : {"Value": CM_due_to_CL_fromFtoE,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE"        : {"Value": CM_due_to_CD_fromFtoE,      "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE"               : {"Value": CM_CG_fromFtoE,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0"        : {"Value": CM_due_to_CL_fromEto0,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0"        : {"Value": CM_due_to_CD_fromEto0,      "Unit": "Non dimensional"},
                    "CM_CG_fromEto0"               : {"Value": CM_CG_fromEto0,             "Unit": "Non dimensional"},

                    "L_wb_from0toS"     : {"Value": L_wb_from0toS,      "Unit": "Non dimensional"},
                    "L_ht_from0toS"     : {"Value": L_ht_from0toS,      "Unit": "Non dimensional"},
                    "L_new_from0toS"    : {"Value": L_new_from0toS,     "Unit": "Non dimensional"},
                    "L_wb_fromStoA"     : {"Value": L_wb_fromStoA,      "Unit": "Non dimensional"},
                    "L_ht_fromStoA"     : {"Value": L_ht_fromStoA,      "Unit": "Non dimensional"},
                    "L_new_fromStoA"    : {"Value": L_new_fromStoA,     "Unit": "Non dimensional"},
                    "L_wb_fromAtoGust1" : {"Value": L_wb_fromAtoGust1,  "Unit": "Non dimensional"},
                    "L_ht_fromAtoGust1" : {"Value": L_ht_fromAtoGust1,  "Unit": "Non dimensional"},
                    "L_new_fromAtoGust1": {"Value": L_new_fromAtoGust1, "Unit": "Non dimensional"},
                    "L_wb_fromGust1toC" : {"Value": L_wb_fromGust1toC,  "Unit": "Non dimensional"},
                    "L_ht_fromGust1toC" : {"Value": L_ht_fromGust1toC,  "Unit": "Non dimensional"},
                    "L_new_fromGust1toC": {"Value": L_new_fromGust1toC, "Unit": "Non dimensional"},
                    "L_wb_fromCtoGust2" : {"Value": L_wb_fromCtoGust2,  "Unit": "Non dimensional"},
                    "L_ht_fromCtoGust2" : {"Value": L_ht_fromCtoGust2,  "Unit": "Non dimensional"},
                    "L_new_fromCtoGust2": {"Value": L_new_fromCtoGust2, "Unit": "Non dimensional"},
                    "L_wb_fromGust2toD" : {"Value": L_wb_fromGust2toD,  "Unit": "Non dimensional"},
                    "L_ht_fromGust2toD" : {"Value": L_ht_fromGust2toD,  "Unit": "Non dimensional"},
                    "L_new_fromGust2toD": {"Value": L_new_fromGust2toD, "Unit": "Non dimensional"},
                    "L_wb_fromDto0"     : {"Value": L_wb_fromDto0,      "Unit": "Non dimensional"},
                    "L_ht_fromDto0"     : {"Value": L_ht_fromDto0,      "Unit": "Non dimensional"},
                    "L_new_fromDto0"    : {"Value": L_new_fromDto0,     "Unit": "Non dimensional"},
                    "L_wb_from0toSinv"  : {"Value": L_wb_from0toSinv,   "Unit": "Non dimensional"},
                    "L_ht_from0toSinv"  : {"Value": L_ht_from0toSinv,   "Unit": "Non dimensional"},
                    "L_new_from0toSinv" : {"Value": L_new_from0toSinv,  "Unit": "Non dimensional"},
                    "L_wb_fromSinvtoG"  : {"Value": L_wb_fromSinvtoG,   "Unit": "Non dimensional"},
                    "L_ht_fromSinvtoG"  : {"Value": L_ht_fromSinvtoG,   "Unit": "Non dimensional"},
                    "L_new_fromSinvtoG" : {"Value": L_ht_fromSinvtoG,   "Unit": "Non dimensional"},
                    "L_wb_fromGtoF"     : {"Value": L_wb_fromGtoF,      "Unit": "Non dimensional"},
                    "L_ht_fromGtoF"     : {"Value": L_ht_fromGtoF,      "Unit": "Non dimensional"},
                    "L_new_fromGtoF"    : {"Value": L_new_fromGtoF,     "Unit": "Non dimensional"},
                    "L_wb_fromFtoE"     : {"Value": L_wb_fromFtoE,      "Unit": "Non dimensional"},
                    "L_ht_fromFtoE"     : {"Value": L_ht_fromFtoE,      "Unit": "Non dimensional"},
                    "L_new_fromFtoE"    : {"Value": L_new_fromFtoE,     "Unit": "Non dimensional"},
                    "L_wb_fromEto0"     : {"Value": L_wb_fromEto0,      "Unit": "Non dimensional"},
                    "L_ht_fromEto0"     : {"Value": L_ht_fromEto0,      "Unit": "Non dimensional"},
                    "L_new_fromEto0"    : {"Value": L_new_fromEto0,     "Unit": "Non dimensional"},
 
                    "alfa_from0toS"        : {"Value": alfa_from0toS,         "Unit": "deg"},
                    "alfa_new_from0toS"    : {"Value": alfa_new_from0toS,     "Unit": "deg"},
                    "alfa_fromStoA"        : {"Value": alfa_fromStoA,         "Unit": "deg"},
                    "alfa_new_fromStoA"    : {"Value": alfa_new_fromStoA,     "Unit": "deg"},
                    "alfa_fromAtoGust1"    : {"Value": alfa_fromAtoGust1,     "Unit": "deg"},
                    "alfa_new_fromAtoGust1": {"Value": alfa_new_fromAtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toC"    : {"Value": alfa_fromGust1toC,     "Unit": "deg"},
                    "alfa_new_fromGust1toC": {"Value": alfa_new_fromGust1toC, "Unit": "deg"},
                    "alfa_fromCtoGust2"    : {"Value": alfa_fromCtoGust2,     "Unit": "deg"},
                    "alfa_new_fromCtoGust2": {"Value": alfa_new_fromCtoGust2, "Unit": "deg"},
                    "alfa_fromGust2toD"    : {"Value": alfa_fromGust2toD,     "Unit": "deg"},
                    "alfa_new_fromGust2toD": {"Value": alfa_new_fromGust2toD, "Unit": "deg"},
                    "alfa_fromDto0"        : {"Value": alfa_fromDto0,         "Unit": "deg"},
                    "alfa_new_fromDto0"    : {"Value": alfa_new_fromDto0,     "Unit": "deg"},
                    "alfa_from0toSinv"     : {"Value": alfa_from0toSinv,      "Unit": "deg"},
                    "alfa_new_from0toSinv" : {"Value": alfa_new_from0toSinv,  "Unit": "deg"},
                    "alfa_fromSinvtoG"     : {"Value": alfa_fromSinvtoG,      "Unit": "deg"},
                    "alfa_new_fromSinvtoG" : {"Value": alfa_new_fromSinvtoG,  "Unit": "deg"},
                    "alfa_fromGtoF"        : {"Value": alfa_fromGtoF,         "Unit": "deg"},
                    "alfa_new_fromGtoF"    : {"Value": alfa_new_fromGtoF,     "Unit": "deg"},
                    "alfa_fromFtoE"        : {"Value": alfa_fromFtoE,         "Unit": "deg"},
                    "alfa_new_fromFtoE"    : {"Value": alfa_new_fromFtoE,     "Unit": "deg"},
                    "alfa_fromEto0"        : {"Value": alfa_fromEto0,         "Unit": "deg"},
                    "alfa_new_fromEto0"    : {"Value": alfa_new_fromEto0,     "Unit": "deg"}
                    }
                }

        elif (pos_case_flag == "Case1") and (neg_case_flag == "Case3_inverted"):
            
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS": {"Value": V_fe_from0toS, "Unit": "m/s"},
                    "n_fe_from0toS": {"Value": n_fe_from0toS, "Unit": "g"},
                    "V_fe_fromStoA": {"Value": V_fe_fromStoA, "Unit": "m/s"},
                    "n_fe_fromStoA": {"Value": n_fe_fromStoA, "Unit": "g"},
                    "V_fe_fromAtoGust1": {"Value": V_fe_fromAtoGust1, "Unit": "m/s"},
                    "n_fe_fromAtoGust1": {"Value": n_fe_fromAtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toC": {"Value": V_fe_fromGust1toC, "Unit": "m/s"},
                    "n_fe_fromGust1toC": {"Value": n_fe_fromGust1toC, "Unit": "g"},
                    "V_fe_fromCtoGust2": {"Value": V_fe_fromCtoGust2, "Unit": "m/s"},
                    "n_fe_fromCtoGust2": {"Value": n_fe_fromCtoGust2, "Unit": "g"},
                    "V_fe_fromGust2toD": {"Value": V_fe_fromGust2toD, "Unit": "m/s"},
                    "n_fe_fromGust2toD": {"Value": n_fe_fromGust2toD, "Unit": "g"},
                    "V_fe_fromDto0":     {"Value": V_fe_fromDto0, "Unit": "m/s"},
                    "n_fe_fromDto0":     {"Value": n_fe_fromDto0, "Unit": "g"},
                    "V_fe_from0toSinv":  {"Value": V_fe_from0toSinv, "Unit": "m/s"},
                    "n_fe_from0toSinv":  {"Value": n_fe_from0toSinv, "Unit": "g"},
                    "V_fe_fromSinvtoG":  {"Value": V_fe_fromSinvtoG, "Unit": "m/s"},
                    "n_fe_fromSinvtoG":  {"Value": n_fe_fromSinvtoG, "Unit": "g"},
                    "V_fe_fromGtoF":     {"Value": V_fe_fromGtoF, "Unit": "m/s"},
                    "n_fe_fromGtoF":     {"Value": n_fe_fromGtoF, "Unit": "g"},
                    "V_fe_fromFtoE"    : {"Value": V_fe_fromFtoE, "Unit": "m/s"},
                    "n_fe_fromFtoE"    : {"Value": n_fe_fromFtoE, "Unit": "g"},
                    "V_fe_fromEto0"    : {"Value": V_fe_fromEto0, "Unit": "m/s"},
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"},

                    "q_from0toS"     : {"Value": q_from0toS,      "Unit": "Pa"},
                    "CD_from0toS"    : {"Value": CD_from0toS,     "Unit": "Non dimensional"},
                    "q_fromStoA"     : {"Value": q_fromStoA,      "Unit": "Pa"},
                    "CD_fromStoA"    : {"Value": CD_fromStoA,     "Unit": "Non dimensional"},
                    "q_fromAtoGust1" : {"Value": q_fromAtoGust1,  "Unit": "Pa"},
                    "CD_fromAtoGust1": {"Value": CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toC" : {"Value": q_fromGust1toC,  "Unit": "Pa"},
                    "CD_fromGust1toC": {"Value": CD_fromGust1toC, "Unit": "Non dimensional"},
                    "q_fromCtoGust2" : {"Value": q_fromCtoGust2,  "Unit": "Pa"},
                    "CD_fromCtoGust2": {"Value": CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "q_fromGust2toD" : {"Value": q_fromGust2toD,  "Unit": "Pa"},
                    "CD_fromGust2toD": {"Value": CD_fromGust2toD, "Unit": "Non dimensional"},
                    "q_fromDto0"     : {"Value": q_fromDto0,      "Unit": "Pa"},
                    "CD_fromDto0"    : {"Value": CD_fromDto0,     "Unit": "Non dimensional"},
                    "q_from0toSinv"  : {"Value": q_from0toSinv,   "Unit": "Pa"},
                    "CD_from0toSinv" : {"Value": CD_from0toSinv,  "Unit": "Non dimensional"},
                    "q_fromSinvtoG"  : {"Value": q_fromSinvtoG,   "Unit": "Pa"},
                    "CD_fromSinvtoG" : {"Value": CD_fromSinvtoG,  "Unit": "Non dimensional"},
                    "q_fromGtoF"     : {"Value": q_fromGtoF,      "Unit": "Pa"},
                    "CD_fromGtoF"    : {"Value": CD_fromGtoF,     "Unit": "Non dimensional"},
                    "q_fromFtoE"     : {"Value": q_fromFtoE,      "Unit": "Pa"},
                    "CD_fromFtoE"    : {"Value": CD_fromFtoE,     "Unit": "Non dimensional"},
                    "q_fromEto0"     : {"Value": q_fromEto0,      "Unit": "Pa"},
                    "CD_fromEto0"    : {"Value": CD_fromEto0,     "Unit": "Non dimensional"},                  

                    "alfa_from0toS"        : {"Value": alfa_from0toS,         "Unit": "deg"},
                    "alfa_new_from0toS"    : {"Value": alfa_new_from0toS,     "Unit": "deg"},
                    "alfa_fromStoA"        : {"Value": alfa_fromStoA,         "Unit": "deg"},
                    "alfa_new_fromStoA"    : {"Value": alfa_new_fromStoA,     "Unit": "deg"},
                    "alfa_fromAtoGust1"    : {"Value": alfa_fromAtoGust1,     "Unit": "deg"},
                    "alfa_new_fromAtoGust1": {"Value": alfa_new_fromAtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toC"    : {"Value": alfa_fromGust1toC,     "Unit": "deg"},
                    "alfa_new_fromGust1toC": {"Value": alfa_new_fromGust1toC, "Unit": "deg"},
                    "alfa_fromCtoGust2"    : {"Value": alfa_fromCtoGust2,     "Unit": "deg"},
                    "alfa_new_fromCtoGust2": {"Value": alfa_new_fromCtoGust2, "Unit": "deg"},
                    "alfa_fromGust2toD"    : {"Value": alfa_fromGust2toD,     "Unit": "deg"},
                    "alfa_new_fromGust2toD": {"Value": alfa_new_fromGust2toD, "Unit": "deg"},
                    "alfa_fromDto0"        : {"Value": alfa_fromDto0,         "Unit": "deg"},
                    "alfa_new_fromDto0"    : {"Value": alfa_new_fromDto0,     "Unit": "deg"},
                    "alfa_from0toSinv"     : {"Value": alfa_from0toSinv,      "Unit": "deg"},
                    "alfa_new_from0toSinv" : {"Value": alfa_new_from0toSinv,  "Unit": "deg"},
                    "alfa_fromSinvtoG"     : {"Value": alfa_fromSinvtoG,      "Unit": "deg"},
                    "alfa_new_fromSinvtoG" : {"Value": alfa_new_fromSinvtoG,  "Unit": "deg"},
                    "alfa_fromGtoF"        : {"Value": alfa_fromGtoF,         "Unit": "deg"},
                    "alfa_new_fromGtoF"    : {"Value": alfa_new_fromGtoF,     "Unit": "deg"},
                    "alfa_fromFtoE"        : {"Value": alfa_fromFtoE,         "Unit": "deg"},
                    "alfa_new_fromFtoE"    : {"Value": alfa_new_fromFtoE,     "Unit": "deg"},
                    "alfa_fromEto0"        : {"Value": alfa_fromEto0,         "Unit": "deg"},
                    "alfa_new_fromEto0"    : {"Value": alfa_new_fromEto0,     "Unit": "deg"},      
    
                    "CL_wb_from0toS"      : {"Value": CL_wb_from0toS,     "Unit": "Non dimensional"},
                    "CL_ht_from0toS"      : {"Value": CL_ht_from0toS,     "Unit": "Non dimensional"},
                    "CL_new_from0toS"     : {"Value": CL_new_from0toS,    "Unit": "Non dimensional"},
                    "CL_wb_fromStoA"      : {"Value": CL_wb_fromStoA,     "Unit": "Non dimensional"},
                    "CL_ht_fromStoA"      : {"Value": CL_ht_fromStoA,     "Unit": "Non dimensional"},
                    "CL_new_fromStoA"     : {"Value": CL_new_fromStoA,    "Unit": "Non dimensional"},
                    "CL_wb_fromAtoGust1"  : {"Value": CL_wb_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_ht_fromAtoGust1"  : {"Value": CL_ht_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_new_fromAtoGust1" : {"Value": CL_new_fromAtoGust1,"Unit": "Non dimensional"},
                    "CL_wb_fromGust1toC"  : {"Value": CL_wb_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toC"  : {"Value": CL_ht_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_new_fromGust1toC" : {"Value": CL_new_fromGust1toC,"Unit": "Non dimensional"},
                    "CL_wb_fromCtoGust2"  : {"Value": CL_wb_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_ht_fromCtoGust2"  : {"Value": CL_ht_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_new_fromCtoGust2" : {"Value": CL_new_fromCtoGust2,"Unit": "Non dimensional"},
                    "CL_wb_fromGust2toD"  : {"Value": CL_wb_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_ht_fromGust2toD"  : {"Value": CL_ht_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_new_fromGust2toD" : {"Value": CL_new_fromGust2toD,"Unit": "Non dimensional"},
                    "CL_wb_fromDto0"      : {"Value": CL_wb_fromDto0,     "Unit": "Non dimensional"},
                    "CL_ht_fromDto0"      : {"Value": CL_ht_fromDto0,     "Unit": "Non dimensional"},
                    "CL_new_fromDto0"     : {"Value": CL_new_fromDto0,    "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv"   : {"Value": CL_wb_from0toSinv,  "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv"   : {"Value": CL_ht_from0toSinv,  "Unit": "Non dimensional"},
                    "CL_new_from0toSinv"  : {"Value": CL_new_from0toSinv, "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG"   : {"Value": CL_wb_fromSinvtoG,  "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG"   : {"Value": CL_ht_fromSinvtoG,  "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG"  : {"Value": CL_new_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_wb_fromGtoF"      : {"Value": CL_wb_fromGtoF,     "Unit": "Non dimensional"},
                    "CL_ht_fromGtoF"      : {"Value": CL_ht_fromGtoF,     "Unit": "Non dimensional"},
                    "CL_new_fromGtoF"     : {"Value": CL_new_fromGtoF,    "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE"      : {"Value": CL_wb_fromFtoE,     "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE"      : {"Value": CL_ht_fromFtoE,     "Unit": "Non dimensional"},
                    "CL_new_fromFtoE"     : {"Value": CL_new_fromFtoE,    "Unit": "Non dimensional"},
                    "CL_wb_fromEto0"      : {"Value": CL_wb_fromEto0,     "Unit": "Non dimensional"},
                    "CL_ht_fromEto0"      : {"Value": CL_ht_fromEto0,     "Unit": "Non dimensional"},
                    "CL_new_fromEto0"     : {"Value": CL_new_fromEto0,    "Unit": "Non dimensional"},
    
                    "L_wb_from0toS"      : {"Value": L_wb_from0toS,     "Unit": "daN"},
                    "L_ht_from0toS"      : {"Value": L_ht_from0toS,     "Unit": "daN"},
                    "L_new_from0toS"     : {"Value": L_new_from0toS,    "Unit": "daN"},
                    "L_wb_fromStoA"      : {"Value": L_wb_fromStoA,     "Unit": "daN"},
                    "L_ht_fromStoA"      : {"Value": L_ht_fromStoA,     "Unit": "daN"},
                    "L_new_fromStoA"     : {"Value": L_new_fromStoA,    "Unit": "daN"},
                    "L_wb_fromAtoGust1"  : {"Value": L_wb_fromAtoGust1, "Unit": "daN"},
                    "L_ht_fromAtoGust1"  : {"Value": L_ht_fromAtoGust1, "Unit": "daN"},
                    "L_new_fromAtoGust1" : {"Value": L_new_fromAtoGust1,"Unit": "daN"},
                    "L_wb_fromGust1toC"  : {"Value": L_wb_fromGust1toC, "Unit": "daN"},
                    "L_ht_fromGust1toC"  : {"Value": L_ht_fromGust1toC, "Unit": "daN"},
                    "L_new_fromGust1toC" : {"Value": L_new_fromGust1toC,"Unit": "daN"},
                    "L_wb_fromCtoGust2"  : {"Value": L_wb_fromCtoGust2, "Unit": "daN"},
                    "L_ht_fromCtoGust2"  : {"Value": L_ht_fromCtoGust2, "Unit": "daN"},
                    "L_new_fromCtoGust2" : {"Value": L_new_fromCtoGust2,"Unit": "daN"},
                    "L_wb_fromGust2toD"  : {"Value": L_wb_fromGust2toD, "Unit": "daN"},
                    "L_ht_fromGust2toD"  : {"Value": L_ht_fromGust2toD, "Unit": "daN"},
                    "L_new_fromGust2toD" : {"Value": L_new_fromGust2toD,"Unit": "daN"},
                    "L_wb_fromDto0"      : {"Value": L_wb_fromDto0,     "Unit": "daN"},
                    "L_ht_fromDto0"      : {"Value": L_ht_fromDto0,     "Unit": "daN"},
                    "L_new_fromDto0"     : {"Value": L_new_fromDto0,    "Unit": "daN"},
                    "L_wb_from0toSinv"   : {"Value": L_wb_from0toSinv,  "Unit": "daN"},
                    "L_ht_from0toSinv"   : {"Value": L_ht_from0toSinv,  "Unit": "daN"},
                    "L_new_from0toSinv"  : {"Value": L_new_from0toSinv, "Unit": "daN"},
                    "L_wb_fromSinvtoG"   : {"Value": L_wb_fromSinvtoG,  "Unit": "daN"},
                    "L_ht_fromSinvtoG"   : {"Value": L_ht_fromSinvtoG,  "Unit": "daN"},
                    "L_new_fromSinvtoG"  : {"Value": L_new_fromSinvtoG, "Unit": "daN"},
                    "L_wb_fromGtoF"      : {"Value": L_wb_fromGtoF,     "Unit": "daN"},
                    "L_ht_fromGtoF"      : {"Value": L_ht_fromGtoF,     "Unit": "daN"},
                    "L_new_fromGtoF"     : {"Value": L_new_fromGtoF,    "Unit": "daN"},
                    "L_wb_fromFtoE"      : {"Value": L_wb_fromFtoE,     "Unit": "daN"},
                    "L_ht_fromFtoE"      : {"Value": L_ht_fromFtoE,     "Unit": "daN"},
                    "L_new_fromFtoE"     : {"Value": L_new_fromFtoE,    "Unit": "daN"},
                    "L_wb_fromEto0"      : {"Value": L_wb_fromEto0,     "Unit": "daN"},
                    "L_ht_fromEto0"      : {"Value": L_ht_fromEto0,     "Unit": "daN"},
                    "L_new_fromEto0"     : {"Value": L_new_fromEto0,    "Unit": "daN"},
    
                    "CM_due_to_CL_from0toS"      : {"Value": CM_due_to_CL_from0toS,     "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS"      : {"Value": CM_due_to_CD_from0toS,     "Unit": "Non dimensional"},
                    "CM_CG_from0toS"             : {"Value": CM_CG_from0toS,            "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA"      : {"Value": CM_due_to_CL_fromStoA,     "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA"      : {"Value": CM_due_to_CD_fromStoA,     "Unit": "Non dimensional"},
                    "CM_CG_fromStoA"             : {"Value": CM_CG_fromStoA,            "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoGust1"  : {"Value": CM_due_to_CL_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoGust1"  : {"Value": CM_due_to_CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_CG_fromAtoGust1"         : {"Value": CM_CG_fromAtoGust1,        "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toC"  : {"Value": CM_due_to_CL_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toC"  : {"Value": CM_due_to_CD_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toC"         : {"Value": CM_CG_fromGust1toC,        "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoGust2"  : {"Value": CM_due_to_CL_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoGust2"  : {"Value": CM_due_to_CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_CG_fromCtoGust2"         : {"Value": CM_CG_fromCtoGust2,        "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust2toD"  : {"Value": CM_due_to_CL_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust2toD"  : {"Value": CM_due_to_CD_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_CG_fromGust2toD"         : {"Value": CM_CG_fromGust2toD,        "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0"      : {"Value": CM_due_to_CL_fromDto0,     "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0"      : {"Value": CM_due_to_CD_fromDto0,     "Unit": "Non dimensional"},
                    "CM_CG_fromDto0"             : {"Value": CM_CG_fromDto0,            "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv"   : {"Value": CM_due_to_CL_from0toSinv,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv"   : {"Value": CM_due_to_CD_from0toSinv,  "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv"          : {"Value": CM_CG_from0toSinv,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG"   : {"Value": CM_due_to_CL_fromSinvtoG,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG"   : {"Value": CM_due_to_CD_fromSinvtoG,  "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG"          : {"Value": CM_CG_fromSinvtoG,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoF"      : {"Value": CM_due_to_CL_fromGtoF,     "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoF"      : {"Value": CM_due_to_CD_fromGtoF,     "Unit": "Non dimensional"},
                    "CM_CG_fromGtoF"             : {"Value": CM_CG_fromGtoF,            "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE"      : {"Value": CM_due_to_CL_fromFtoE,     "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE"      : {"Value": CM_due_to_CD_fromFtoE,     "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE"             : {"Value": CM_CG_fromFtoE,            "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0"      : {"Value": CM_due_to_CL_fromEto0,     "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0"      : {"Value": CM_due_to_CD_fromEto0,     "Unit": "Non dimensional"},
                    "CM_CG_fromEto0"             : {"Value": CM_CG_fromEto0,            "Unit": "Non dimensional"}
                    }
                }

        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case1_inverted"):
            
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS":     {"Value": V_fe_from0toS, "Unit": "m/s"},
                    "n_fe_from0toS":     {"Value": n_fe_from0toS, "Unit": "g"},
                    "V_fe_fromStoA":     {"Value": V_fe_fromStoA, "Unit": "m/s"},
                    "n_fe_fromStoA":     {"Value": n_fe_fromStoA, "Unit": "g"},
                    "V_fe_fromAtoC":     {"Value": V_fe_fromAtoC, "Unit": "m/s"},
                    "n_fe_fromAtoC":     {"Value": n_fe_fromAtoC, "Unit": "g"},
                    "V_fe_fromCtoD":     {"Value": V_fe_fromCtoD, "Unit": "m/s"},
                    "n_fe_fromCtoD":     {"Value": n_fe_fromCtoD, "Unit": "g"},
                    "V_fe_fromDto0":     {"Value": V_fe_fromDto0, "Unit": "m/s"},
                    "n_fe_fromDto0":     {"Value": n_fe_fromDto0, "Unit": "g"},
                    "V_fe_from0toSinv":  {"Value": V_fe_from0toSinv, "Unit": "m/s"},
                    "n_fe_from0toSinv":  {"Value": n_fe_from0toSinv, "Unit": "g"},
                    "V_fe_fromSinvtoG":  {"Value": V_fe_fromSinvtoG, "Unit": "m/s"},
                    "n_fe_fromSinvtoG":  {"Value": n_fe_fromSinvtoG, "Unit": "g"},
                    "V_fe_fromGtoGust1": {"Value": V_fe_fromGtoGust1, "Unit": "m/s"},
                    "n_fe_fromGtoGust1": {"Value": n_fe_fromGtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toF": {"Value": V_fe_fromGust1toF, "Unit": "m/s"},
                    "n_fe_fromGust1toF": {"Value": n_fe_fromGust1toF, "Unit": "g"},
                    "V_fe_fromFtoE"    : {"Value": V_fe_fromFtoE, "Unit": "m/s"},
                    "n_fe_fromFtoE"    : {"Value": n_fe_fromFtoE, "Unit": "g"},
                    "V_fe_fromEto0"    : {"Value": V_fe_fromEto0, "Unit": "m/s"},
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"},
 
                    "q_from0toS"     : {"Value": q_from0toS,      "Unit": "Pa"},
                    "CD_from0toS"    : {"Value": CD_from0toS,     "Unit": "Non dimensional"},
                    "q_fromStoA"     : {"Value": q_fromStoA,      "Unit": "Pa"},
                    "CD_fromStoA"    : {"Value": CD_fromStoA,     "Unit": "Non dimensional"},
                    "q_fromAtoC"     : {"Value": q_fromAtoC,      "Unit": "Pa"},
                    "CD_fromAtoC"    : {"Value": CD_fromAtoC,     "Unit": "Non dimensional"},
                    "q_fromCtoD"     : {"Value": q_fromCtoD,      "Unit": "Pa"},
                    "CD_fromCtoD"    : {"Value": CD_fromCtoD,     "Unit": "Non dimensional"},
                    "q_fromDto0"     : {"Value": q_fromDto0,      "Unit": "Pa"},
                    "CD_fromDto0"    : {"Value": CD_fromDto0,     "Unit": "Non dimensional"},
                    "q_from0toSinv"  : {"Value": q_from0toSinv,   "Unit": "Pa"},
                    "CD_from0toSinv" : {"Value": CD_from0toSinv,  "Unit": "Non dimensional"},
                    "q_fromSinvtoG"  : {"Value": q_fromSinvtoG,   "Unit": "Pa"},
                    "CD_fromSinvtoG" : {"Value": CD_fromSinvtoG,  "Unit": "Non dimensional"},
                    "q_fromGtoGust1" : {"Value": q_fromGtoGust1,  "Unit": "Pa"},
                    "CD_fromGtoGust1": {"Value": CD_fromGtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toF" : {"Value": q_fromGust1toF,  "Unit": "Pa"},
                    "CD_fromGust1toF": {"Value": CD_fromGust1toF, "Unit": "Non dimensional"},
                    "q_fromFtoE"     : {"Value": q_fromFtoE,      "Unit": "Pa"},
                    "CD_fromFtoE"    : {"Value": CD_fromFtoE,     "Unit": "Non dimensional"},
                    "q_fromEto0"     : {"Value": q_fromEto0,      "Unit": "Pa"},
                    "CD_fromEto0"    : {"Value": CD_fromEto0,     "Unit": "Non dimensional"},

                    "alfa_from0toS"        : {"Value": alfa_from0toS,         "Unit": "deg"},
                    "alfa_new_from0toS"    : {"Value": alfa_new_from0toS,     "Unit": "deg"},
                    "alfa_fromStoA"        : {"Value": alfa_fromStoA,         "Unit": "deg"},
                    "alfa_new_fromStoA"    : {"Value": alfa_new_fromStoA,     "Unit": "deg"},
                    "alfa_fromAtoC"        : {"Value": alfa_fromAtoC,         "Unit": "deg"},
                    "alfa_new_fromAtoC"    : {"Value": alfa_new_fromAtoC,     "Unit": "deg"},
                    "alfa_fromCtoD"        : {"Value": alfa_fromCtoD,         "Unit": "deg"},
                    "alfa_new_fromCtoD"    : {"Value": alfa_new_fromCtoD,     "Unit": "deg"},
                    "alfa_fromDto0"        : {"Value": alfa_fromDto0,         "Unit": "deg"},
                    "alfa_new_fromDto0"    : {"Value": alfa_new_fromDto0,     "Unit": "deg"},
                    "alfa_from0toSinv"     : {"Value": alfa_from0toSinv,      "Unit": "deg"},
                    "alfa_new_from0toSinv" : {"Value": alfa_new_from0toSinv,  "Unit": "deg"},
                    "alfa_fromSinvtoG"     : {"Value": alfa_fromSinvtoG,      "Unit": "deg"},
                    "alfa_new_fromSinvtoG" : {"Value": alfa_new_fromSinvtoG,  "Unit": "deg"},
                    "alfa_fromGtoGust1"    : {"Value": alfa_fromGtoGust1,     "Unit": "deg"},
                    "alfa_new_fromGtoGust1": {"Value": alfa_new_fromGtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toF"    : {"Value": alfa_fromGust1toF,     "Unit": "deg"},
                    "alfa_new_fromGust1toF": {"Value": alfa_new_fromGust1toF, "Unit": "deg"},
                    "alfa_fromFtoE"        : {"Value": alfa_fromFtoE,         "Unit": "deg"},
                    "alfa_new_fromFtoE"    : {"Value": alfa_new_fromFtoE,     "Unit": "deg"},
                    "alfa_fromEto0"        : {"Value": alfa_fromEto0,         "Unit": "deg"},
                    "alfa_new_fromEto0"    : {"Value": alfa_new_fromEto0,     "Unit": "deg"},

                    "CL_wb_from0toS"     : {"Value": CL_wb_from0toS,      "Unit": "Non dimensional"},
                    "CL_ht_from0toS"     : {"Value": CL_ht_from0toS,      "Unit": "Non dimensional"},
                    "CL_new_from0toS"    : {"Value": CL_new_from0toS,     "Unit": "Non dimensional"},
                    "CL_wb_fromStoA"     : {"Value": CL_wb_fromStoA,      "Unit": "Non dimensional"},
                    "CL_ht_fromStoA"     : {"Value": CL_ht_fromStoA,      "Unit": "Non dimensional"},
                    "CL_new_fromStoA"    : {"Value": CL_new_fromStoA,     "Unit": "Non dimensional"},
                    "CL_wb_fromAtoC"     : {"Value": CL_wb_fromAtoC,      "Unit": "Non dimensional"},
                    "CL_ht_fromAtoC"     : {"Value": CL_ht_fromAtoC,      "Unit": "Non dimensional"},
                    "CL_new_fromAtoC"    : {"Value": CL_new_fromAtoC,     "Unit": "Non dimensional"},
                    "CL_wb_fromCtoD"     : {"Value": CL_wb_fromCtoD,      "Unit": "Non dimensional"},
                    "CL_ht_fromCtoD"     : {"Value": CL_ht_fromCtoD,      "Unit": "Non dimensional"},
                    "CL_new_fromCtoD"    : {"Value": CL_new_fromCtoD,     "Unit": "Non dimensional"},
                    "CL_wb_fromDto0"     : {"Value": CL_wb_fromDto0,      "Unit": "Non dimensional"},
                    "CL_ht_fromDto0"     : {"Value": CL_ht_fromDto0,      "Unit": "Non dimensional"},
                    "CL_new_fromDto0"    : {"Value": CL_new_fromDto0,     "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv"  : {"Value": CL_wb_from0toSinv,   "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv"  : {"Value": CL_ht_from0toSinv,   "Unit": "Non dimensional"},
                    "CL_new_from0toSinv" : {"Value": CL_new_from0toSinv,  "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG"  : {"Value": CL_wb_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG"  : {"Value": CL_ht_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG" : {"Value": CL_new_fromSinvtoG,  "Unit": "Non dimensional"},
                    "CL_wb_fromGtoGust1" : {"Value": CL_wb_fromGtoGust1,  "Unit": "Non dimensional"},
                    "CL_ht_fromGtoGust1" : {"Value": CL_ht_fromGtoGust1,  "Unit": "Non dimensional"},
                    "CL_new_fromGtoGust1": {"Value": CL_new_fromGtoGust1, "Unit": "Non dimensional"},
                    "CL_wb_fromGust1toF" : {"Value": CL_wb_fromGust1toF,  "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toF" : {"Value": CL_ht_fromGust1toF,  "Unit": "Non dimensional"},
                    "CL_new_fromGust1toF": {"Value": CL_new_fromGust1toF, "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE"     : {"Value": CL_wb_fromFtoE,      "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE"     : {"Value": CL_ht_fromFtoE,      "Unit": "Non dimensional"},
                    "CL_new_fromFtoE"    : {"Value": CL_new_fromFtoE,     "Unit": "Non dimensional"},
                    "CL_wb_fromEto0"     : {"Value": CL_wb_fromEto0,      "Unit": "Non dimensional"},
                    "CL_ht_fromEto0"     : {"Value": CL_ht_fromEto0,      "Unit": "Non dimensional"},
                    "CL_new_fromEto0"    : {"Value": CL_new_fromEto0,     "Unit": "Non dimensional"},

                    "L_wb_from0toS"     : {"Value": L_wb_from0toS,      "Unit": "daN"},
                    "L_ht_from0toS"     : {"Value": L_ht_from0toS,      "Unit": "daN"},
                    "L_new_from0toS"    : {"Value": L_new_from0toS,     "Unit": "daN"},
                    "L_wb_fromStoA"     : {"Value": L_wb_fromStoA,      "Unit": "daN"},
                    "L_ht_fromStoA"     : {"Value": L_ht_fromStoA,      "Unit": "daN"},
                    "L_new_fromStoA"    : {"Value": L_new_fromStoA,     "Unit": "daN"},
                    "L_wb_fromAtoC"     : {"Value": L_wb_fromAtoC,      "Unit": "daN"},
                    "L_ht_fromAtoC"     : {"Value": L_ht_fromAtoC,      "Unit": "daN"},
                    "L_new_fromAtoC"    : {"Value": L_new_fromAtoC,     "Unit": "daN"},
                    "L_wb_fromCtoD"     : {"Value": L_wb_fromCtoD,      "Unit": "daN"},
                    "L_ht_fromCtoD"     : {"Value": L_ht_fromCtoD,      "Unit": "daN"},
                    "L_new_fromCtoD"    : {"Value": L_new_fromCtoD,     "Unit": "daN"},
                    "L_wb_fromDto0"     : {"Value": L_wb_fromDto0,      "Unit": "daN"},
                    "L_ht_fromDto0"     : {"Value": L_ht_fromDto0,      "Unit": "daN"},
                    "L_new_fromDto0"    : {"Value": L_new_fromDto0,     "Unit": "daN"},
                    "L_wb_from0toSinv"  : {"Value": L_wb_from0toSinv,   "Unit": "daN"},
                    "L_ht_from0toSinv"  : {"Value": L_ht_from0toSinv,   "Unit": "daN"},
                    "L_new_from0toSinv" : {"Value": L_new_from0toSinv,  "Unit": "daN"},
                    "L_wb_fromSinvtoG"  : {"Value": L_wb_fromSinvtoG,   "Unit": "daN"},
                    "L_ht_fromSinvtoG"  : {"Value": L_ht_fromSinvtoG,   "Unit": "daN"},
                    "L_new_fromSinvtoG" : {"Value": L_new_fromSinvtoG,  "Unit": "daN"},
                    "L_wb_fromGtoGust1" : {"Value": L_wb_fromGtoGust1,  "Unit": "daN"},
                    "L_ht_fromGtoGust1" : {"Value": L_ht_fromGtoGust1,  "Unit": "daN"},
                    "L_new_fromGtoGust1": {"Value": L_new_fromGtoGust1, "Unit": "daN"},
                    "L_wb_fromGust1toF" : {"Value": L_wb_fromGust1toF,  "Unit": "daN"},
                    "L_ht_fromGust1toF" : {"Value": L_ht_fromGust1toF,  "Unit": "daN"},
                    "L_new_fromGust1toF": {"Value": L_new_fromGust1toF, "Unit": "daN"},
                    "L_wb_fromFtoE"     : {"Value": L_wb_fromFtoE,      "Unit": "daN"},
                    "L_ht_fromFtoE"     : {"Value": L_ht_fromFtoE,      "Unit": "daN"},
                    "L_new_fromFtoE"    : {"Value": L_new_fromFtoE,     "Unit": "daN"},
                    "L_wb_fromEto0"     : {"Value": L_wb_fromEto0,      "Unit": "daN"},
                    "L_ht_fromEto0"     : {"Value": L_ht_fromEto0,      "Unit": "daN"},
                    "L_new_fromEto0"    : {"Value": L_new_fromEto0,     "Unit": "daN"},

                    "CM_due_to_CL_from0toS"     : {"Value": CM_due_to_CL_from0toS,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS"     : {"Value": CM_due_to_CD_from0toS,      "Unit": "Non dimensional"},
                    "CM_CG_from0toS"            : {"Value": CM_CG_from0toS,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA"     : {"Value": CM_due_to_CL_fromStoA,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA"     : {"Value": CM_due_to_CD_fromStoA,      "Unit": "Non dimensional"},
                    "CM_CG_fromStoA"            : {"Value": CM_CG_fromStoA,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoC"     : {"Value": CM_due_to_CL_fromAtoC,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoC"     : {"Value": CM_due_to_CD_fromAtoC,      "Unit": "Non dimensional"},
                    "CM_CG_fromAtoC"            : {"Value": CM_CG_fromAtoC,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoD"     : {"Value": CM_due_to_CL_fromCtoD,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoD"     : {"Value": CM_due_to_CD_fromCtoD,      "Unit": "Non dimensional"},
                    "CM_CG_fromCtoD"            : {"Value": CM_CG_fromCtoD,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0"     : {"Value": CM_due_to_CL_fromDto0,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0"     : {"Value": CM_due_to_CD_fromDto0,      "Unit": "Non dimensional"},
                    "CM_CG_fromDto0"            : {"Value": CM_CG_fromDto0,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv"  : {"Value": CM_due_to_CL_from0toSinv,   "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv"  : {"Value": CM_due_to_CD_from0toSinv,   "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv"         : {"Value": CM_CG_from0toSinv,          "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG"  : {"Value": CM_due_to_CL_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG"  : {"Value": CM_due_to_CD_fromSinvtoG,   "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG"         : {"Value": CM_CG_fromSinvtoG,          "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoGust1" : {"Value": CM_due_to_CL_fromGtoGust1,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoGust1" : {"Value": CM_due_to_CD_fromGtoGust1,  "Unit": "Non dimensional"},
                    "CM_CG_fromGtoGust1"        : {"Value": CM_CG_fromGtoGust1,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toF" : {"Value": CM_due_to_CL_fromGust1toF,  "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toF" : {"Value": CM_due_to_CD_fromGust1toF,  "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toF"        : {"Value": CM_CG_fromGust1toF,         "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE"     : {"Value": CM_due_to_CL_fromFtoE,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE"     : {"Value": CM_due_to_CD_fromFtoE,      "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE"            : {"Value": CM_CG_fromFtoE,             "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0"     : {"Value": CM_due_to_CL_fromEto0,      "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0"     : {"Value": CM_due_to_CD_fromEto0,      "Unit": "Non dimensional"},
                    "CM_CG_fromEto0"            : {"Value": CM_CG_fromEto0,             "Unit": "Non dimensional"}
                    }
                }

        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case2_inverted"):
            
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS":     {"Value": V_fe_from0toS, "Unit": "m/s"},
                    "n_fe_from0toS":     {"Value": n_fe_from0toS, "Unit": "g"},
                    "V_fe_fromStoA":     {"Value": V_fe_fromStoA, "Unit": "m/s"},
                    "n_fe_fromStoA":     {"Value": n_fe_fromStoA, "Unit": "g"},
                    "V_fe_fromAtoC":     {"Value": V_fe_fromAtoC, "Unit": "m/s"},
                    "n_fe_fromAtoC":     {"Value": n_fe_fromAtoC, "Unit": "g"},
                    "V_fe_fromCtoD":     {"Value": V_fe_fromCtoD, "Unit": "m/s"},
                    "n_fe_fromCtoD":     {"Value": n_fe_fromCtoD, "Unit": "g"},
                    "V_fe_fromDto0":     {"Value": V_fe_fromDto0, "Unit": "m/s"},
                    "n_fe_fromDto0":     {"Value": n_fe_fromDto0, "Unit": "g"},
                    "V_fe_from0toSinv":  {"Value": V_fe_from0toSinv, "Unit": "m/s"},
                    "n_fe_from0toSinv":  {"Value": n_fe_from0toSinv, "Unit": "g"},
                    "V_fe_fromSinvtoG":  {"Value": V_fe_fromSinvtoG, "Unit": "m/s"},
                    "n_fe_fromSinvtoG":  {"Value": n_fe_fromSinvtoG, "Unit": "g"},
                    "V_fe_fromGtoF":     {"Value": V_fe_fromGtoF, "Unit": "m/s"},
                    "n_fe_fromGtoF":     {"Value": n_fe_fromGtoF, "Unit": "g"},
                    "V_fe_fromFtoE":     {"Value": V_fe_fromFtoE, "Unit": "m/s"},
                    "n_fe_fromFtoE":     {"Value": n_fe_fromFtoE, "Unit": "g"},
                    "V_fe_fromEto0":     {"Value": V_fe_fromEto0, "Unit": "m/s"},
                    "n_fe_fromEto0":     {"Value": n_fe_fromEto0, "Unit": "g"},

                    "q_from0toS":     {"Value": q_from0toS, "Unit": "Pa"},
                    "CD_from0toS":     {"Value": CD_from0toS, "Unit": "Non dimensional"},
                    "q_fromStoA":     {"Value": q_fromStoA, "Unit": "Pa"},
                    "CD_fromStoA":     {"Value": CD_fromStoA, "Unit": "Non dimensional"},
                    "q_fromAtoC":     {"Value": q_fromAtoC, "Unit": "Pa"},
                    "CD_fromAtoC":     {"Value": CD_fromAtoC, "Unit": "Non dimensional"},
                    "q_fromCtoD":     {"Value": q_fromCtoD, "Unit": "Pa"},
                    "CD_fromCtoD":     {"Value": CD_fromCtoD, "Unit": "Non dimensional"},
                    "q_fromDto0":     {"Value": q_fromDto0, "Unit": "Pa"},
                    "CD_fromDto0":     {"Value": CD_fromDto0, "Unit": "Non dimensional"},
                    "q_from0toSinv":  {"Value": q_from0toSinv, "Unit": "Pa"},
                    "CD_from0toSinv":  {"Value": CD_from0toSinv, "Unit": "Non dimensional"},
                    "q_fromSinvtoG":  {"Value": q_fromSinvtoG, "Unit": "Pa"},
                    "CD_fromSinvtoG":  {"Value": CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "q_fromGtoF":     {"Value": q_fromGtoF, "Unit": "Pa"},
                    "CD_fromGtoF":     {"Value": CD_fromGtoF, "Unit": "Non dimensional"},
                    "q_fromFtoE":     {"Value": q_fromFtoE, "Unit": "Pa"},
                    "CD_fromFtoE":     {"Value": CD_fromFtoE, "Unit": "Non dimensional"},
                    "q_fromEto0":     {"Value": q_fromEto0, "Unit": "Pa"},
                    "CD_fromEto0":     {"Value": CD_fromEto0, "Unit": "Non dimensional"},

                    "alfa_from0toS":     {"Value": alfa_from0toS, "Unit": "deg"},
                    "alfa_new_from0toS":     {"Value": alfa_new_from0toS, "Unit": "deg"},
                    "alfa_fromStoA":     {"Value": alfa_fromStoA, "Unit": "deg"},
                    "alfa_new_fromStoA":     {"Value": alfa_new_fromStoA, "Unit": "deg"},
                    "alfa_fromAtoC":     {"Value": alfa_fromAtoC, "Unit": "deg"},
                    "alfa_new_fromAtoC":     {"Value": alfa_new_fromAtoC, "Unit": "deg"},
                    "alfa_fromCtoD":     {"Value": alfa_fromCtoD, "Unit": "deg"},
                    "alfa_new_fromCtoD":     {"Value": alfa_new_fromCtoD, "Unit": "deg"},
                    "alfa_fromDto0":     {"Value": alfa_fromDto0, "Unit": "deg"},
                    "alfa_new_fromDto0":     {"Value": alfa_new_fromDto0, "Unit": "deg"},
                    "alfa_from0toSinv":  {"Value": alfa_from0toSinv, "Unit": "deg"},
                    "alfa_new_from0toSinv":  {"Value": alfa_new_from0toSinv, "Unit": "deg"},
                    "alfa_fromSinvtoG":  {"Value": alfa_fromSinvtoG, "Unit": "deg"},
                    "alfa_new_fromSinvtoG":  {"Value": alfa_new_fromSinvtoG, "Unit": "deg"},
                    "alfa_fromGtoF":     {"Value": alfa_fromGtoF, "Unit": "deg"},
                    "alfa_new_fromGtoF":     {"Value": alfa_new_fromGtoF, "Unit": "deg"},
                    "alfa_fromFtoE":     {"Value": alfa_fromFtoE, "Unit": "deg"},
                    "alfa_new_fromFtoE":     {"Value": alfa_new_fromFtoE, "Unit": "deg"},
                    "alfa_fromEto0":     {"Value": alfa_fromEto0, "Unit": "deg"},
                    "alfa_new_fromEto0":     {"Value": alfa_new_fromEto0, "Unit": "deg"},

                    "CL_wb_from0toS":     {"Value": CL_wb_from0toS, "Unit": "Non dimensional"},
                    "CL_ht_from0toS":     {"Value": CL_ht_from0toS, "Unit": "Non dimensional"},
                    "CL_new_from0toS":     {"Value": CL_new_from0toS, "Unit": "Non dimensional"},
                    "CL_wb_fromStoA":     {"Value": CL_wb_fromStoA, "Unit": "Non dimensional"},
                    "CL_ht_fromStoA":     {"Value": CL_ht_fromStoA, "Unit": "Non dimensional"},
                    "CL_new_fromStoA":     {"Value": CL_new_fromStoA, "Unit": "Non dimensional"},
                    "CL_wb_fromAtoC":     {"Value": CL_wb_fromAtoC, "Unit": "Non dimensional"},
                    "CL_ht_fromAtoC":     {"Value": CL_ht_fromAtoC, "Unit": "Non dimensional"},
                    "CL_new_fromAtoC":     {"Value": CL_new_fromAtoC, "Unit": "Non dimensional"},
                    "CL_wb_fromCtoD":     {"Value": CL_wb_fromCtoD, "Unit": "Non dimensional"},
                    "CL_ht_fromCtoD":     {"Value": CL_ht_fromCtoD, "Unit": "Non dimensional"},
                    "CL_new_fromCtoD":     {"Value": CL_new_fromCtoD, "Unit": "Non dimensional"},
                    "CL_wb_fromDto0":     {"Value": CL_wb_fromDto0, "Unit": "Non dimensional"},
                    "CL_ht_fromDto0":     {"Value": CL_ht_fromDto0, "Unit": "Non dimensional"},
                    "CL_new_fromDto0":     {"Value": CL_new_fromDto0, "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv":  {"Value": CL_wb_from0toSinv, "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv":  {"Value": CL_ht_from0toSinv, "Unit": "Non dimensional"},
                    "CL_new_from0toSinv":  {"Value": CL_new_from0toSinv, "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG":  {"Value": CL_wb_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG":  {"Value": CL_ht_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG":  {"Value": CL_new_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_wb_fromGtoF":     {"Value": CL_wb_fromGtoF, "Unit": "Non dimensional"},
                    "CL_ht_fromGtoF":     {"Value": CL_ht_fromGtoF, "Unit": "Non dimensional"},
                    "CL_new_fromGtoF":     {"Value": CL_new_fromGtoF, "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE":     {"Value": CL_wb_fromFtoE, "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE":     {"Value": CL_ht_fromFtoE, "Unit": "Non dimensional"},
                    "CL_new_fromFtoE":     {"Value": CL_new_fromFtoE, "Unit": "Non dimensional"},
                    "CL_wb_fromEto0":     {"Value": CL_wb_fromEto0, "Unit": "Non dimensional"},
                    "CL_ht_fromEto0":     {"Value": CL_ht_fromEto0, "Unit": "Non dimensional"},
                    "CL_new_fromEto0":     {"Value": CL_new_fromEto0, "Unit": "Non dimensional"},    

                    "L_wb_from0toS":     {"Value": L_wb_from0toS, "Unit": "daN"},
                    "L_ht_from0toS":     {"Value": L_ht_from0toS, "Unit": "daN"},
                    "L_new_from0toS":     {"Value": L_new_from0toS, "Unit": "daN"},
                    "L_wb_fromStoA":     {"Value": L_wb_fromStoA, "Unit": "daN"},
                    "L_ht_fromStoA":     {"Value": L_ht_fromStoA, "Unit": "daN"},
                    "L_new_fromStoA":     {"Value": L_new_fromStoA, "Unit": "daN"},
                    "L_wb_fromAtoC":     {"Value": L_wb_fromAtoC, "Unit": "daN"},
                    "L_ht_fromAtoC":     {"Value": L_ht_fromAtoC, "Unit": "daN"},
                    "L_new_fromAtoC":     {"Value": L_new_fromAtoC, "Unit": "daN"},
                    "L_wb_fromCtoD":     {"Value": L_wb_fromCtoD, "Unit": "daN"},
                    "L_ht_fromCtoD":     {"Value": L_ht_fromCtoD, "Unit": "daN"},
                    "L_new_fromCtoD":     {"Value": L_new_fromCtoD, "Unit": "daN"},
                    "L_wb_fromDto0":     {"Value": L_wb_fromDto0, "Unit": "daN"},
                    "L_ht_fromDto0":     {"Value": L_ht_fromDto0, "Unit": "daN"},
                    "L_new_fromDto0":     {"Value": L_new_fromDto0, "Unit": "daN"},
                    "L_wb_from0toSinv":  {"Value": L_wb_from0toSinv, "Unit": "daN"},
                    "L_ht_from0toSinv":  {"Value": L_ht_from0toSinv, "Unit": "daN"},
                    "L_new_from0toSinv":  {"Value": L_new_from0toSinv, "Unit": "daN"},
                    "L_wb_fromSinvtoG":  {"Value": L_wb_fromSinvtoG, "Unit": "daN"},
                    "L_ht_fromSinvtoG":  {"Value": L_ht_fromSinvtoG, "Unit": "daN"},
                    "L_new_fromSinvtoG":  {"Value": L_new_fromSinvtoG, "Unit": "daN"},
                    "L_wb_fromGtoF":     {"Value": L_wb_fromGtoF, "Unit": "daN"},
                    "L_ht_fromGtoF":     {"Value": L_ht_fromGtoF, "Unit": "daN"},
                    "L_new_fromGtoF":     {"Value": L_new_fromGtoF, "Unit": "daN"},
                    "L_wb_fromFtoE":     {"Value": L_wb_fromFtoE, "Unit": "daN"},
                    "L_ht_fromFtoE":     {"Value": L_ht_fromFtoE, "Unit": "daN"},
                    "L_new_fromFtoE":     {"Value": L_new_fromFtoE, "Unit": "daN"},
                    "L_wb_fromEto0":     {"Value": L_wb_fromEto0, "Unit": "daN"},
                    "L_ht_fromEto0":     {"Value": L_ht_fromEto0, "Unit": "daN"},
                    "L_new_fromEto0":     {"Value": L_new_fromEto0, "Unit": "daN"},

                    "CM_due_to_CL_from0toS":     {"Value": CM_due_to_CL_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS":     {"Value": CM_due_to_CD_from0toS, "Unit": "Non dimensional"},
                    "CM_CG_from0toS":     {"Value": CM_CG_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA":     {"Value": CM_due_to_CL_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA":     {"Value": CM_due_to_CD_fromStoA, "Unit": "Non dimensional"},
                    "CM_CG_fromStoA":     {"Value": CM_CG_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoC":     {"Value": CM_due_to_CL_fromAtoC, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoC":     {"Value": CM_due_to_CD_fromAtoC, "Unit": "Non dimensional"},
                    "CM_CG_fromAtoC":     {"Value": CM_CG_fromAtoC, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoD":     {"Value": CM_due_to_CL_fromCtoD, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoD":     {"Value": CM_due_to_CD_fromCtoD, "Unit": "Non dimensional"},
                    "CM_CG_fromCtoD":     {"Value": CM_CG_fromCtoD, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0":     {"Value": CM_due_to_CL_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0":     {"Value": CM_due_to_CD_fromDto0, "Unit": "Non dimensional"},
                    "CM_CG_fromDto0":     {"Value": CM_CG_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv":  {"Value": CM_due_to_CL_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv":  {"Value": CM_due_to_CD_from0toSinv, "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv":  {"Value": CM_CG_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG":  {"Value": CM_due_to_CL_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG":  {"Value": CM_due_to_CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG":  {"Value": CM_CG_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoF":     {"Value": CM_due_to_CL_fromGtoF, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoF":     {"Value": CM_due_to_CD_fromGtoF, "Unit": "Non dimensional"},
                    "CM_CG_fromGtoF":     {"Value": CM_CG_fromGtoF, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE":     {"Value": CM_due_to_CL_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE":     {"Value": CM_due_to_CD_fromFtoE, "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE":     {"Value": CM_CG_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0":     {"Value": CM_due_to_CL_fromEto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0":     {"Value": CM_due_to_CD_fromEto0, "Unit": "Non dimensional"},
                    "CM_CG_fromEto0":     {"Value": CM_CG_fromEto0, "Unit": "Non dimensional"}  

                    }
                }

        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case3_inverted"):
            
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS":     {"Value": V_fe_from0toS, "Unit": "m/s"},
                    "n_fe_from0toS":     {"Value": n_fe_from0toS, "Unit": "g"},
                    "V_fe_fromStoA":     {"Value": V_fe_fromStoA, "Unit": "m/s"},
                    "n_fe_fromStoA":     {"Value": n_fe_fromStoA, "Unit": "g"},
                    "V_fe_fromAtoC":     {"Value": V_fe_fromAtoC, "Unit": "m/s"},
                    "n_fe_fromAtoC":     {"Value": n_fe_fromAtoC, "Unit": "g"},
                    "V_fe_fromCtoD":     {"Value": V_fe_fromCtoD, "Unit": "m/s"},
                    "n_fe_fromCtoD":     {"Value": n_fe_fromCtoD, "Unit": "g"},
                    "V_fe_fromDto0":     {"Value": V_fe_fromDto0, "Unit": "m/s"},
                    "n_fe_fromDto0":     {"Value": n_fe_fromDto0, "Unit": "g"},
                    "V_fe_from0toSinv":  {"Value": V_fe_from0toSinv, "Unit": "m/s"},
                    "n_fe_from0toSinv":  {"Value": n_fe_from0toSinv, "Unit": "g"},
                    "V_fe_fromSinvtoG":  {"Value": V_fe_fromSinvtoG, "Unit": "m/s"},
                    "n_fe_fromSinvtoG":  {"Value": n_fe_fromSinvtoG, "Unit": "g"},
                    "V_fe_fromGtoF":     {"Value": V_fe_fromGtoF, "Unit": "m/s"},
                    "n_fe_fromGtoF":     {"Value": n_fe_fromGtoF, "Unit": "g"},
                    "V_fe_fromFtoE"    : {"Value": V_fe_fromFtoE, "Unit": "m/s"},
                    "n_fe_fromFtoE"    : {"Value": n_fe_fromFtoE, "Unit": "g"},
                    "V_fe_fromEto0"    : {"Value": V_fe_fromEto0, "Unit": "m/s"},
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"},

                    "q_from0toS":     {"Value": q_from0toS, "Unit": "Pa"},
                    "CD_from0toS":     {"Value": CD_from0toS, "Unit": "Non dimensional"},
                    "q_fromStoA":     {"Value": q_fromStoA, "Unit": "Pa"},
                    "CD_fromStoA":     {"Value": CD_fromStoA, "Unit": "Non dimensional"},
                    "q_fromAtoC":     {"Value": q_fromAtoC, "Unit": "Pa"},
                    "CD_fromAtoC":     {"Value": CD_fromAtoC, "Unit": "Non dimensional"},
                    "q_fromCtoD":     {"Value": q_fromCtoD, "Unit": "Pa"},
                    "CD_fromCtoD":     {"Value": CD_fromCtoD, "Unit": "Non dimensional"},
                    "q_fromDto0":     {"Value": q_fromDto0, "Unit": "Pa"},
                    "CD_fromDto0":     {"Value": CD_fromDto0, "Unit": "Non dimensional"},
                    "q_from0toSinv":  {"Value": q_from0toSinv, "Unit": "Pa"},
                    "CD_from0toSinv":  {"Value": CD_from0toSinv, "Unit": "Non dimensional"},
                    "q_fromSinvtoG":  {"Value": q_fromSinvtoG, "Unit": "Pa"},
                    "CD_fromSinvtoG":  {"Value": CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "q_fromGtoF":     {"Value": q_fromGtoF, "Unit": "Pa"},
                    "CD_fromGtoF":     {"Value": CD_fromGtoF, "Unit": "Non dimensional"},
                    "q_fromFtoE"    : {"Value": q_fromFtoE, "Unit": "Pa"},
                    "CD_fromFtoE"    : {"Value": CD_fromFtoE, "Unit": "Non dimensional"},
                    "q_fromEto0"    : {"Value": q_fromEto0, "Unit": "Pa"},
                    "CD_fromEto0"    : {"Value": CD_fromEto0, "Unit": "Non dimensional"},
                    
                    "alfa_from0toS":     {"Value": alfa_from0toS, "Unit": "deg"},
                    "alfa_new_from0toS":     {"Value": alfa_new_from0toS, "Unit": "deg"},
                    "alfa_fromStoA":     {"Value": alfa_fromStoA, "Unit": "deg"},
                    "alfa_new_fromStoA":     {"Value": alfa_new_fromStoA, "Unit": "deg"},
                    "alfa_fromAtoC":     {"Value": alfa_fromAtoC, "Unit": "deg"},
                    "alfa_new_fromAtoC":     {"Value": alfa_new_fromAtoC, "Unit": "deg"},
                    "alfa_fromCtoD":     {"Value": alfa_fromCtoD, "Unit": "deg"},
                    "alfa_new_fromCtoD":     {"Value": alfa_new_fromCtoD, "Unit": "deg"},
                    "alfa_fromDto0":     {"Value": alfa_fromDto0, "Unit": "deg"},
                    "alfa_new_fromDto0":     {"Value": alfa_new_fromDto0, "Unit": "deg"},
                    "alfa_from0toSinv":  {"Value": alfa_from0toSinv, "Unit": "deg"},
                    "alfa_new_from0toSinv":  {"Value": alfa_new_from0toSinv, "Unit": "deg"},
                    "alfa_fromSinvtoG":  {"Value": alfa_fromSinvtoG, "Unit": "deg"},
                    "alfa_new_fromSinvtoG":  {"Value": alfa_new_fromSinvtoG, "Unit": "deg"},
                    "alfa_fromGtoF":     {"Value": alfa_fromGtoF, "Unit": "deg"},
                    "alfa_new_fromGtoF":     {"Value": alfa_new_fromGtoF, "Unit": "deg"},
                    "alfa_fromFtoE"    : {"Value": alfa_fromFtoE, "Unit": "deg"},
                    "alfa_new_fromFtoE"    : {"Value": alfa_new_fromFtoE, "Unit": "deg"},
                    "alfa_fromEto0"    : {"Value": alfa_fromEto0, "Unit": "deg"},
                    "alfa_new_fromEto0"    : {"Value": alfa_new_fromEto0, "Unit": "deg"},
                    
                    "CL_wb_from0toS":     {"Value": CL_wb_from0toS, "Unit": "Non dimensional"},
                    "CL_ht_from0toS":     {"Value": CL_ht_from0toS, "Unit": "Non dimensional"},
                    "CL_new_from0toS":     {"Value": CL_new_from0toS, "Unit": "Non dimensional"},
                    "CL_wb_fromStoA":     {"Value": CL_wb_fromStoA, "Unit": "Non dimensional"},
                    "CL_ht_fromStoA":     {"Value": CL_ht_fromStoA, "Unit": "Non dimensional"},
                    "CL_new_fromStoA":     {"Value": CL_new_fromStoA, "Unit": "Non dimensional"},
                    "CL_wb_fromAtoC":     {"Value": CL_wb_fromAtoC, "Unit": "Non dimensional"},
                    "CL_ht_fromAtoC":     {"Value": CL_ht_fromAtoC, "Unit": "Non dimensional"},
                    "CL_new_fromAtoC":     {"Value": CL_new_fromAtoC, "Unit": "Non dimensional"},
                    "CL_wb_fromCtoD":     {"Value": CL_wb_fromCtoD, "Unit": "Non dimensional"},
                    "CL_ht_fromCtoD":     {"Value": CL_ht_fromCtoD, "Unit": "Non dimensional"},
                    "CL_new_fromCtoD":     {"Value": CL_new_fromCtoD, "Unit": "Non dimensional"},
                    "CL_wb_fromDto0":     {"Value": CL_wb_fromDto0, "Unit": "Non dimensional"},
                    "CL_ht_fromDto0":     {"Value": CL_ht_fromDto0, "Unit": "Non dimensional"},
                    "CL_new_fromDto0":     {"Value": CL_new_fromDto0, "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv":  {"Value": CL_wb_from0toSinv, "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv":  {"Value": CL_ht_from0toSinv, "Unit": "Non dimensional"},
                    "CL_new_from0toSinv":  {"Value": CL_new_from0toSinv, "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG":  {"Value": CL_wb_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG":  {"Value": CL_ht_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG":  {"Value": CL_new_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_wb_fromGtoF":     {"Value": CL_wb_fromGtoF, "Unit": "Non dimensional"},
                    "CL_ht_fromGtoF":     {"Value": CL_ht_fromGtoF, "Unit": "Non dimensional"},
                    "CL_new_fromGtoF":     {"Value": CL_new_fromGtoF, "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE"    : {"Value": CL_wb_fromFtoE, "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE"    : {"Value": CL_ht_fromFtoE, "Unit": "Non dimensional"},
                    "CL_new_fromFtoE"    : {"Value": CL_new_fromFtoE, "Unit": "Non dimensional"},
                    "CL_wb_fromEto0"    : {"Value": CL_wb_fromEto0, "Unit": "Non dimensional"},
                    "CL_ht_fromEto0"    : {"Value": CL_ht_fromEto0, "Unit": "Non dimensional"},
                    "CL_new_fromEto0"    : {"Value": CL_new_fromEto0, "Unit": "Non dimensional"},
                    
                    "L_wb_from0toS":     {"Value": L_wb_from0toS, "Unit": "daN"},
                    "L_ht_from0toS":     {"Value": L_ht_from0toS, "Unit": "daN"},
                    "L_new_from0toS":     {"Value": L_new_from0toS, "Unit": "daN"},
                    "L_wb_fromStoA":     {"Value": L_wb_fromStoA, "Unit": "daN"},
                    "L_ht_fromStoA":     {"Value": L_ht_fromStoA, "Unit": "daN"},
                    "L_new_fromStoA":     {"Value": L_new_fromStoA, "Unit": "daN"},
                    "L_wb_fromAtoC":     {"Value": L_wb_fromAtoC, "Unit": "daN"},
                    "L_ht_fromAtoC":     {"Value": L_ht_fromAtoC, "Unit": "daN"},
                    "L_new_fromAtoC":     {"Value": L_new_fromAtoC, "Unit": "daN"},
                    "L_wb_fromCtoD":     {"Value": L_wb_fromCtoD, "Unit": "daN"},
                    "L_ht_fromCtoD":     {"Value": L_ht_fromCtoD, "Unit": "daN"},
                    "L_new_fromCtoD":     {"Value": L_new_fromCtoD, "Unit": "daN"},
                    "L_wb_fromDto0":     {"Value": L_wb_fromDto0, "Unit": "daN"},
                    "L_ht_fromDto0":     {"Value": L_ht_fromDto0, "Unit": "daN"},
                    "L_new_fromDto0":     {"Value": L_new_fromDto0, "Unit": "daN"},
                    "L_wb_from0toSinv":  {"Value": L_wb_from0toSinv, "Unit": "daN"},
                    "L_ht_from0toSinv":  {"Value": L_ht_from0toSinv, "Unit": "daN"},
                    "L_new_from0toSinv":  {"Value": L_new_from0toSinv, "Unit": "daN"},
                    "L_wb_fromSinvtoG":  {"Value": L_wb_fromSinvtoG, "Unit": "daN"},
                    "L_ht_fromSinvtoG":  {"Value": L_ht_fromSinvtoG, "Unit": "daN"},
                    "L_new_fromSinvtoG":  {"Value": L_new_fromSinvtoG, "Unit": "daN"},
                    "L_wb_fromGtoF":     {"Value": L_wb_fromGtoF, "Unit": "daN"},
                    "L_ht_fromGtoF":     {"Value": L_ht_fromGtoF, "Unit": "daN"},
                    "L_new_fromGtoF":     {"Value": L_new_fromGtoF, "Unit": "daN"},
                    "L_wb_fromFtoE"    : {"Value": L_wb_fromFtoE, "Unit": "daN"},
                    "L_ht_fromFtoE"    : {"Value": L_ht_fromFtoE, "Unit": "daN"},
                    "L_new_fromFtoE"    : {"Value": L_new_fromFtoE, "Unit": "daN"},
                    "L_wb_fromEto0"    : {"Value": L_wb_fromEto0, "Unit": "daN"},
                    "L_ht_fromEto0"    : {"Value": L_ht_fromEto0, "Unit": "daN"},
                    "L_new_fromEto0"    : {"Value": L_new_fromEto0, "Unit": "daN"},
                                        
                    "CM_due_to_CL_from0toS":     {"Value": CM_due_to_CL_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS":     {"Value": CM_due_to_CD_from0toS, "Unit": "Non dimensional"},
                    "CM_CG_from0toS":     {"Value": CM_CG_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA":     {"Value": CM_due_to_CL_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA":     {"Value": CM_due_to_CD_fromStoA, "Unit": "Non dimensional"},
                    "CM_CG_fromStoA":     {"Value": CM_CG_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoC":     {"Value": CM_due_to_CL_fromAtoC, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoC":     {"Value": CM_due_to_CD_fromAtoC, "Unit": "Non dimensional"},
                    "CM_CG_fromAtoC":     {"Value": CM_CG_fromAtoC, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoD":     {"Value": CM_due_to_CL_fromCtoD, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoD":     {"Value": CM_due_to_CD_fromCtoD, "Unit": "Non dimensional"},
                    "CM_CG_fromCtoD":     {"Value": CM_CG_fromCtoD, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0":     {"Value": CM_due_to_CL_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0":     {"Value": CM_due_to_CD_fromDto0, "Unit": "Non dimensional"},
                    "CM_CG_fromDto0":     {"Value": CM_CG_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv":  {"Value": CM_due_to_CL_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv":  {"Value": CM_due_to_CD_from0toSinv, "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv":  {"Value": CM_CG_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG":  {"Value": CM_due_to_CL_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG":  {"Value": CM_due_to_CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG":  {"Value": CM_CG_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoF":     {"Value": CM_due_to_CL_fromGtoF, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoF":     {"Value": CM_due_to_CD_fromGtoF, "Unit": "Non dimensional"},
                    "CM_CG_fromGtoF":     {"Value": CM_CG_fromGtoF, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE"    : {"Value": CM_due_to_CL_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE"    : {"Value": CM_due_to_CD_fromFtoE, "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE"    : {"Value": CM_CG_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0"    : {"Value": CM_due_to_CL_fromEto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0"    : {"Value": CM_due_to_CD_fromEto0, "Unit": "Non dimensional"},
                    "CM_CG_fromEto0"    : {"Value": CM_CG_fromEto0, "Unit": "Non dimensional"}
                    
                    }
                }

        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case1_inverted"):
            
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS": {"Value": V_fe_from0toS, "Unit": "m/s"},
                    "n_fe_from0toS": {"Value": n_fe_from0toS, "Unit": "g"},
                    "V_fe_fromStoA": {"Value": V_fe_fromStoA, "Unit": "m/s"},
                    "n_fe_fromStoA": {"Value": n_fe_fromStoA, "Unit": "g"},
                    "V_fe_fromAtoGust1": {"Value": V_fe_fromAtoGust1, "Unit": "m/s"},
                    "n_fe_fromAtoGust1": {"Value": n_fe_fromAtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toC": {"Value": V_fe_fromGust1toC, "Unit": "m/s"},
                    "n_fe_fromGust1toC": {"Value": n_fe_fromGust1toC, "Unit": "g"},
                    "V_fe_fromCtoGust2": {"Value": V_fe_fromCtoGust2, "Unit": "m/s"},
                    "n_fe_fromCtoGust2": {"Value": n_fe_fromCtoGust2, "Unit": "g"},
                    "V_fe_fromGust2toD": {"Value": V_fe_fromGust2toD, "Unit": "m/s"},
                    "n_fe_fromGust2toD": {"Value": n_fe_fromGust2toD, "Unit": "g"},
                    "V_fe_fromDto0":     {"Value": V_fe_fromDto0, "Unit": "m/s"},
                    "n_fe_fromDto0":     {"Value": n_fe_fromDto0, "Unit": "g"},
                    "V_fe_from0toSinv":  {"Value": V_fe_from0toSinv, "Unit": "m/s"},
                    "n_fe_from0toSinv":  {"Value": n_fe_from0toSinv, "Unit": "g"},
                    "V_fe_fromSinvtoG":  {"Value": V_fe_fromSinvtoG, "Unit": "m/s"},
                    "n_fe_fromSinvtoG":  {"Value": n_fe_fromSinvtoG, "Unit": "g"},
                    "V_fe_fromGtoGust1": {"Value": V_fe_fromGtoGust1, "Unit": "m/s"},
                    "n_fe_fromGtoGust1": {"Value": n_fe_fromGtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toF": {"Value": V_fe_fromGust1toF, "Unit": "m/s"},
                    "n_fe_fromGust1toF": {"Value": n_fe_fromGust1toF, "Unit": "g"},
                    "V_fe_fromFtoE"    : {"Value": V_fe_fromFtoE, "Unit": "m/s"},
                    "n_fe_fromFtoE"    : {"Value": n_fe_fromFtoE, "Unit": "g"},
                    "V_fe_fromEto0"    : {"Value": V_fe_fromEto0, "Unit": "m/s"},
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"},

                    "q_from0toS": {"Value": q_from0toS, "Unit": "Pa"},
                    "CD_from0toS": {"Value": CD_from0toS, "Unit": "Non dimensional"},
                    "q_fromStoA": {"Value": q_fromStoA, "Unit": "Pa"},
                    "CD_fromStoA": {"Value": CD_fromStoA, "Unit": "Non dimensional"},
                    "q_fromAtoGust1": {"Value": q_fromAtoGust1, "Unit": "Pa"},
                    "CD_fromAtoGust1": {"Value": CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toC": {"Value": q_fromGust1toC, "Unit": "Pa"},
                    "CD_fromGust1toC": {"Value": CD_fromGust1toC, "Unit": "Non dimensional"},
                    "q_fromCtoGust2": {"Value": q_fromCtoGust2, "Unit": "Pa"},
                    "CD_fromCtoGust2": {"Value": CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "q_fromGust2toD": {"Value": q_fromGust2toD, "Unit": "Pa"},
                    "CD_fromGust2toD": {"Value": CD_fromGust2toD, "Unit": "Non dimensional"},
                    "q_fromDto0":     {"Value": q_fromDto0, "Unit": "Pa"},
                    "CD_fromDto0":     {"Value": CD_fromDto0, "Unit": "Non dimensional"},
                    "q_from0toSinv":  {"Value": q_from0toSinv, "Unit": "Pa"},
                    "CD_from0toSinv":  {"Value": CD_from0toSinv, "Unit": "Non dimensional"},
                    "q_fromSinvtoG":  {"Value": q_fromSinvtoG, "Unit": "Pa"},
                    "CD_fromSinvtoG":  {"Value": CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "q_fromGtoGust1": {"Value": q_fromGtoGust1, "Unit": "Pa"},
                    "CD_fromGtoGust1": {"Value": CD_fromGtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toF": {"Value": q_fromGust1toF, "Unit": "Pa"},
                    "CD_fromGust1toF": {"Value": CD_fromGust1toF, "Unit": "Non dimensional"},
                    "q_fromFtoE"    : {"Value": q_fromFtoE, "Unit": "Pa"},
                    "CD_fromFtoE"    : {"Value": CD_fromFtoE, "Unit": "Non dimensional"},
                    "q_fromEto0"    : {"Value": q_fromEto0, "Unit": "Pa"},
                    "CD_fromEto0"    : {"Value": CD_fromEto0, "Unit": "Non dimensional"},

                    "alfa_from0toS": {"Value": alfa_from0toS, "Unit": "deg"},
                    "alfa_new_from0toS": {"Value": alfa_new_from0toS, "Unit": "deg"},
                    "alfa_fromStoA": {"Value": alfa_fromStoA, "Unit": "deg"},
                    "alfa_new_fromStoA": {"Value": alfa_new_fromStoA, "Unit": "deg"},
                    "alfa_fromAtoGust1": {"Value": alfa_fromAtoGust1, "Unit": "deg"},
                    "alfa_new_fromAtoGust1": {"Value": alfa_new_fromAtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toC": {"Value": alfa_fromGust1toC, "Unit": "deg"},
                    "alfa_new_fromGust1toC": {"Value": alfa_new_fromGust1toC, "Unit": "deg"},
                    "alfa_fromCtoGust2": {"Value": alfa_fromCtoGust2, "Unit": "deg"},
                    "alfa_new_fromCtoGust2": {"Value": alfa_new_fromCtoGust2, "Unit": "deg"},
                    "alfa_fromGust2toD": {"Value": alfa_fromGust2toD, "Unit": "deg"},
                    "alfa_new_fromGust2toD": {"Value": alfa_new_fromGust2toD, "Unit": "deg"},
                    "alfa_fromDto0":     {"Value": alfa_fromDto0, "Unit": "deg"},
                    "alfa_new_fromDto0":     {"Value": alfa_new_fromDto0, "Unit": "deg"},
                    "alfa_from0toSinv":  {"Value": alfa_from0toSinv, "Unit": "deg"},
                    "alfa_new_from0toSinv":  {"Value": alfa_new_from0toSinv, "Unit": "deg"},
                    "alfa_fromSinvtoG":  {"Value": alfa_fromSinvtoG, "Unit": "deg"},
                    "alfa_new_fromSinvtoG":  {"Value": alfa_new_fromSinvtoG, "Unit": "deg"},
                    "alfa_fromGtoGust1": {"Value": alfa_fromGtoGust1, "Unit": "deg"},
                    "alfa_new_fromGtoGust1": {"Value": alfa_new_fromGtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toF": {"Value": alfa_fromGust1toF, "Unit": "deg"},
                    "alfa_new_fromGust1toF": {"Value": alfa_new_fromGust1toF, "Unit": "deg"},
                    "alfa_fromFtoE"    : {"Value": alfa_fromFtoE, "Unit": "deg"},
                    "alfa_new_fromFtoE"    : {"Value": alfa_new_fromFtoE, "Unit": "deg"},
                    "alfa_fromEto0"    : {"Value": alfa_fromEto0, "Unit": "deg"},
                    "alfa_new_fromEto0"    : {"Value": alfa_new_fromEto0, "Unit": "deg"},

                    "CL_wb_from0toS": {"Value": CL_wb_from0toS, "Unit": "Non dimensional"},
                    "CL_ht_from0toS": {"Value": CL_ht_from0toS, "Unit": "Non dimensional"},
                    "CL_new_from0toS": {"Value": CL_new_from0toS, "Unit": "Non dimensional"},
                    "CL_wb_fromStoA": {"Value": CL_wb_fromStoA, "Unit": "Non dimensional"},
                    "CL_ht_fromStoA": {"Value": CL_ht_fromStoA, "Unit": "Non dimensional"},
                    "CL_new_fromStoA": {"Value": CL_new_fromStoA, "Unit": "Non dimensional"},
                    "CL_wb_fromAtoGust1": {"Value": CL_wb_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_ht_fromAtoGust1": {"Value": CL_ht_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_new_fromAtoGust1": {"Value": CL_new_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_wb_fromGust1toC": {"Value": CL_wb_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toC": {"Value": CL_ht_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_new_fromGust1toC": {"Value": CL_new_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_wb_fromCtoGust2": {"Value": CL_wb_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_ht_fromCtoGust2": {"Value": CL_ht_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_new_fromCtoGust2": {"Value": CL_new_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_wb_fromGust2toD": {"Value": CL_wb_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_ht_fromGust2toD": {"Value": CL_ht_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_new_fromGust2toD": {"Value": CL_new_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_wb_fromDto0":     {"Value": CL_wb_fromDto0, "Unit": "Non dimensional"},
                    "CL_ht_fromDto0":     {"Value": CL_ht_fromDto0, "Unit": "Non dimensional"},
                    "CL_new_fromDto0":     {"Value": CL_new_fromDto0, "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv":  {"Value": CL_wb_from0toSinv, "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv":  {"Value": CL_ht_from0toSinv, "Unit": "Non dimensional"},
                    "CL_new_from0toSinv":  {"Value": CL_new_from0toSinv, "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG":  {"Value": CL_wb_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG":  {"Value": CL_ht_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG":  {"Value": CL_new_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_wb_fromGtoGust1": {"Value": CL_wb_fromGtoGust1, "Unit": "Non dimensional"},
                    "CL_ht_fromGtoGust1": {"Value": CL_ht_fromGtoGust1, "Unit": "Non dimensional"},
                    "CL_new_fromGtoGust1": {"Value": CL_new_fromGtoGust1, "Unit": "Non dimensional"},
                    "CL_wb_fromGust1toF": {"Value": CL_wb_fromGust1toF, "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toF": {"Value": CL_ht_fromGust1toF, "Unit": "Non dimensional"},
                    "CL_new_fromGust1toF": {"Value": CL_new_fromGust1toF, "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE"    : {"Value": CL_wb_fromFtoE, "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE"    : {"Value": CL_ht_fromFtoE, "Unit": "Non dimensional"},
                    "CL_new_fromFtoE"    : {"Value": CL_new_fromFtoE, "Unit": "Non dimensional"},
                    "CL_wb_fromEto0"    : {"Value": CL_wb_fromEto0, "Unit": "Non dimensional"},
                    "CL_ht_fromEto0"    : {"Value": CL_ht_fromEto0, "Unit": "Non dimensional"},
                    "CL_new_fromEto0"    : {"Value": CL_new_fromEto0, "Unit": "Non dimensional"},

                    "CM_due_to_CL_from0toS": {"Value": CM_due_to_CL_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS": {"Value": CM_due_to_CD_from0toS, "Unit": "Non dimensional"},
                    "CM_CG_from0toS": {"Value": CM_CG_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA": {"Value": CM_due_to_CL_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA": {"Value": CM_due_to_CD_fromStoA, "Unit": "Non dimensional"},
                    "CM_CG_fromStoA": {"Value": CM_CG_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoGust1": {"Value": CM_due_to_CL_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoGust1": {"Value": CM_due_to_CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_CG_fromAtoGust1": {"Value": CM_CG_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toC": {"Value": CM_due_to_CL_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toC": {"Value": CM_due_to_CD_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toC": {"Value": CM_CG_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoGust2": {"Value": CM_due_to_CL_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoGust2": {"Value": CM_due_to_CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_CG_fromCtoGust2": {"Value": CM_CG_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust2toD": {"Value": CM_due_to_CL_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust2toD": {"Value": CM_due_to_CD_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_CG_fromGust2toD": {"Value": CM_CG_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0":     {"Value": CM_due_to_CL_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0":     {"Value": CM_due_to_CD_fromDto0, "Unit": "Non dimensional"},
                    "CM_CG_fromDto0":     {"Value": CM_CG_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv":  {"Value": CM_due_to_CL_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv":  {"Value": CM_due_to_CD_from0toSinv, "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv":  {"Value": CM_CG_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG":  {"Value": CM_due_to_CL_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG":  {"Value": CM_due_to_CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG":  {"Value": CM_CG_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoGust1": {"Value": CM_due_to_CL_fromGtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoGust1": {"Value": CM_due_to_CD_fromGtoGust1, "Unit": "Non dimensional"},
                    "CM_CG_fromGtoGust1": {"Value": CM_CG_fromGtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toF": {"Value": CM_due_to_CL_fromGust1toF, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toF": {"Value": CM_due_to_CD_fromGust1toF, "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toF": {"Value": CM_CG_fromGust1toF, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE"    : {"Value": CM_due_to_CL_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE"    : {"Value": CM_due_to_CD_fromFtoE, "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE"    : {"Value": CM_CG_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0"    : {"Value": CM_due_to_CL_fromEto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0"    : {"Value": CM_due_to_CD_fromEto0, "Unit": "Non dimensional"},
                    "CM_CG_fromEto0"    : {"Value": CM_CG_fromEto0, "Unit": "Non dimensional"},

                    "L_wb_from0toS": {"Value": L_wb_from0toS, "Unit": "daN"},
                    "L_ht_from0toS": {"Value": L_ht_from0toS, "Unit": "daN"},
                    "L_new_from0toS": {"Value": L_new_from0toS, "Unit": "daN"},
                    "L_wb_fromStoA": {"Value": L_wb_fromStoA, "Unit": "daN"},
                    "L_ht_fromStoA": {"Value": L_ht_fromStoA, "Unit": "daN"},
                    "L_new_fromStoA": {"Value": L_new_fromStoA, "Unit": "daN"},
                    "L_wb_fromAtoGust1": {"Value": L_wb_fromAtoGust1, "Unit": "daN"},
                    "L_ht_fromAtoGust1": {"Value": L_ht_fromAtoGust1, "Unit": "daN"},
                    "L_new_fromAtoGust1": {"Value": L_new_fromAtoGust1, "Unit": "daN"},
                    "L_wb_fromGust1toC": {"Value": L_wb_fromGust1toC, "Unit": "daN"},
                    "L_ht_fromGust1toC": {"Value": L_ht_fromGust1toC, "Unit": "daN"},
                    "L_new_fromGust1toC": {"Value": L_new_fromGust1toC, "Unit": "daN"},
                    "L_wb_fromCtoGust2": {"Value": L_wb_fromCtoGust2, "Unit": "daN"},
                    "L_ht_fromCtoGust2": {"Value": L_ht_fromCtoGust2, "Unit": "daN"},
                    "L_new_fromCtoGust2": {"Value": L_new_fromCtoGust2, "Unit": "daN"},
                    "L_wb_fromGust2toD": {"Value": L_wb_fromGust2toD, "Unit": "daN"},
                    "L_ht_fromGust2toD": {"Value": L_ht_fromGust2toD, "Unit": "daN"},
                    "L_new_fromGust2toD": {"Value": L_new_fromGust2toD, "Unit": "daN"},
                    "L_wb_fromDto0":     {"Value": L_wb_fromDto0, "Unit": "daN"},
                    "L_ht_fromDto0":     {"Value": L_ht_fromDto0, "Unit": "daN"},
                    "L_new_fromDto0":     {"Value": L_new_fromDto0, "Unit": "daN"},
                    "L_wb_from0toSinv":  {"Value": L_wb_from0toSinv, "Unit": "daN"},
                    "L_ht_from0toSinv":  {"Value": L_ht_from0toSinv, "Unit": "daN"},
                    "L_new_from0toSinv":  {"Value": L_new_from0toSinv, "Unit": "daN"},
                    "L_wb_fromSinvtoG":  {"Value": L_wb_fromSinvtoG, "Unit": "daN"},
                    "L_ht_fromSinvtoG":  {"Value": L_ht_fromSinvtoG, "Unit": "daN"},
                    "L_new_fromSinvtoG":  {"Value": L_new_fromSinvtoG, "Unit": "daN"},
                    "L_wb_fromGtoGust1": {"Value": L_wb_fromGtoGust1, "Unit": "daN"},
                    "L_ht_fromGtoGust1": {"Value": L_ht_fromGtoGust1, "Unit": "daN"},
                    "L_new_fromGtoGust1": {"Value": L_new_fromGtoGust1, "Unit": "daN"},
                    "L_wb_fromGust1toF": {"Value": L_wb_fromGust1toF, "Unit": "daN"},
                    "L_ht_fromGust1toF": {"Value": L_ht_fromGust1toF, "Unit": "daN"},
                    "L_new_fromGust1toF": {"Value": L_new_fromGust1toF, "Unit": "daN"},
                    "L_wb_fromFtoE"    : {"Value": L_wb_fromFtoE, "Unit": "daN"},
                    "L_ht_fromFtoE"    : {"Value": L_ht_fromFtoE, "Unit": "daN"},
                    "L_new_fromFtoE"    : {"Value": L_new_fromFtoE, "Unit": "daN"},
                    "L_wb_fromEto0"    : {"Value": L_wb_fromEto0, "Unit": "daN"},
                    "L_ht_fromEto0"    : {"Value": L_ht_fromEto0, "Unit": "daN"},
                    "L_new_fromEto0"    : {"Value": L_new_fromEto0, "Unit": "daN"},
                    }
                }

        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case2_inverted"):
            
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS": {"Value": V_fe_from0toS, "Unit": "m/s"},
                    "n_fe_from0toS": {"Value": n_fe_from0toS, "Unit": "g"},
                    "V_fe_fromStoA": {"Value": V_fe_fromStoA, "Unit": "m/s"},
                    "n_fe_fromStoA": {"Value": n_fe_fromStoA, "Unit": "g"},
                    "V_fe_fromAtoGust1": {"Value": V_fe_fromAtoGust1, "Unit": "m/s"},
                    "n_fe_fromAtoGust1": {"Value": n_fe_fromAtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toC": {"Value": V_fe_fromGust1toC, "Unit": "m/s"},
                    "n_fe_fromGust1toC": {"Value": n_fe_fromGust1toC, "Unit": "g"},
                    "V_fe_fromCtoGust2": {"Value": V_fe_fromCtoGust2, "Unit": "m/s"},
                    "n_fe_fromCtoGust2": {"Value": n_fe_fromCtoGust2, "Unit": "g"},
                    "V_fe_fromGust2toD": {"Value": V_fe_fromGust2toD, "Unit": "m/s"},
                    "n_fe_fromGust2toD": {"Value": n_fe_fromGust2toD, "Unit": "g"},
                    "V_fe_fromDto0":     {"Value": V_fe_fromDto0, "Unit": "m/s"},
                    "n_fe_fromDto0":     {"Value": n_fe_fromDto0, "Unit": "g"},
                    "V_fe_from0toSinv":  {"Value": V_fe_from0toSinv, "Unit": "m/s"},
                    "n_fe_from0toSinv":  {"Value": n_fe_from0toSinv, "Unit": "g"},
                    "V_fe_fromSinvtoG":  {"Value": V_fe_fromSinvtoG, "Unit": "m/s"},
                    "n_fe_fromSinvtoG":  {"Value": n_fe_fromSinvtoG, "Unit": "g"},
                    "V_fe_fromGtoF":     {"Value": V_fe_fromGtoF, "Unit": "m/s"},
                    "n_fe_fromGtoF":     {"Value": n_fe_fromGtoF, "Unit": "g"},
                    "V_fe_fromFtoE":     {"Value": V_fe_fromFtoE, "Unit": "m/s"},
                    "n_fe_fromFtoE":     {"Value": n_fe_fromFtoE, "Unit": "g"},
                    "V_fe_fromEto0":     {"Value": V_fe_fromEto0, "Unit": "m/s"},
                    "n_fe_fromEto0":     {"Value": n_fe_fromEto0, "Unit": "g"},

                    "q_from0toS": {"Value": q_from0toS, "Unit": "Pa"},
                    "CD_from0toS": {"Value": CD_from0toS, "Unit": "Non dimensional"},
                    "q_fromStoA": {"Value": q_fromStoA, "Unit": "Pa"},
                    "CD_fromStoA": {"Value": CD_fromStoA, "Unit": "Non dimensional"},
                    "q_fromAtoGust1": {"Value": q_fromAtoGust1, "Unit": "Pa"},
                    "CD_fromAtoGust1": {"Value": CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toC": {"Value": q_fromGust1toC, "Unit": "Pa"},
                    "CD_fromGust1toC": {"Value": CD_fromGust1toC, "Unit": "Non dimensional"},
                    "q_fromCtoGust2": {"Value": q_fromCtoGust2, "Unit": "Pa"},
                    "CD_fromCtoGust2": {"Value": CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "q_fromGust2toD": {"Value": q_fromGust2toD, "Unit": "Pa"},
                    "CD_fromGust2toD": {"Value": CD_fromGust2toD, "Unit": "Non dimensional"},
                    "q_fromDto0":     {"Value": q_fromDto0, "Unit": "Pa"},
                    "CD_fromDto0":     {"Value": CD_fromDto0, "Unit": "Non dimensional"},
                    "q_from0toSinv":  {"Value": q_from0toSinv, "Unit": "Pa"},
                    "CD_from0toSinv":  {"Value": CD_from0toSinv, "Unit": "Non dimensional"},
                    "q_fromSinvtoG":  {"Value": q_fromSinvtoG, "Unit": "Pa"},
                    "CD_fromSinvtoG":  {"Value": CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "q_fromGtoF":     {"Value": q_fromGtoF, "Unit": "Pa"},
                    "CD_fromGtoF":     {"Value": CD_fromGtoF, "Unit": "Non dimensional"},
                    "q_fromFtoE":     {"Value": q_fromFtoE, "Unit": "Pa"},
                    "CD_fromFtoE":     {"Value": CD_fromFtoE, "Unit": "Non dimensional"},
                    "q_fromEto0":     {"Value": q_fromEto0, "Unit": "Pa"},
                    "CD_fromEto0":     {"Value": CD_fromEto0, "Unit": "Non dimensional"},

                    "alfa_from0toS": {"Value": alfa_from0toS, "Unit": "deg"},
                    "alfa_new_from0toS": {"Value": alfa_new_from0toS, "Unit": "deg"},
                    "alfa_fromStoA": {"Value": alfa_fromStoA, "Unit": "deg"},
                    "alfa_new_fromStoA": {"Value": alfa_new_fromStoA, "Unit": "deg"},
                    "alfa_fromAtoGust1": {"Value": alfa_fromAtoGust1, "Unit": "deg"},
                    "alfa_new_fromAtoGust1": {"Value": alfa_new_fromAtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toC": {"Value": alfa_fromGust1toC, "Unit": "deg"},
                    "alfa_new_fromGust1toC": {"Value": alfa_new_fromGust1toC, "Unit": "deg"},
                    "alfa_fromCtoGust2": {"Value": alfa_fromCtoGust2, "Unit": "deg"},
                    "alfa_new_fromCtoGust2": {"Value": alfa_new_fromCtoGust2, "Unit": "deg"},
                    "alfa_fromGust2toD": {"Value": alfa_fromGust2toD, "Unit": "deg"},
                    "alfa_new_fromGust2toD": {"Value": alfa_new_fromGust2toD, "Unit": "deg"},
                    "alfa_fromDto0":     {"Value": alfa_fromDto0, "Unit": "deg"},
                    "alfa_new_fromDto0":     {"Value": alfa_new_fromDto0, "Unit": "deg"},
                    "alfa_from0toSinv":  {"Value": alfa_from0toSinv, "Unit": "deg"},
                    "alfa_new_from0toSinv":  {"Value": alfa_new_from0toSinv, "Unit": "deg"},
                    "alfa_fromSinvtoG":  {"Value": alfa_fromSinvtoG, "Unit": "deg"},
                    "alfa_new_fromSinvtoG":  {"Value": alfa_new_fromSinvtoG, "Unit": "deg"},
                    "alfa_fromGtoF":     {"Value": alfa_fromGtoF, "Unit": "deg"},
                    "alfa_new_fromGtoF":     {"Value": alfa_new_fromGtoF, "Unit": "deg"},
                    "alfa_fromFtoE":     {"Value": alfa_fromFtoE, "Unit": "deg"},
                    "alfa_new_fromFtoE":     {"Value": alfa_new_fromFtoE, "Unit": "deg"},
                    "alfa_fromEto0":     {"Value": alfa_fromEto0, "Unit": "deg"},
                    "alfa_new_fromEto0":     {"Value": alfa_new_fromEto0, "Unit": "deg"},

                    "CL_wb_from0toS": {"Value": CL_wb_from0toS, "Unit": "Non dimensional"},
                    "CL_ht_from0toS": {"Value": CL_ht_from0toS, "Unit": "Non dimensional"},
                    "CL_new_from0toS": {"Value": CL_new_from0toS, "Unit": "Non dimensional"},
                    "CL_wb_fromStoA": {"Value": CL_wb_fromStoA, "Unit": "Non dimensional"},
                    "CL_ht_fromStoA": {"Value": CL_ht_fromStoA, "Unit": "Non dimensional"},
                    "CL_new_fromStoA": {"Value": CL_new_fromStoA, "Unit": "Non dimensional"},
                    "CL_wb_fromAtoGust1": {"Value": CL_wb_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_ht_fromAtoGust1": {"Value": CL_ht_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_new_fromAtoGust1": {"Value": CL_new_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_wb_fromGust1toC": {"Value": CL_wb_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toC": {"Value": CL_ht_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_new_fromGust1toC": {"Value": CL_new_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_wb_fromCtoGust2": {"Value": CL_wb_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_ht_fromCtoGust2": {"Value": CL_ht_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_new_fromCtoGust2": {"Value": CL_new_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_wb_fromGust2toD": {"Value": CL_wb_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_ht_fromGust2toD": {"Value": CL_ht_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_new_fromGust2toD": {"Value": CL_new_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_wb_fromDto0":     {"Value": CL_wb_fromDto0, "Unit": "Non dimensional"},
                    "CL_ht_fromDto0":     {"Value": CL_ht_fromDto0, "Unit": "Non dimensional"},
                    "CL_new_fromDto0":     {"Value": CL_new_fromDto0, "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv":  {"Value": CL_wb_from0toSinv, "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv":  {"Value": CL_ht_from0toSinv, "Unit": "Non dimensional"},
                    "CL_new_from0toSinv":  {"Value": CL_new_from0toSinv, "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG":  {"Value": CL_wb_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG":  {"Value": CL_ht_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG":  {"Value": CL_new_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_wb_fromGtoF":     {"Value": CL_wb_fromGtoF, "Unit": "Non dimensional"},
                    "CL_ht_fromGtoF":     {"Value": CL_ht_fromGtoF, "Unit": "Non dimensional"},
                    "CL_new_fromGtoF":     {"Value": CL_new_fromGtoF, "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE":     {"Value": CL_wb_fromFtoE, "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE":     {"Value": CL_ht_fromFtoE, "Unit": "Non dimensional"},
                    "CL_new_fromFtoE":     {"Value": CL_new_fromFtoE, "Unit": "Non dimensional"},
                    "CL_wb_fromEto0":     {"Value": CL_wb_fromEto0, "Unit": "Non dimensional"},
                    "CL_ht_fromEto0":     {"Value": CL_ht_fromEto0, "Unit": "Non dimensional"},
                    "CL_new_fromEto0":     {"Value": CL_new_fromEto0, "Unit": "Non dimensional"},

                    "CM_due_to_CL_from0toS": {"Value": CM_due_to_CL_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS": {"Value": CM_due_to_CD_from0toS, "Unit": "Non dimensional"},
                    "CM_CG_from0toS": {"Value": CM_CG_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA": {"Value": CM_due_to_CL_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA": {"Value": CM_due_to_CD_fromStoA, "Unit": "Non dimensional"},
                    "CM_CG_fromStoA": {"Value": CM_CG_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoGust1": {"Value": CM_due_to_CL_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoGust1": {"Value": CM_due_to_CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_CG_fromAtoGust1": {"Value": CM_CG_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toC": {"Value": CM_due_to_CL_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toC": {"Value": CM_due_to_CD_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toC": {"Value": CM_CG_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoGust2": {"Value": CM_due_to_CL_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoGust2": {"Value": CM_due_to_CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_CG_fromCtoGust2": {"Value": CM_CG_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust2toD": {"Value": CM_due_to_CL_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust2toD": {"Value": CM_due_to_CD_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_CG_fromGust2toD": {"Value": CM_CG_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0":     {"Value": CM_due_to_CL_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0":     {"Value": CM_due_to_CD_fromDto0, "Unit": "Non dimensional"},
                    "CM_CG_fromDto0":     {"Value": CM_CG_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv":  {"Value": CM_due_to_CL_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv":  {"Value": CM_due_to_CD_from0toSinv, "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv":  {"Value": CM_CG_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG":  {"Value": CM_due_to_CL_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG":  {"Value": CM_due_to_CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG":  {"Value": CM_CG_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoF":     {"Value": CM_due_to_CL_fromGtoF, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoF":     {"Value": CM_due_to_CD_fromGtoF, "Unit": "Non dimensional"},
                    "CM_CG_fromGtoF":     {"Value": CM_CG_fromGtoF, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE":     {"Value": CM_due_to_CL_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE":     {"Value": CM_due_to_CD_fromFtoE, "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE":     {"Value": CM_CG_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0":     {"Value": CM_due_to_CL_fromEto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0":     {"Value": CM_due_to_CD_fromEto0, "Unit": "Non dimensional"},
                    "CM_CG_fromEto0":     {"Value": CM_CG_fromEto0, "Unit": "Non dimensional"},
                    }
                }

        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case3_inverted"):
            
            
            Balancing_loads = {
                "Balancing_loads": {
                    "V_fe_from0toS": {"Value": V_fe_from0toS, "Unit": "m/s"},
                    "n_fe_from0toS": {"Value": n_fe_from0toS, "Unit": "g"},
                    "V_fe_fromStoA": {"Value": V_fe_fromStoA, "Unit": "m/s"},
                    "n_fe_fromStoA": {"Value": n_fe_fromStoA, "Unit": "g"},
                    "V_fe_fromAtoGust1": {"Value": V_fe_fromAtoGust1, "Unit": "m/s"},
                    "n_fe_fromAtoGust1": {"Value": n_fe_fromAtoGust1, "Unit": "g"},
                    "V_fe_fromGust1toC": {"Value": V_fe_fromGust1toC, "Unit": "m/s"},
                    "n_fe_fromGust1toC": {"Value": n_fe_fromGust1toC, "Unit": "g"},
                    "V_fe_fromCtoGust2": {"Value": V_fe_fromCtoGust2, "Unit": "m/s"},
                    "n_fe_fromCtoGust2": {"Value": n_fe_fromCtoGust2, "Unit": "g"},
                    "V_fe_fromGust2toD": {"Value": V_fe_fromGust2toD, "Unit": "m/s"},
                    "n_fe_fromGust2toD": {"Value": n_fe_fromGust2toD, "Unit": "g"},
                    "V_fe_fromDto0":     {"Value": V_fe_fromDto0, "Unit": "m/s"},
                    "n_fe_fromDto0":     {"Value": n_fe_fromDto0, "Unit": "g"},
                    "V_fe_from0toSinv":  {"Value": V_fe_from0toSinv, "Unit": "m/s"},
                    "n_fe_from0toSinv":  {"Value": n_fe_from0toSinv, "Unit": "g"},
                    "V_fe_fromSinvtoG":  {"Value": V_fe_fromSinvtoG, "Unit": "m/s"},
                    "n_fe_fromSinvtoG":  {"Value": n_fe_fromSinvtoG, "Unit": "g"},
                    "V_fe_fromGtoF":     {"Value": V_fe_fromGtoF, "Unit": "m/s"},
                    "n_fe_fromGtoF":     {"Value": n_fe_fromGtoF, "Unit": "g"},
                    "V_fe_fromFtoE"    : {"Value": V_fe_fromFtoE, "Unit": "m/s"},
                    "n_fe_fromFtoE"    : {"Value": n_fe_fromFtoE, "Unit": "g"},
                    "V_fe_fromEto0"    : {"Value": V_fe_fromEto0, "Unit": "m/s"},
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"},

                    "q_from0toS": {"Value": q_from0toS, "Unit": "Pa"},
                    "CD_from0toS": {"Value": CD_from0toS, "Unit": "Non dimensional"},
                    "q_fromStoA": {"Value": q_fromStoA, "Unit": "Pa"},
                    "CD_fromStoA": {"Value": CD_fromStoA, "Unit": "Non dimensional"},
                    "q_fromAtoGust1": {"Value": q_fromAtoGust1, "Unit": "Pa"},
                    "CD_fromAtoGust1": {"Value": CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "q_fromGust1toC": {"Value": q_fromGust1toC, "Unit": "Pa"},
                    "CD_fromGust1toC": {"Value": CD_fromGust1toC, "Unit": "Non dimensional"},
                    "q_fromCtoGust2": {"Value": q_fromCtoGust2, "Unit": "Pa"},
                    "CD_fromCtoGust2": {"Value": CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "q_fromGust2toD": {"Value": q_fromGust2toD, "Unit": "Pa"},
                    "CD_fromGust2toD": {"Value": CD_fromGust2toD, "Unit": "Non dimensional"},
                    "q_fromDto0":     {"Value": q_fromDto0, "Unit": "Pa"},
                    "CD_fromDto0":     {"Value": CD_fromDto0, "Unit": "Non dimensional"},
                    "q_from0toSinv":  {"Value": q_from0toSinv, "Unit": "Pa"},
                    "CD_from0toSinv":  {"Value": CD_from0toSinv, "Unit": "Non dimensional"},
                    "q_fromSinvtoG":  {"Value": q_fromSinvtoG, "Unit": "Pa"},
                    "CD_fromSinvtoG":  {"Value": CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "q_fromGtoF":     {"Value": q_fromGtoF, "Unit": "Pa"},
                    "CD_fromGtoF":     {"Value": CD_fromGtoF, "Unit": "Non dimensional"},
                    "q_fromFtoE"    : {"Value": q_fromFtoE, "Unit": "Pa"},
                    "CD_fromFtoE"    : {"Value": CD_fromFtoE, "Unit": "Non dimensional"},
                    "q_fromEto0"    : {"Value": q_fromEto0, "Unit": "Pa"},
                    "CD_fromEto0"    : {"Value": CD_fromEto0, "Unit": "Non dimensional"},

                    "alfa_from0toS": {"Value": alfa_from0toS, "Unit": "deg"},
                    "alfa_new_from0toS": {"Value": alfa_new_from0toS, "Unit": "deg"},
                    "alfa_fromStoA": {"Value": alfa_fromStoA, "Unit": "deg"},
                    "alfa_new_fromStoA": {"Value": alfa_new_fromStoA, "Unit": "deg"},
                    "alfa_fromAtoGust1": {"Value": alfa_fromAtoGust1, "Unit": "deg"},
                    "alfa_new_fromAtoGust1": {"Value": alfa_new_fromAtoGust1, "Unit": "deg"},
                    "alfa_fromGust1toC": {"Value": alfa_fromGust1toC, "Unit": "deg"},
                    "alfa_new_fromGust1toC": {"Value": alfa_new_fromGust1toC, "Unit": "deg"},
                    "alfa_fromCtoGust2": {"Value": alfa_fromCtoGust2, "Unit": "deg"},
                    "alfa_new_fromCtoGust2": {"Value": alfa_new_fromCtoGust2, "Unit": "deg"},
                    "alfa_fromGust2toD": {"Value": alfa_fromGust2toD, "Unit": "deg"},
                    "alfa_new_fromGust2toD": {"Value": alfa_new_fromGust2toD, "Unit": "deg"},
                    "alfa_fromDto0":     {"Value": alfa_fromDto0, "Unit": "deg"},
                    "alfa_new_fromDto0":     {"Value": alfa_new_fromDto0, "Unit": "deg"},
                    "alfa_from0toSinv":  {"Value": alfa_from0toSinv, "Unit": "deg"},
                    "alfa_new_from0toSinv":  {"Value": alfa_new_from0toSinv, "Unit": "deg"},
                    "alfa_fromSinvtoG":  {"Value": alfa_fromSinvtoG, "Unit": "deg"},
                    "alfa_new_fromSinvtoG":  {"Value": alfa_new_fromSinvtoG, "Unit": "deg"},
                    "alfa_fromGtoF":     {"Value": alfa_fromGtoF, "Unit": "deg"},
                    "alfa_new_fromGtoF":     {"Value": alfa_new_fromGtoF, "Unit": "deg"},
                    "alfa_fromFtoE"    : {"Value": alfa_fromFtoE, "Unit": "deg"},
                    "alfa_new_fromFtoE"    : {"Value": alfa_new_fromFtoE, "Unit": "deg"},
                    "alfa_fromEto0"    : {"Value": alfa_fromEto0, "Unit": "deg"},
                    "alfa_new_fromEto0"    : {"Value": alfa_new_fromEto0, "Unit": "deg"},

                    "CL_wb_from0toS": {"Value": CL_wb_from0toS, "Unit": "Non dimensional"},
                    "CL_ht_from0toS": {"Value": CL_ht_from0toS, "Unit": "Non dimensional"},
                    "CL_new_from0toS": {"Value": CL_new_from0toS, "Unit": "Non dimensional"},
                    "CL_wb_fromStoA": {"Value": CL_wb_fromStoA, "Unit": "Non dimensional"},
                    "CL_ht_fromStoA": {"Value": CL_ht_fromStoA, "Unit": "Non dimensional"},
                    "CL_new_fromStoA": {"Value": CL_new_fromStoA, "Unit": "Non dimensional"},
                    "CL_wb_fromAtoGust1": {"Value": CL_wb_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_ht_fromAtoGust1": {"Value": CL_ht_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_new_fromAtoGust1": {"Value": CL_new_fromAtoGust1, "Unit": "Non dimensional"},
                    "CL_wb_fromGust1toC": {"Value": CL_wb_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_ht_fromGust1toC": {"Value": CL_ht_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_new_fromGust1toC": {"Value": CL_new_fromGust1toC, "Unit": "Non dimensional"},
                    "CL_wb_fromCtoGust2": {"Value": CL_wb_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_ht_fromCtoGust2": {"Value": CL_ht_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_new_fromCtoGust2": {"Value": CL_new_fromCtoGust2, "Unit": "Non dimensional"},
                    "CL_wb_fromGust2toD": {"Value": CL_wb_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_ht_fromGust2toD": {"Value": CL_ht_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_new_fromGust2toD": {"Value": CL_new_fromGust2toD, "Unit": "Non dimensional"},
                    "CL_wb_fromDto0":     {"Value": CL_wb_fromDto0, "Unit": "Non dimensional"},
                    "CL_ht_fromDto0":     {"Value": CL_ht_fromDto0, "Unit": "Non dimensional"},
                    "CL_new_fromDto0":     {"Value": CL_new_fromDto0, "Unit": "Non dimensional"},
                    "CL_wb_from0toSinv":  {"Value": CL_wb_from0toSinv, "Unit": "Non dimensional"},
                    "CL_ht_from0toSinv":  {"Value": CL_ht_from0toSinv, "Unit": "Non dimensional"},
                    "CL_new_from0toSinv":  {"Value": CL_new_from0toSinv, "Unit": "Non dimensional"},
                    "CL_wb_fromSinvtoG":  {"Value": CL_wb_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_ht_fromSinvtoG":  {"Value": CL_ht_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_new_fromSinvtoG":  {"Value": CL_new_fromSinvtoG, "Unit": "Non dimensional"},
                    "CL_wb_fromGtoF":     {"Value": CL_wb_fromGtoF, "Unit": "Non dimensional"},
                    "CL_ht_fromGtoF":     {"Value": CL_ht_fromGtoF, "Unit": "Non dimensional"},
                    "CL_new_fromGtoF":     {"Value": CL_new_fromGtoF, "Unit": "Non dimensional"},
                    "CL_wb_fromFtoE"    : {"Value": CL_wb_fromFtoE, "Unit": "Non dimensional"},
                    "CL_ht_fromFtoE"    : {"Value": CL_ht_fromFtoE, "Unit": "Non dimensional"},
                    "CL_new_fromFtoE"    : {"Value": CL_new_fromFtoE, "Unit": "Non dimensional"},
                    "CL_wb_fromEto0"    : {"Value": CL_wb_fromEto0, "Unit": "Non dimensional"},
                    "CL_ht_fromEto0"    : {"Value": CL_ht_fromEto0, "Unit": "Non dimensional"},
                    "CL_new_fromEto0"    : {"Value": CL_new_fromEto0, "Unit": "Non dimensional"},

                    "L_wb_from0toS": {"Value": L_wb_from0toS, "Unit": "daN"},
                    "L_ht_from0toS": {"Value": L_ht_from0toS, "Unit": "daN"},
                    "L_new_from0toS": {"Value": L_new_from0toS, "Unit": "daN"},
                    "L_wb_fromStoA": {"Value": L_wb_fromStoA, "Unit": "daN"},
                    "L_ht_fromStoA": {"Value": L_ht_fromStoA, "Unit": "daN"},
                    "L_new_fromStoA": {"Value": L_new_fromStoA, "Unit": "daN"},
                    "L_wb_fromAtoGust1": {"Value": L_wb_fromAtoGust1, "Unit": "daN"},
                    "L_ht_fromAtoGust1": {"Value": L_ht_fromAtoGust1, "Unit": "daN"},
                    "L_new_fromAtoGust1": {"Value": L_new_fromAtoGust1, "Unit": "daN"},
                    "L_wb_fromGust1toC": {"Value": L_wb_fromGust1toC, "Unit": "daN"},
                    "L_ht_fromGust1toC": {"Value": L_ht_fromGust1toC, "Unit": "daN"},
                    "L_new_fromGust1toC": {"Value": L_new_fromGust1toC, "Unit": "daN"},
                    "L_wb_fromCtoGust2": {"Value": L_wb_fromCtoGust2, "Unit": "daN"},
                    "L_ht_fromCtoGust2": {"Value": L_ht_fromCtoGust2, "Unit": "daN"},
                    "L_new_fromCtoGust2": {"Value": L_new_fromCtoGust2, "Unit": "daN"},
                    "L_wb_fromGust2toD": {"Value": L_wb_fromGust2toD, "Unit": "daN"},
                    "L_ht_fromGust2toD": {"Value": L_ht_fromGust2toD, "Unit": "daN"},
                    "L_new_fromGust2toD": {"Value": L_new_fromGust2toD, "Unit": "daN"},
                    "L_wb_fromDto0":     {"Value": L_wb_fromDto0, "Unit": "daN"},
                    "L_ht_fromDto0":     {"Value": L_ht_fromDto0, "Unit": "daN"},
                    "L_new_fromDto0":     {"Value": L_new_fromDto0, "Unit": "daN"},
                    "L_wb_from0toSinv":  {"Value": L_wb_from0toSinv, "Unit": "daN"},
                    "L_ht_from0toSinv":  {"Value": L_ht_from0toSinv, "Unit": "daN"},
                    "L_new_from0toSinv":  {"Value": L_new_from0toSinv, "Unit": "daN"},
                    "L_wb_fromSinvtoG":  {"Value": L_wb_fromSinvtoG, "Unit": "daN"},
                    "L_ht_fromSinvtoG":  {"Value": L_ht_fromSinvtoG, "Unit": "daN"},
                    "L_new_fromSinvtoG":  {"Value": L_new_fromSinvtoG, "Unit": "daN"},
                    "L_wb_fromGtoF":     {"Value": L_wb_fromGtoF, "Unit": "daN"},
                    "L_ht_fromGtoF":     {"Value": L_ht_fromGtoF, "Unit": "daN"},
                    "L_new_fromGtoF":     {"Value": L_new_fromGtoF, "Unit": "daN"},
                    "L_wb_fromFtoE"    : {"Value": L_wb_fromFtoE, "Unit": "daN"},
                    "L_ht_fromFtoE"    : {"Value": L_ht_fromFtoE, "Unit": "daN"},
                    "L_new_fromFtoE"    : {"Value": L_new_fromFtoE, "Unit": "daN"},
                    "L_wb_fromEto0"    : {"Value": L_wb_fromEto0, "Unit": "daN"},
                    "L_ht_fromEto0"    : {"Value": L_ht_fromEto0, "Unit": "daN"},
                    "L_new_fromEto0"    : {"Value": L_new_fromEto0, "Unit": "daN"},

                    "CM_due_to_CL_from0toS": {"Value": CM_due_to_CL_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toS": {"Value": CM_due_to_CD_from0toS, "Unit": "Non dimensional"},
                    "CM_CG_from0toS": {"Value": CM_CG_from0toS, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromStoA": {"Value": CM_due_to_CL_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromStoA": {"Value": CM_due_to_CD_fromStoA, "Unit": "Non dimensional"},
                    "CM_CG_fromStoA": {"Value": CM_CG_fromStoA, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromAtoGust1": {"Value": CM_due_to_CL_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromAtoGust1": {"Value": CM_due_to_CD_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_CG_fromAtoGust1": {"Value": CM_CG_fromAtoGust1, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust1toC": {"Value": CM_due_to_CL_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust1toC": {"Value": CM_due_to_CD_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_CG_fromGust1toC": {"Value": CM_CG_fromGust1toC, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromCtoGust2": {"Value": CM_due_to_CL_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromCtoGust2": {"Value": CM_due_to_CD_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_CG_fromCtoGust2": {"Value": CM_CG_fromCtoGust2, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGust2toD": {"Value": CM_due_to_CL_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGust2toD": {"Value": CM_due_to_CD_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_CG_fromGust2toD": {"Value": CM_CG_fromGust2toD, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromDto0":     {"Value": CM_due_to_CL_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromDto0":     {"Value": CM_due_to_CD_fromDto0, "Unit": "Non dimensional"},
                    "CM_CG_fromDto0":     {"Value": CM_CG_fromDto0, "Unit": "Non dimensional"},
                    "CM_due_to_CL_from0toSinv":  {"Value": CM_due_to_CL_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CD_from0toSinv":  {"Value": CM_due_to_CD_from0toSinv, "Unit": "Non dimensional"},
                    "CM_CG_from0toSinv":  {"Value": CM_CG_from0toSinv, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromSinvtoG":  {"Value": CM_due_to_CL_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromSinvtoG":  {"Value": CM_due_to_CD_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_CG_fromSinvtoG":  {"Value": CM_CG_fromSinvtoG, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromGtoF":     {"Value": CM_due_to_CL_fromGtoF, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromGtoF":     {"Value": CM_due_to_CD_fromGtoF, "Unit": "Non dimensional"},
                    "CM_CG_fromGtoF":     {"Value": CM_CG_fromGtoF, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromFtoE"    : {"Value": CM_due_to_CL_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromFtoE"    : {"Value": CM_due_to_CD_fromFtoE, "Unit": "Non dimensional"},
                    "CM_CG_fromFtoE"    : {"Value": CM_CG_fromFtoE, "Unit": "Non dimensional"},
                    "CM_due_to_CL_fromEto0"    : {"Value": CM_due_to_CL_fromEto0, "Unit": "Non dimensional"},
                    "CM_due_to_CD_fromEto0"    : {"Value": CM_due_to_CD_fromEto0, "Unit": "Non dimensional"},
                    "CM_CG_fromEto0"    : {"Value": CM_CG_fromEto0, "Unit": "Non dimensional"},
                    }
                }
       
        return Balancing_loads                                                                      