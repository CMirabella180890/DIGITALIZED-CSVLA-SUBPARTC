# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 09:02:38 2022

@author: claum
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 18:37:28 2021

@author: claum
"""
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

class csvla: 
    """
    CS - VLA Amdt1 -  05 Mar 2009
    
      CS - VLA 1 Applicability
      This airworthiness code is applicable to aeroplanes with a single
      engine, spark or compression-ignition, having not more than two
      seats, with a Max Certificated Takeoff Weight of not more than 750 kg
      and a stalling speed in the landing configuration of not more than 83
      km/h (45 KCAS) to be approved for day-VFR only. 
      
      CS - VLA 3 Aeroplane categories
      This CS - VLA applies to aeroplanes intended for non - aerobatic
      operation only. Non - aerobatic operation includes 
      (a) Any manoeuvre incident to normal flying;
      (b) Stalls, excepts whip stalls; and 
      (c) Lazy eights, chandelles, and steep turns, in which the angle of
          bank is not more than 60 degrees. 
      CS - VLA 301 Loads 
      (a) Strength requrements are specified in terms of limit loads,
          the maximum loads to be expected in service, and ultimate
          loads, limit loads multiplied by prescribed factors of
          safety. Unless otherwise provided, prescribed loads are
          limit loads. 
      (b) Unless otherwise provided, the air, groun, and water loads
          must be placed in equilibrium with inertia forces,
          considering each item of mass in the aeroplane. These loads
          must be distributed to conservatively approximate or closely
          represent actual conditions. 
      (c) If deflections under load would significantly change the
          distribution of external or internal loads, this
          redistribution must be taken into account.
      (d) Simplified structural design criteria given in this Subpart
          C and its appendices may be used only for aeroplanes with
          conventional configurations. If Appendix A is used, the
          entire appendix must be substituted for the corresponding
          paragraphs of this subpart, i.e. CS - VLA 321 to 459, see
          also CS - VLA 301 (d).  
    """
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # LOAD FACTOR CALCULATOR
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def calcn(nmax, n_Elements):
        """
            n = calcn(cls, nmax)
             Function that calculates load factors values along the
             stall curve for flight envelope calculation.
            
             CS - VLA 337 Limit manoeuvring load factors
               (a) The positive limit manoeuvring load factor n may not be
                   less than 3.8. 
               (b) The negative limit manoeuvring load factor may not be
                   less than -1.5.
              
              INPUT
              nmax = Appliable limit load factor
              OUTPUT 
              n    = Vector of load factor values
        """
        if nmax > 0.0:
            n = np.linspace(1.0, nmax, n_Elements) 
        elif nmax < 0.0: 
            n = np.linspace(-1.0, nmax, n_Elements)
        return n
    # calcn = staticmethod(calcn)
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # STALL SPEED CALCULATOR
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def calcvs(rho, WS, CLmax, n):
        """
            VS = calcvs(rho, WS, CLmax, n)
            Stall speed values correspondent to limit load factors from 1g 
            to limit g's. 

            CS - VLA 335 Design Airspeeds (1)(i)
            VS is a computed stalling speed with flaps retracted at the de-
            sign weight, normally based on the maximum aeroplane normal 
            force coefficients, CNA. 
            CS - VLA 335 Design Airspeeds (1)(ii)
            n is the limit manoeuvring load factor used in design. 

            INPUT 
            rho   = Density at the selected altitude [kg/m**3]
            WS    = Wing loading in [Pa]
            CLmax = Applicable maximum lift coefficient [Non dim.]
            n     = Vector of limit load factor values [g's].
            OUTPUT 
            VS    = Vector of stall speed values. 
        """
        VS = np.sqrt( WS * (2/rho) * (1/CLmax) * n )
        
        return VS
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # STALL SPEED CALCULATOR
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def calcvd(vcmin, vc):
        """
            VD = calcvd(vcmin, vc)
            Design dive speed CS - VLA (b)(1)(2) VD 
            Vd must not be less than (1.25)*Vc; and, with Vc min, the 
            required minimum design cruising speed, Vd may not be less
            than (1.40)*Vc_min. 

            INPUT 
            vcimn = The required minimum design cruising speed, VD may 
                    not be less than (1.40)*Vc_min 
            vc    = Design cruise speed, which may not be less than 
                    (4.7)*sqrt(W/S)
                    where 
                    W/S = Wing loading in [Pa]
                    Vc  = Cruise speed in [m/s]
            OUTPUT 
            VD    = Design dive speed
        """
        vd1 = 1.4*vcmin 
        vd2 = 1.25*vc

        if vd1 > vd2: 
            VD = vd1
        elif vd2 > vd1: 
            VD = vd2

        return VD
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # STALL SPEED CALCULATOR
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def calcvc(WS, VH):
        """
            VC = calcvc(WS, VH)
            Design cruise speed CS - VLA (b)(1)(2) VC 
            Vc (in [m/s]) may not be less than (2.4)*sqrt(W/S), where 
            W/S = Wing loading in [Pa] 
            Vc  = Cruise speed in [m/s]

            Vc need not to be more than (0.9)*Vh at sea level, where 
            Vh  = Max continous power max horizontal speed in [m/s]

            INPUT 
            WS = Wing loading in [Pa] 
            vh = Max continous power max horizontal speed in [m/s]

            OUTPUT
            VD = Design dive speed 
        """
        VD = (2.4)*np.sqrt(WS)
        x  = (0.9)*VH 
        
        if VH != 0: 
            if VD > x:
                VD = x
            elif VD < x:
                VD = VD
        elif VH == 0: 
            VD = VD
            

        return VD
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # V - n DIAGRAM 
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def V_n_diagram(n_from0toS, n_fromStoA, n_fromAtoC, n_fromCtoD, n_fromDto0,\
                    V_from0toS, V_fromStoA, V_fromAtoC, V_fromCtoD, V_fromDto0,\
                    n_from0toSinv, n_fromSinvtoG, n_fromGtoF, n_fromFto0,\
                    V_from0toSinv, V_fromSinvtoG, V_fromGtoF, V_fromFto0,\
                    Reg, Aircraft_name, n,\
                    nS, VS, nA, VA, nC, VC, nD, VD,\
                    nSinv, VSinv, nG, VG, nF, VF):
        """
        

        Parameters
        ----------
        n_from0toS : TYPE
            DESCRIPTION.
        n_fromStoA : TYPE
            DESCRIPTION.
        n_fromAtoC : TYPE
            DESCRIPTION.
        n_fromCtoD : TYPE
            DESCRIPTION.
        n_fromDto0 : TYPE
            DESCRIPTION.
        V_from0toS : TYPE
            DESCRIPTION.
        V_fromStoA : TYPE
            DESCRIPTION.
        V_fromAtoC : TYPE
            DESCRIPTION.
        V_fromCtoD : TYPE
            DESCRIPTION.
        V_fromDto0 : TYPE
            DESCRIPTION.
        n_from0toSinv : TYPE
            DESCRIPTION.
        n_fromSinvtoG : TYPE
            DESCRIPTION.
        n_fromGtoF : TYPE
            DESCRIPTION.
        n_fromFto0 : TYPE
            DESCRIPTION.
        V_from0toSinv : TYPE
            DESCRIPTION.
        V_fromSinvtoG : TYPE
            DESCRIPTION.
        V_fromGtoF : TYPE
            DESCRIPTION.
        V_fromFto0 : TYPE
            DESCRIPTION.
        Reg : TYPE
            DESCRIPTION.
        Aircraft_name : TYPE
            DESCRIPTION.
        n : TYPE
            DESCRIPTION.
        nS : TYPE
            DESCRIPTION.
        VS : TYPE
            DESCRIPTION.
        nA : TYPE
            DESCRIPTION.
        VA : TYPE
            DESCRIPTION.
        nC : TYPE
            DESCRIPTION.
        VC : TYPE
            DESCRIPTION.
        nD : TYPE
            DESCRIPTION.
        VD : TYPE
            DESCRIPTION.
        nSinv : TYPE
            DESCRIPTION.
        VSinv : TYPE
            DESCRIPTION.
        nG : TYPE
            DESCRIPTION.
        VG : TYPE
            DESCRIPTION.
        nF : TYPE
            DESCRIPTION.
        VF : TYPE
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
        fig  = plt.figure()
        plt.plot(V_from0toS,    n_from0toS,    color="red", linewidth=1.0, linestyle = "solid", zorder=1)
        plt.plot(V_fromStoA,    n_fromStoA,    color="red", linewidth=1.0, linestyle = "solid", zorder=1)
        plt.plot(V_fromAtoC,    n_fromAtoC,    color="red", linewidth=1.0, linestyle = "solid", zorder=1)
        plt.plot(V_fromCtoD,    n_fromCtoD,    color="red", linewidth=1.0, linestyle = "solid", zorder=1)
        plt.plot(V_fromDto0,    n_fromDto0,    color="red", linewidth=1.0, linestyle = "solid", zorder=1)
        plt.plot(V_from0toSinv, n_from0toSinv, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
        plt.plot(V_fromSinvtoG, n_fromSinvtoG, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
        plt.plot(V_fromGtoF,    n_fromGtoF,    color="red", linewidth=1.0, linestyle = "solid", zorder=1)
        plt.plot(V_fromFto0,    n_fromFto0,    color="red", linewidth=1.0, linestyle = "solid", zorder=1)
          
        plt.scatter(VS, nS, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
        # fontsize or size
        # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
        # 'x-large', 'xx-large'}
        # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
        # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
        # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
        label=r"Point S"
        plt.text(VS, nS, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VA, nA, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
        label=r"Point A"
        plt.text(VA, nA, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VC, nC, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
        label=r"Point C"
        plt.text(VC, nC, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VD, nD, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
        label=r"Point D"
        plt.text(VD-2.0, nD, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VSinv, nSinv, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point S inverted", zorder=2)
        label=r"Point S inverted"
        plt.text(VSinv, nSinv, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VG, nG, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
        label=r"Point G"
        plt.text(VG, nG, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VF, nF, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
        label=r"Point F"
        plt.text(VF, nF, label, fontdict=None, fontsize='x-small', fontweight='bold')
        
        plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
        plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
        # plt.ylim(nmin - 0.5, nmax + 0.5)
        # plt.xlim(0.0, VD + 10.0)
        # plt.legend(loc="best")
        plt.title(r'V - n diagram per ' + Reg + ' - ' + Aircraft_name) 
        plt.grid(True, linestyle='-.', which="both")
        plt.minorticks_on()
        name_figure = 'FlightEnvelope' + str(n) + '.pdf'
        plt.savefig(name_figure, bbox_inches='tight')
        plt.show()
        
        return fig
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # V - n DIAGRAM 
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AEROPLANE MASS FACTOR FUNCTION 
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def calcmug(WS, MGC, a, rho, g):
        """
        Aeroplane mass ratio = mug = calcmug(WS, MGC, a, rho, g)
        Function that calculates the aeroplane mass ratio mug. See page 41 
        of EASA CS - VLA Easy Access. In particular, 
        
        CS - VLA 341 Gust load factors
        In the abscence of a more rational analysis, the gust load factors
        may be computed as follows: 
            
                      0.5 * rho0 * V * a * Kg * Ude
            n = 1.0 + -----------------------------
                               ( W / S )
            where 
                  0.88 * mu_g
            Kg = ------------- = gust alleviation factor 
                  5.3 + mu_g
                  
                       2 * (M / S)
            mu_g = ------------------- = aeroplane mass ratio 
                   rho * MGC * CL_alfa
                   
            U_de    = derived gust velocities referred to in CS - VLA 333(c) [m/s]
            rho0    = density of air at sea level [kg/m^3]
            rho     = density of air at operative altitude [kg/m^3]
            M/S     = wing loading [kg/m^2]
            MGC     = mean geometric chord [m]
            g       = acceleration due to gravity [m/s^2]
            V       = aeroplane equivalent speed [m/s]
            CL_alfa = slope of the aeroplane normal force coefficient curve CNA
                      in [1/rad] if the gust loads are applied to the wings and
                      horizontal tail surfaces simultaneously by a rational me-
                      thod. The wing lift curve slope CL_alfa in [1/rad] may be
                      used when the gust load is applied to the wings only and 
                      the horizontal tail gust loads are treated as a separate 
                      condition.

        Returns
        -------
        mu_g : float 
            Aeroplane mass ratio.

        """
        x    = 2 * (WS / g)
        y    = rho * MGC * a 
        mu_g = x / y 
        
        return mu_g
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AEROPLANE GUST ALLEVIATION FACTOR 
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def calc_kg(mu_g):
        """
        Aeroplane mass ratio = mug = calcmug(WS, MGC, a, rho, g)
        Function that calculates the aeroplane mass ratio mug. See page 41 
        of EASA CS - VLA Easy Access. In particular, 
        
        CS - VLA 341 Gust load factors
        In the abscence of a more rational analysis, the gust load factors
        may be computed as follows: 
            
                      0.5 * rho0 * V * a * Kg * Ude
            n = 1.0 + -----------------------------
                               ( W / S )
            where 
                  0.88 * mu_g
            Kg = ------------- = gust alleviation factor 
                  5.3 + mu_g
                  
                       2 * (M / S)
            mu_g = ------------------- = aeroplane mass ratio 
                   rho * MGC * CL_alfa
                   
            U_de    = derived gust velocities referred to in CS - VLA 333(c) [m/s]
            rho0    = density of air at sea level [kg/m^3]
            rho     = density of air at operative altitude [kg/m^3]
            M/S     = wing loading [kg/m^2]
            MGC     = mean geometric chord [m]
            g       = acceleration due to gravity [m/s^2]
            V       = aeroplane equivalent speed [m/s]
            CL_alfa = slope of the aeroplane normal force coefficient curve CNA
                      in [1/rad] if the gust loads are applied to the wings and
                      horizontal tail surfaces simultaneously by a rational me-
                      thod. The wing lift curve slope CL_alfa in [1/rad] may be
                      used when the gust load is applied to the wings only and 
                      the horizontal tail gust loads are treated as a separate 
                      condition.

        Returns
        -------
        mu_g : float 
            Aeroplane mass ratio.

        """
        Kg = ( 0.88 * mu_g ) / ( 5.3 + mu_g )
        
        return Kg
        
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # GUST LOAD FACTORS 
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def calc_n_gust(rho0, VC, VD, a, kg, Ude, WS, case): 
        """
        Gust load factors = n_gust = calc_n_gust(rho0, V, a, kg, Ude, WS)
        Function that calculates the relevant gust load factors. See page 41 
        of EASA CS - VLA Easy Access. In particular, 
        
        CS - VLA 341 Gust load factors
        In the abscence of a more rational analysis, the gust load factors
        may be computed as follows: 
            
                      0.5 * rho0 * V * a * Kg * Ude
            n = 1.0 + -----------------------------
                               ( W / S )
            where 
                  0.88 * mu_g
            Kg = ------------- = gust alleviation factor 
                  5.3 + mu_g
                  
                       2 * (M / S)
            mu_g = ------------------- = aeroplane mass ratio 
                   rho * MGC * CL_alfa
                   
            U_de    = derived gust velocities referred to in CS - VLA 333(c) [m/s]
            rho0    = density of air at sea level [kg/m^3]
            rho     = density of air at operative altitude [kg/m^3]
            M/S     = wing loading [kg/m^2]
            MGC     = mean geometric chord [m]
            g       = acceleration due to gravity [m/s^2]
            V       = aeroplane equivalent speed [m/s]
            CL_alfa = slope of the aeroplane normal force coefficient curve CNA
                      in [1/rad] if the gust loads are applied to the wings and
                      horizontal tail surfaces simultaneously by a rational me-
                      thod. The wing lift curve slope CL_alfa in [1/rad] may be
                      used when the gust load is applied to the wings only and 
                      the horizontal tail gust loads are treated as a separate 
                      condition.

        Parameters
        ----------
        rho0 : float
            Air density at sea level in [kg/m^3].
        VC : float
            Cruise airspeed vector in [m/s]. 
        VD : float
            Dive airspeed vector in [m/s]. 
        a : float
            Lift curve slope in [1/rad].
        kg : float
            Gust alleviation factor.
        Ude : float
            Corresponding gust reference speed.
        WS : float
            Wing loading in [Pa].
        case : flag
            Gust load calculation flag: 
                case1 --> 'positive'
                case2 --> 'negative'                    

        Returns
        -------
        n_gust : TYPE
            Gust load factors.

        """
        if (Ude == 15.24):
            V = VC
        elif (Ude == 7.62):
            V = VD
        
        if (case == "positive"):
            numerator   = 0.5 * rho0 * V * a * kg * Ude
            denominator = WS
            n_gust      = 1.0 + (numerator / denominator)
        elif (case == "negative"):
            numerator   = 0.5 * rho0 * V * a * kg * Ude
            denominator = WS
            n_gust      = 1.0 - (numerator / denominator)
        
        return n_gust
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # GUST ENVELOPE DIAGRAM 
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def gust_envelope_diagram(n_from0toS, n_fromStoA, n_fromAtoC, n_fromCtoD, n_fromDto0,\
                    V_from0toS, V_fromStoA, V_fromAtoC, V_fromCtoD, V_fromDto0,\
                    n_from0toSinv, n_fromSinvtoG, n_fromGtoF, n_fromFto0,\
                    V_from0toSinv, V_fromSinvtoG, V_fromGtoF, V_fromFto0,\
                    Reg, Aircraft_name, n,\
                    nS, VS, nA, VA, nC, VC, nD, VD,\
                    nSinv, VSinv, nG, VG, nF, VF,\
                    ng_pos_cruise, ng_neg_cruise, ng_pos_dive, ng_neg_dive,\
                    ng_fromCtoD, ng_fromFtoE, vg_cruise, vg_dive):
        """
        

        Parameters
        ----------
        n_from0toS : TYPE
            DESCRIPTION.
        n_fromStoA : TYPE
            DESCRIPTION.
        n_fromAtoC : TYPE
            DESCRIPTION.
        n_fromCtoD : TYPE
            DESCRIPTION.
        n_fromDto0 : TYPE
            DESCRIPTION.
        V_from0toS : TYPE
            DESCRIPTION.
        V_fromStoA : TYPE
            DESCRIPTION.
        V_fromAtoC : TYPE
            DESCRIPTION.
        V_fromCtoD : TYPE
            DESCRIPTION.
        V_fromDto0 : TYPE
            DESCRIPTION.
        n_from0toSinv : TYPE
            DESCRIPTION.
        n_fromSinvtoG : TYPE
            DESCRIPTION.
        n_fromGtoF : TYPE
            DESCRIPTION.
        n_fromFto0 : TYPE
            DESCRIPTION.
        V_from0toSinv : TYPE
            DESCRIPTION.
        V_fromSinvtoG : TYPE
            DESCRIPTION.
        V_fromGtoF : TYPE
            DESCRIPTION.
        V_fromFto0 : TYPE
            DESCRIPTION.
        Reg : TYPE
            DESCRIPTION.
        Aircraft_name : TYPE
            DESCRIPTION.
        n : TYPE
            DESCRIPTION.
        nS : TYPE
            DESCRIPTION.
        VS : TYPE
            DESCRIPTION.
        nA : TYPE
            DESCRIPTION.
        VA : TYPE
            DESCRIPTION.
        nC : TYPE
            DESCRIPTION.
        VC : TYPE
            DESCRIPTION.
        nD : TYPE
            DESCRIPTION.
        VD : TYPE
            DESCRIPTION.
        nSinv : TYPE
            DESCRIPTION.
        VSinv : TYPE
            DESCRIPTION.
        nG : TYPE
            DESCRIPTION.
        VG : TYPE
            DESCRIPTION.
        nF : TYPE
            DESCRIPTION.
        VF : TYPE
            DESCRIPTION.
        ng_pos_cruise : TYPE
            DESCRIPTION.
        ng_neg_cruise : TYPE
            DESCRIPTION.
        ng_pos_dive : TYPE
            DESCRIPTION.
        ng_neg_dive : TYPE
            DESCRIPTION.
        ng_fromCtoD : TYPE
            DESCRIPTION.
        ng_fromFtoE : TYPE
            DESCRIPTION.
        vg_cruise : TYPE
            DESCRIPTION.
        vg_dive : TYPE
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
        fig  = plt.figure()
        plt.plot(V_from0toS,    n_from0toS,    color="red", linewidth=0.8, linestyle = "solid", zorder=1)
        plt.plot(V_fromStoA,    n_fromStoA,    color="red", linewidth=0.8, linestyle = "solid", zorder=1)
        plt.plot(V_fromAtoC,    n_fromAtoC,    color="red", linewidth=0.8, linestyle = "solid", zorder=1)
        plt.plot(V_fromCtoD,    n_fromCtoD,    color="red", linewidth=0.8, linestyle = "solid", zorder=1)
        plt.plot(V_fromDto0,    n_fromDto0,    color="red", linewidth=0.8, linestyle = "solid", zorder=1)
        plt.plot(V_from0toSinv, n_from0toSinv, color="red", linewidth=0.8, linestyle = "solid", zorder=1)
        plt.plot(V_fromSinvtoG, n_fromSinvtoG, color="red", linewidth=0.8, linestyle = "solid", zorder=1)
        plt.plot(V_fromGtoF,    n_fromGtoF,    color="red", linewidth=0.8, linestyle = "solid", zorder=1)
        plt.plot(V_fromFto0,    n_fromFto0,    color="red", linewidth=0.8, linestyle = "solid", zorder=1)

        plt.plot(vg_cruise,    ng_pos_cruise,    color="black", linewidth=1.0, linestyle = "dashed", zorder=1)
        plt.plot(vg_cruise,    ng_neg_cruise,    color="black", linewidth=1.0, linestyle = "dashed", zorder=1)
        plt.plot(vg_dive,      ng_pos_dive,      color="black", linewidth=1.0, linestyle = "dashed", zorder=1)
        plt.plot(vg_dive,      ng_neg_dive,      color="black", linewidth=1.0, linestyle = "dashed", zorder=1)

        plt.plot(V_fromCtoD,      ng_fromCtoD,      color="black", linewidth=1.0, linestyle = "dashed", zorder=1)
        plt.plot(V_fromFto0,      ng_fromFtoE,      color="black", linewidth=1.0, linestyle = "dashed", zorder=1)
          
        plt.scatter(VS, nS, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
        # fontsize or size
        # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
        # 'x-large', 'xx-large'}
        # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
        # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
        # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
        label=r"Point S"
        plt.text(VS, nS, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VA, nA, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
        label=r"Point A"
        plt.text(VA, nA, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VC, nC, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
        label=r"Point C"
        plt.text(VC, nC, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VD, nD, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
        label=r"Point D"
        plt.text(VD-2.0, nD, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VSinv, nSinv, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point S inverted", zorder=2)
        label=r"Point S inverted"
        plt.text(VSinv, nSinv, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VG, nG, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
        label=r"Point G"
        plt.text(VG, nG, label, fontdict=None, fontsize='x-small', fontweight='bold')
        plt.scatter(VF, nF, s=None, c='black', marker='.', cmap=None, norm=None,\
                    vmin=None, vmax=None, alpha=1, linewidths=None,\
                    edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
        label=r"Point F"
        plt.text(VF, nF, label, fontdict=None, fontsize='x-small', fontweight='bold')
        
        plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
        plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
        # plt.ylim(nmin - 0.5, nmax + 0.5)
        # plt.xlim(0.0, VD + 10.0)
        # plt.legend(loc="best")
        plt.title(r'Gust envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
        plt.grid(True, linestyle='-.', which="both")
        plt.minorticks_on()
        name_figure = 'GustEnvelope' + str(n) + '.pdf'
        plt.savefig(name_figure, bbox_inches='tight')
        plt.show()
        
        return fig    
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # V - n DIAGRAM 
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def final_envelope(V_fe_from0toS, n_fe_from0toS, V_fe_fromStoA,\
                       n_fe_fromStoA, n_fe_fromAtoGust1, V_fe_fromAtoGust1,\
                       n_fe_fromGust1toC, V_fe_fromGust1toC, n_fe_fromCtoGust2,\
                       V_fe_fromCtoGust2, n_fe_fromGust2toD, V_fe_fromGust2toD,\
                       n_fe_fromDto0, V_fe_fromDto0, n_fe_from0toSinv,\
                       V_fe_from0toSinv, n_fe_fromSinvtoG, V_fe_fromSinvtoG,\
                       n_fe_fromGtoGust1, V_fe_fromGtoGust1, n_fe_fromGust1toF,\
                       V_fe_fromGust1toF, n_fe_fromFtoE, V_fe_fromFtoE,\
                       n_fe_fromEto0, V_fe_fromEto0, n_fe_fromGtoF,\
                       V_fe_fromGtoF, n_fe_fromAtoC, V_fe_fromAtoC,\
                       n_fe_fromCtoD, V_fe_fromCtoD, Reg, Aircraft_name, n,\
                       pos_case_flag, neg_case_flag):
        """
        

        Parameters
        ----------
        V_fe_from0toS : TYPE
            DESCRIPTION.
        n_fe_from0toS : TYPE
            DESCRIPTION.
        V_fe_fromStoA : TYPE
            DESCRIPTION.
        n_fe_fromStoA : TYPE
            DESCRIPTION.
        n_fe_fromAtoGust1 : TYPE
            DESCRIPTION.
        V_fe_fromAtoGust1 : TYPE
            DESCRIPTION.
        n_fe_fromGust1toC : TYPE
            DESCRIPTION.
        V_fe_fromGust1toC : TYPE
            DESCRIPTION.
        n_fe_fromCtoGust2 : TYPE
            DESCRIPTION.
        V_fe_fromCtoGust2 : TYPE
            DESCRIPTION.
        n_fe_fromGust2toD : TYPE
            DESCRIPTION.
        V_fe_fromGust2toD : TYPE
            DESCRIPTION.
        n_fe_fromDto0 : TYPE
            DESCRIPTION.
        V_fe_fromDto0 : TYPE
            DESCRIPTION.
        n_fe_from0toSinv : TYPE
            DESCRIPTION.
        V_fe_from0toSinv : TYPE
            DESCRIPTION.
        n_fe_fromSinvtoG : TYPE
            DESCRIPTION.
        V_fe_fromSinvtoG : TYPE
            DESCRIPTION.
        n_fe_fromGtoGust1 : TYPE
            DESCRIPTION.
        V_fe_fromGtoGust1 : TYPE
            DESCRIPTION.
        n_fe_fromGust1toF : TYPE
            DESCRIPTION.
        V_fe_fromGust1toF : TYPE
            DESCRIPTION.
        n_fe_fromFtoE : TYPE
            DESCRIPTION.
        V_fe_fromFtoE : TYPE
            DESCRIPTION.
        n_fe_fromEto0 : TYPE
            DESCRIPTION.
        V_fe_fromEto0 : TYPE
            DESCRIPTION.
        n_fe_fromGtoF : TYPE
            DESCRIPTION.
        V_fe_fromGtoF : TYPE
            DESCRIPTION.
        n_fe_fromAtoC : TYPE
            DESCRIPTION.
        V_fe_fromAtoC : TYPE
            DESCRIPTION.
        n_fe_fromCtoD : TYPE
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
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, n_fe_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, n_fe_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, n_fe_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, n_fe_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, n_fe_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, n_fe_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, n_fe_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], n_fe_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], n_fe_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], n_fe_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()    
            
        elif (pos_case_flag == "Case1") and (neg_case_flag == "Case2_inverted"):
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, n_fe_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, n_fe_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, n_fe_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, n_fe_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     n_fe_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, n_fe_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], n_fe_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()                
            
        elif (pos_case_flag == "Case1") and (neg_case_flag == "Case3_inverted"):
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, n_fe_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, n_fe_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, n_fe_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, n_fe_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     n_fe_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, n_fe_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], n_fe_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()                
            
        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case1_inverted"):
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     n_fe_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     n_fe_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, n_fe_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, n_fe_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], n_fe_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], n_fe_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, n_fe_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], n_fe_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], n_fe_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], n_fe_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()    
            
        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case2_inverted"):
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     n_fe_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     n_fe_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     n_fe_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], n_fe_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], n_fe_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, n_fe_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], n_fe_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()                
            
        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case3_inverted"):
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoC,     n_fe_fromAtoC,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoD,     n_fe_fromCtoD,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     n_fe_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromAtoC[-1], n_fe_fromAtoC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromAtoC[-1], n_fe_fromAtoC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromCtoD[-1]-2.0, n_fe_fromCtoD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromCtoD[-1], n_fe_fromCtoD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()                
            
        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case1_inverted"):
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, n_fe_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, n_fe_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, n_fe_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, n_fe_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoGust1, n_fe_fromGtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toF, n_fe_fromGust1toF, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, n_fe_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], n_fe_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGust1toF[-1], n_fe_fromGust1toF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toF[-1], n_fe_fromGust1toF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()    
            
        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case2_inverted"):
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, n_fe_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, n_fe_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, n_fe_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, n_fe_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     n_fe_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, n_fe_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], n_fe_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()                
            
        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case3_inverted"):
            
            fig  = plt.figure()
            plt.plot(V_fe_from0toS,     n_fe_from0toS,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromStoA,     n_fe_fromStoA,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromAtoGust1, n_fe_fromAtoGust1, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust1toC, n_fe_fromGust1toC, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromCtoGust2, n_fe_fromCtoGust2, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGust2toD, n_fe_fromGust2toD, color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromDto0,     n_fe_fromDto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            
            plt.plot(V_fe_from0toSinv,  n_fe_from0toSinv,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromSinvtoG,  n_fe_fromSinvtoG,  color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromGtoF,     n_fe_fromGtoF,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromFtoE,     n_fe_fromFtoE,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)
            plt.plot(V_fe_fromEto0,     n_fe_fromEto0,     color="red", linewidth=1.0, linestyle = "solid", zorder=1)

            label=r"Point S"    
            plt.text(V_fe_from0toS[-1], n_fe_from0toS[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')          
            plt.scatter(V_fe_from0toS[-1], n_fe_from0toS[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point S", zorder=2)
            # fontsize or size
            # float or {'xx-small', 'x-small', 'small', 'medium', 'large',
            # 'x-large', 'xx-large'}
            # fontweight or weight {a numeric value in range 0-1000, 'ultralight', 
            # 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold',
            # 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
            label=r"Point A"
            plt.text(V_fe_fromStoA[-1], n_fe_fromStoA[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromStoA[-1], n_fe_fromStoA[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point A", zorder=2)
            label=r"Point C"
            plt.text(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust1toC[-1], n_fe_fromGust1toC[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point C", zorder=2)
            label=r"Point D"
            plt.text(V_fe_fromGust2toD[-1]-2.0, n_fe_fromGust2toD[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGust2toD[-1], n_fe_fromGust2toD[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point D", zorder=2)
            label=r"Point S inverted"
            plt.text(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_from0toSinv[-1], n_fe_from0toSinv[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point G", zorder=2)
            label=r"Point G"
            plt.text(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromSinvtoG[-1], n_fe_fromSinvtoG[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point F"
            plt.text(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromGtoF[-1], n_fe_fromGtoF[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            label=r"Point E"
            plt.text(V_fe_fromFtoE[-1]-2.0, n_fe_fromFtoE[-1], label, fontdict=None, fontsize='x-small', fontweight='bold')
            plt.scatter(V_fe_fromFtoE[-1], n_fe_fromFtoE[-1], s=None, c='black', marker='.', cmap=None, norm=None,\
                        vmin=None, vmax=None, alpha=1, linewidths=None,\
                        edgecolors=None, plotnonfinite=False, data=None, label="Point F", zorder=2)
            
            plt.xlabel(r'\textsc{Airspeed} ~ $V$ ~ $[m/s]$')  # x-label to the axes.
            plt.ylabel(r'\textsc{Load Factor} ~ $n$ ~ $[g]$') # y-label to the axes.
            # plt.ylim(nmin - 0.5, nmax + 0.5)
            # plt.xlim(0.0, VD + 10.0)
            # plt.legend(loc="best")
            plt.title(r'Final envelope diagram per ' + Reg + ' - ' + Aircraft_name) 
            plt.grid(True, linestyle='-.', which="both")
            plt.minorticks_on()
            name_figure = 'FinalEnvelope' + str(n) + '.pdf'
            plt.savefig(name_figure, bbox_inches='tight')
            plt.show()                
            
        
        return fig
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # FINAL ENVELOPE DIAGRAM - STORING POINTS
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def final_envelope_store_points(V_fe_from0toS, n_fe_from0toS, V_fe_fromStoA,\
                       n_fe_fromStoA, n_fe_fromAtoGust1, V_fe_fromAtoGust1,\
                       n_fe_fromGust1toC, V_fe_fromGust1toC, n_fe_fromCtoGust2,\
                       V_fe_fromCtoGust2, n_fe_fromGust2toD, V_fe_fromGust2toD,\
                       n_fe_fromDto0, V_fe_fromDto0, n_fe_from0toSinv,\
                       V_fe_from0toSinv, n_fe_fromSinvtoG, V_fe_fromSinvtoG,\
                       n_fe_fromGtoGust1, V_fe_fromGtoGust1, n_fe_fromGust1toF,\
                       V_fe_fromGust1toF, n_fe_fromFtoE, V_fe_fromFtoE,\
                       n_fe_fromEto0, V_fe_fromEto0, n_fe_fromGtoF,\
                       V_fe_fromGtoF, n_fe_fromAtoC, V_fe_fromAtoC,\
                       n_fe_fromCtoD, V_fe_fromCtoD,\
                       pos_case_flag, neg_case_flag):
        
        if   (pos_case_flag == "Case1") and (neg_case_flag == "Case1_inverted"):

            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }

        elif (pos_case_flag == "Case1") and (neg_case_flag == "Case2_inverted"):
            
            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0":     {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }

        elif (pos_case_flag == "Case1") and (neg_case_flag == "Case3_inverted"):
            
            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }

        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case1_inverted"):
            
            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }

        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case2_inverted"):
            
            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0":     {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }

        elif (pos_case_flag == "Case2") and (neg_case_flag == "Case3_inverted"):
            
            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }

        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case1_inverted"):
            
            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }

        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case2_inverted"):
            
            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0":     {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }

        elif (pos_case_flag == "Case3") and (neg_case_flag == "Case3_inverted"):
            
            Final_envelope = {
                "Final_envelope": {
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
                    "n_fe_fromEto0"    : {"Value": n_fe_fromEto0, "Unit": "g"}
                    }
                }
       
        return Final_envelope