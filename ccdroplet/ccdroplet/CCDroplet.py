## CCDroplet.py
## Burning Droplet in a Combustion Chamber model
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160303
##                                                                               
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, or   
## (at your option) any later version.                                        
##     
## This program is distributed in the hope that it will be useful,               
## but WITHOUT ANY WARRANTY; without even the implied warranty of              
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 
## GNU General Public License for more details.                             
##                                                                          
## Governing the ethics of using this program, I default to the Caltech Honor Code: 
## ``No member of the Caltech community shall take unfair advantage of        
## any other member of the Caltech community.''                               
##                                                         
## linkedin     : ernestyalumni                                                    
## wordpress    : ernestyalumni                                                    
############################################################################ 
import time

import cantera as ct
import sympy

from collections import namedtuple

import numpy as np
from scipy.optimize import newton
from scipy import integrate

try:
    import LiquidVaporEq
except:
    print "LiquidVaporEq not here!"

# Wrap up given physical parameters and inlet conditions into namedtuple classes

Chamber_Params = namedtuple('Chamber_parameters',['A_totinj','A_cc','l_cc'])
# 'A_totinj' - total cross sectional area for LIQUID FUEL only (through injection plate) 
# 'A_cc'     - combustion chamber cross-sectional area
# 'l_cc'     - length of combustion chamber

Inlet_Conds    = namedtuple('Inlet_Conditions',['TP','phi_g','phi_overall','v_d','D'])

Reaction_Tuple = namedtuple('Reaction_Tuple',['Fuelstr','Oxstr','n_F','n_Ox','Prod_str'])

# Be prepared to wrap up mass flow quantities, derived from inlet conditions and 
# Be prepared to wrap up results, namely 
Flow_Tuple     = namedtuple('Flow_Tuple',['dotm_l','dotm_Ox','dotm_g','v_g'])
T_gDv_d_Tuple  = namedtuple('T_gDv_d',['T_g','D','v_d'])


# ODE related only
NUM_STEPS = 1000
EPSILON = 0.001

###########################################################################
##### LOX-methane rocket combustion chamber engine
###########################################################################
# Then recall the equivalence ratio relation (cf. https://en.wikipedia.org/wiki/Air-fuel_ratio)
# phi = ( n_Fu/n_Ox )/( n_Fu/n_Ox)_st, so (n_Fu/n_Ox)=phi*(n_Fu/n_Ox)_st = 0.45*0.5

# cf. Turns Combustion pp. 400 Ch. 10
# cf. Priem and Heidmann [39] "Propellant Vaporization as a Design Criterion for Rocket-Engine Combustion Chambers" NASA Technical Report R-67, 1960
# Dipprey [7] "Liquid Rocket Engines" Chemistry in Space Research (R.F. Landel and A. Rembaum, eds.), Elsevier, NY pp. 464-597, 1972


def compute_flows(gas,reaction_tuple,chamber_params,inlet_conds,flow_inlet,F_over_O,rho_l,D,T_g):
    """
    compute_flows = compute_flows(gas,reaction_tuple,chamber_params,inlet_conds,flow_inlet,F_over_O,rho_l,D)
    """
    P        = inlet_conds.TP[1]
    D_0      = inlet_conds.D
    A_cc     = chamber_params.A_cc

    DOTm_l0  = flow_inlet.dotm_l
    DOTm_Ox0 = flow_inlet.dotm_Ox
    DOTm_g0  = flow_inlet.dotm_g
    V_g0     = flow_inlet.v_g

    dotm_lx  = DOTm_l0*(D/D_0)**3
    dotm_gx  = DOTm_g0 + DOTm_l0 - dotm_lx
    phi_gx   = (dotm_gx/DOTm_Ox0-1)*1/F_over_O

    gas.X    = phi_to_X(phi_gx,reaction_tuple)
    MW_gx    = gas.mean_molecular_weight # MW_gx = MW_g(phi_g(x)) amu
 
    v_gx     = dotm_gx*ct.gas_constant*T_g/(MW_gx*P*A_cc)

    return dotm_lx, dotm_gx, phi_gx, v_gx


def ode_creator(gas, reaction_tuple, inlet_conds, chamber_params, flow_inlet, h_fg,F_over_O_Phi1_mass,T_boil,rho_l,h_l):
    def ode_rhs(x, y):
        # unpack arguments
        P = inlet_conds.TP[1]

        dotm_l0  = flow_inlet.dotm_l
        dotm_Ox0 = flow_inlet.dotm_Ox
        dotm_g0  = flow_inlet.dotm_g 

        T_g, Dsq, v_d = y
        if Dsq >0.:
            Ddrop = np.sqrt(Dsq)
            if Ddrop >0.:
                # compute flow rates
                dotm_lx, dotm_gx, phi_gx, v_g = compute_flows(gas,reaction_tuple,chamber_params,inlet_conds,flow_inlet,F_over_O_Phi1_mass,rho_l,Ddrop,T_g)


                # Get gas properties from Cantera
                gas.X = phi_to_X(phi_gx,reaction_tuple)
                gas.TP = T_g,P                
                h_g = gas.enthalpy_mass # J/kg
                rho_g = gas.density

                
                # Calculate Re
                mu_g = gas.viscosity # [Pa-s] = [kg/m-s]
                Re_Dg = Ddrop*np.abs(v_d-v_g)*rho_g/mu_g
                
                # Calculate C_D
                C_D = 24./Re_Dg + 6./(1+np.sqrt(Re_Dg)) + 0.4
                
                # evaporation constant K
                EVAPK = evap_K(gas,reaction_tuple,T_g,P,phi_gx,h_fg,T_boil,rho_l,T_infty=T_g, dotm_Ox=dotm_Ox0, dotm_g=dotm_gx,dotm_l=dotm_lx)
                ddotm_gx_over_dx = 3/2.*dotm_l0/(inlet_conds.D)**3*Ddrop/v_d*EVAPK

                dphi_gx_over_dxRHS = 1/F_over_O_Phi1_mass*1/dotm_Ox0*ddotm_gx_over_dx


                dDsq_over_dxRHS = -EVAPK/v_d
                dv_d_over_dxRHS = (3.*C_D*rho_g*(v_g-v_d)*np.abs(v_g-v_d))/(4.*rho_l*Ddrop*v_d)

                
                dh_g_over_dphiforRHS = dh_g_over_dphi(gas, reaction_tuple, T_g,P,phi_gx )    
                dh_g_over_dT_gforRHS = dh_g_over_dT_g(gas, reaction_tuple, T_g,P,phi_gx )    

                dT_g_over_dxRHS = ((-1./dotm_gx)*ddotm_gx_over_dx*(h_g-h_l)-dh_g_over_dphiforRHS*dphi_gx_over_dxRHS)/dh_g_over_dT_gforRHS

                n = len(y) # should be 3
                dydt = np.zeros((n,1))

                dydt[0] = dT_g_over_dxRHS                
                dydt[1] = dDsq_over_dxRHS
                dydt[2] = dv_d_over_dxRHS

                return dydt
            else:
                n = len(y) # should be 3
                dydt = np.zeros((n,1))
                return dydt
        else:
            print "Check Dsq, i.e. D^2 (droplet diameter squared)"
            n = len(y) # should be 3
            dydt = np.zeros((n,1))
            return dydt
    return ode_rhs


def evap_K(gas,reaction_tuple,T_g,P,phi_g,h_fg,T_boil, rho_l ,T_infty, dotm_Ox, dotm_g,dotm_l):
    """
    evap_K = evap_K(gas, reaction_tuple, T_g,P,phi_g, h_fg, T_boil,rho_l)
    
    Calculates evaporation constant K = K(T_g,P,phi_g)

    INPUTS/Parameters
    -----------------
    T_g          - <float>

    P            - <float>
    pressure P in pascal

    phi_g        - <float>

    h_fg         - <float>
    [h_fg] = J/kg
    
    """

    X_gx = phi_to_X(phi_g,reaction_tuple) # "X_g(x)" i.e. phi_g(x) but converted to Cantera's .X form for mole fractions
    gas.X = X_gx
    gas.TP = T_g,P

    # use as a guess for T_f to be the adibatic flame temperature (cf. Polk, Turns)
    gas.equilibrate('HP')
    T_fguess = gas.T


    # Calculate average temperature in the 
    # "inner zone" (gaseous phase of fuel and oxidizer)
    avgT = (T_boil + T_fguess)*0.5

    # Set the state of the fuel and oxidizer at avgT and P
    gas.X      = reaction_tuple.Fuelstr+':1'
    gas.TP     = avgT, P
    k_gFuel    = gas.thermal_conductivity
    C_PFuel    = gas.cp
    Deltah_Fc  = gas.enthalpy_mass

    gas.X      = reaction_tuple.Oxstr+':1'
    gas.TP     = avgT, P
    k_gOx      = gas.thermal_conductivity
    Deltah_Oxc = gas.enthalpy_mass

    # Calculate k_g and C_P using Law and Williams' (1972) empirical laws for alkane droplets
    C_P = C_PFuel
    k_g = k_gFuel*0.4 + k_gOx*0.6

    # rho_l = rho_l(T_s) = rho_l(T_b) 
    # i.e. the droplet is assumed to have uniform temperature; that uniform temperature is the temperature at the surface, which is the boiling temperature at some particular pressure p!  Then you can look up the density for the fuel as a LIQUID

    gas.X = reaction_tuple.Prod_str
    gas.TP = avgT,P
    Deltah_Prc = gas.enthalpy_mass

    # NU i.e. \nu
    gas.X = X_gx

    Ox_index = gas.species_index(reaction_tuple.Oxstr)
    F_index  = gas.species_index(reaction_tuple.Fuelstr)

    NU = dotm_Ox/(dotm_l+(dotm_g-dotm_Ox))

    Deltah_c = Deltah_Fc + NU*Deltah_Oxc - (NU+1)*Deltah_Prc

    T_s     = T_boil


    B_OQ = ((Deltah_c/NU)+C_P*(T_infty-T_s))/(0+h_fg)

    return 8.* k_g/ (rho_l*C_P)*np.log( 1. + B_OQ)

def dh_g_over_dphi(gas,reaction_tuple,T_g,P,phi,epsilon = EPSILON ):
    """
    dh_g_over_dphi = dh_g_over_dphi(T_g,p,phi, Deltaphi)
    
    Calculates dh_g/dphi(T_g,P,\phi)

    INPUTS/Parameters
    -----------------
    gas      - <cantera._cantera.Solution or Cantera Phase Object>
    
    T_g      - <float>

    g_p      - <float>
     
    P        - <float>
    
    phi      - <float>

    reaction_tuple - <Reaction_Tuple, which is an instance of namedtuple>

    epsilon  - <float>
    \epsilon := \Delta \phi / \phi i.e. epsilon is the ratio between the 
    Delta phi change in phi, over phi 

    OUTPUT/Result
    -------------
    dh_g_over_dphi_result - <float>
    """
        
    phi_1 = phi * (1 + epsilon) # phi_1 slightly above phi
    phi_2 = phi * (1 - epsilon) # phi_2 slightly below phi
    gas.X = phi_to_X(phi_1, reaction_tuple)
    gas.TP = T_g,P
    gas.equilibrate('TP')
    h_1 = gas.enthalpy_mass # J/kg
    gas.X = phi_to_X(phi_2, reaction_tuple)
    gas.TP = T_g,P
    gas.equilibrate('TP')
    h_2 = gas.enthalpy_mass # J/kg
    return (h_1 - h_2)/(phi_1-phi_2)

def dh_g_over_dT_g(gas,reaction_tuple,T_g,P,phi,epsilon=EPSILON):
    """
    dh_g_over_dT_g_v2 
    This is hd_g_over_dT_g_v2, except equilibrate is used, though I 
    disagree that equilibrate is needed

    INPUTS/Parameters
    -----------------
    gas      - <cantera._cantera.Solution or Cantera Phase Object>

    T_g      - <float>
    
    P        - <float>

    phi      - <float>

    reaction_tuple - <Reaction_Tuple, which is an instance of namedtuple>

    epsilon  - <float>
    \epsilon := \Delta T_g / T_g i.e. epsilon is the ratio between the 
    Delta T_g change in T_g, over T_g 
    """
    T_g1   = T_g * (1 + epsilon )
    T_g2   = T_g * (1 - epsilon )
    gas.X  = phi_to_X(phi,reaction_tuple)
    gas.TP = T_g1, P
    h_1 = gas.enthalpy_mass # J/kg
    gas.X  = phi_to_X(phi,reaction_tuple)
    gas.TP = T_g2,P
    h_2 = gas.enthalpy_mass # J/kg
    return (h_1-h_2)/(T_g1-T_g2)


def dh_g_over_dT_g_v2(gas,reaction_tuple,T_g,P,phi, epsilon=0.00001):
    """
    dh_g_over_dT_g_v2 
    This is hd_g_over_dT_g_v2, except equilibrate is used, though I 
    disagree that equilibrate is needed

    INPUTS/Parameters
    -----------------
    gas      - <cantera._cantera.Solution or Cantera Phase Object>

    T_g      - <float>
    
    P        - <float>

    phi      - <float>

    reaction_tuple - <Reaction_Tuple, which is an instance of namedtuple>

    epsilon  - <float>
    \epsilon := \Delta T_g / T_g i.e. epsilon is the ratio between the 
    Delta T_g change in T_g, over T_g 
    """
    T_g1   = T_g * (1 + epsilon )
    T_g2   = T_g * (1 - epsilon )
    gas.X  = phi_to_X(phi,reaction_tuple)
    gas.TP = T_g1, P
    gas.equilibrate('TP')
    h_1 = gas.enthalpy_mass # J/kg
    gas.X  = phi_to_X(phi,reaction_tuple)
    gas.TP = T_g2,P
    gas.equilibrate('TP')
    h_2 = gas.enthalpy_mass # J/kg
    return (h_1-h_2)/(T_g1-T_g2)


def get_initial_flow_vals(gas,reaction_tuple,chamber_params,inlet_conds,rho_l,F_over_O):
    """
    get_initial_flow_vals = get_initial_flow_vals(gas,reaction_tuple,chamber_params,inlet_conds,rho_l,F_over_O)
    """
    P = inlet_conds.TP[1]

    DOTm_l0  = rho_l*inlet_conds.v_d*chamber_params.A_totinj # \dot{m}_l(0)
    
    DOTm_Ox0 = DOTm_l0/( (inlet_conds.phi_overall - inlet_conds.phi_g)*F_over_O)

    DOTm_g0  = (inlet_conds.phi_overall*F_over_O+1)*DOTm_Ox0 - DOTm_l0 
    

    # Equilibrate gas with inlet stoichiometry to get inlet T_g
    gas.X  = phi_to_X(inlet_conds.phi_g, reaction_tuple)
    gas.TP = inlet_conds.TP

    MW_g0  = gas.mean_molecular_weight # MW_g0 = MW_g(phi_g(0)) amu

    gas.equilibrate('HP') # adiabatic flame temperature
    T_g0 = gas.T # inlet T_g

    D_0    = inlet_conds.D
    v_d0   = inlet_conds.v_d

    v_g0   = DOTm_g0*ct.gas_constant*T_g0/(MW_g0*P*chamber_params.A_cc)
    

    FLOW_INLET = Flow_Tuple(dotm_l=DOTm_l0,dotm_Ox=DOTm_Ox0,dotm_g=DOTm_g0,v_g=v_g0)
    INITIAL_CONDS = T_gDv_d_Tuple(T_g=T_g0,D=D_0,v_d=v_d0)
    return FLOW_INLET, INITIAL_CONDS

def _process_liquid_Droplet(gas,liquid,reaction_tuple,T_boilo,P, Hrefchoice=None):
    """
    _process_liquid_Droplet

    _process_liquid_Droplet does 3 things:
    1. grab the Clausius-Clapeyron relation (CCEqn) 
    2. Calculate the new T_boil at the desired pressure P, and use that to
    calculate enthalpies
    3. Armed with (T_boil(P),P) \in \Sigma^liquid, grab liquid density from Cantera
    """
    # let's try to use Cantera as much as possible
    gas.X = reaction_tuple.Fuelstr+':1'
    gas.TP = T_boilo, ct.one_atm

    # Calculating L, latent heat of vaporization, heat required to put a 
    # liquid molecule into gas phase
    # complications arise from 1. difference choices of standard enthalpy that 
    # people choose (e.g. Cantera vs. NIST), so hence the following, accounting 
    # for these differences
    h_g      = gas.enthalpy_mass
    h_g_mole = gas.enthalpy_mole

    liquid.TP = T_boilo, ct.one_atm
    if liquid.X > 0.0:  # making sure you have purely liquid right below boiling
        h_g_vap = liquid.enthalpy_mass
        h_g_vap_mole = liquid.enthalpy_mole
        k = 1 # get the enthalpy of the liquid at temperature slightly below T_boil so that there's no more vapor, i.e. vapor fraction = 0
        deltaT = 0.01
        while liquid.X > 0.0:
            liquid.TP = T_boilo-k*deltaT, ct.one_atm
            k += 1
        h_l = liquid.enthalpy_mass
        h_l_mole = liquid.enthalpy_mole
    else:
        h_l = liquid.enthalpy_mass
        h_l_mole = liquid.enthalpy_mole
        k=1 # get the enthalpy of the liquid at temperature slightly above T_boil so that it's all vapor, i.e. vapor fraction = 1
        deltaT=0.01
        while liquid.X <= 0.0:
            liquid.TP = T_boilo+k*deltaT, ct.one_atm
            k += 1
        h_g_vap = liquid.enthalpy_mass
        h_g_vap_mole = liquid.enthalpy_mole
    if Hrefchoice == None:
        Delta_vapHchoice = raw_input("For ENTHALPY of VAPORIZATION, \Delta_{vap}H or Delta_vap H, type 'y' if you want to calculate Delta_vap H from the difference, at T_boil, from the file for gas with the built-in Phase Object for liquid ('Pure Fluid'), or from an arbitrary temperature above and below the T_boil off the liquid ('Pure Fluid') object: the danger is that the choice of reference for enthalpy of an EXTERNAL file (not built into Cantera) is DIFFERENT from reference enthalpy chosen for Cantera (which is more honest, as it's referenced off of the elements alone). (y/n)?")
    elif Hrefchoice != 'y' and Hrefchoice != 'n':
        Delta_vapHchoice = 'y'
    else:
        Delta_vapHchoice = Hrefchoice
    if Delta_vapHchoice =='y':
        h_fg = h_g - h_l
        h_fg_mole = h_g_mole-h_l_mole # J/kmol
    else:
        h_fg = h_g_vap - h_l
        h_fg_mole = h_g_vap_mole - h_l_mole
        # sanity check
    print "Enthalpy of vaporization per kg mass: %s \n " , h_fg
    print "Enthalpy of vaporization per kmol   : %s \n " , h_fg_mole

    # armed with h_fg_mole, latent heat of vaporization, we can get the 
    # Clausius-Clapeyron relation, CCeqn
    CCeqn =LiquidVaporEq.fit_ClausClapEqn(T_boilo,h_fg_mole)
    T_boil_lambd = sympy.lambdify(LiquidVaporEq.T,CCeqn.rhs-P)
    T_boil = newton( T_boil_lambd, T_boilo) # K

    ## Done with grabbing and using the Clausius-Clapeyron relation and 
    ## computing the new T_boil = T_boil(P) at the new desired pressure P.  
    gas.X = reaction_tuple.Fuelstr+':1'
    gas.TP = T_boil, P
    h_g      = gas.enthalpy_mass

    liquid.TP = T_boil, P
    if liquid.X > 0.0:  # making sure you have purely liquid right below boiling
        h_g_vap = liquid.enthalpy_mass
        k = 1 # get the enthalpy of the liquid at temperature slightly below T_boil so that there's no more vapor, i.e. vapor fraction = 0
        deltaT = 0.01
        while liquid.X > 0.0:
            liquid.TP = T_boil-k*deltaT, P
            k += 1
        h_l = liquid.enthalpy_mass # J/kg
        rho_l = liquid.density # kg/m^3
    else:
        rho_l = liquid.density # kg/m^3
        h_l = liquid.enthalpy_mass # J/kg
        k=1 # get the enthalpy of the liquid at temperature slightly above T_boil so that it's all vapor, i.e. vapor fraction = 1
        deltaT=0.01
        while liquid.X <= 0.0:
            liquid.TP = T_boil+k*deltaT, P
            k += 1
        h_g_vap = liquid.enthalpy_mass

    # Note: Use Cantera for fluids as much as possible, 
    #  if fluid/liquid phase is available.
    # Otherwise, external files for the gas phase, other than the 
    #  built-in Cantera data files such as "gri30.cti", 
    #  may have reference temperatures, reference enthalpies that are 
    #  NOT the same as Cantera.  

    # case of, we trust that gas and liquid both que off the same reference temperature, reference enthalpy
    if Delta_vapHchoice =='y': 
        h_fg = h_g-h_l

    # case of, we DON'T trust that gas and liquid both que off the same reference temperature, reference enthalpy, so we'll have to subtract away this H(To,p) difference
    else: 
        h_fg = h_g_vap - h_l # at new T_boil for P 
        
        # now we have to fix the reference enthalpy that the given gas cues off of 
        # to the reference enthalpy that the liquid cues off of
        # to obtain a h_l
        H_1 = h_g_vap
        H_2 = h_g
        Deltah_12 = h_g_vap - h_g
        h_l = h_l - Deltah_12
    return h_fg, T_boil, rho_l, h_l



def main(gas,liquid,species_name,reaction_tuple,chamber_params,inlet_conds,Hrefchoice=None,T_boilo=None):
    """
    main - main(gas,liquid,species_name,reaction_tuple,chamber_params,inlet_conds,Hrefchoice=None,T_boilo=None)

    INPUTS/Parameters
    -----------------
    gas                      - <cantera.Solution or Cantera Phase object>
    liquid                   - <cantera.Solution or Cantera Phase object>
    reaction_tuple           - <namedtuple defined as Reaction_tuple>
     Reaction_tuple(Fuelstr=<string>,
      Oxstr=<string>,
      n_F=<int>,  # number of moles of Fuel in stoichiometric reaction
      n_Ox=<int>, # number of moles of oxidizer in stoichiometric reaction 
      Prod_str=<string>) 
    chamber_params = CHAMBER_PARAMS - <namedtuple defined as Chamber_Params>
     Chamber_Params(A_totinj=<float>,
      A_cc=<float>
      l_cc=<float>
    inlet_conds    = INLET_CONDS - <namedtuple defined as Inlet_Conds>
     Inlet_Conds(TP=(<float>,<float>),
      phi_g=<float>,
      phi_overall=<float>,
      v_d=<float>,
      D=<float>)

    """
    P = inlet_conds.TP[1]

    Fu_ind    = gas.species_index(reaction_tuple.Fuelstr)
    Ox_ind    = gas.species_index(reaction_tuple.Oxstr) 

    # Grab Clausius-Clapeyron relation from NIST Chemistry Webbook
    if T_boilo == None:
        try:
            dat = LiquidVaporEq.cleaned_Phase_data(species_name)
            print "Done with accessing online the NIST Chemistry Website. \n"
            T_boil0 = dat.T_boil['Value'] # T_boil0, T_boil at "standard" pressure reference of ct.one_atm, 1 atm
        except:
            print "Check if species name (species_name) is in the NIST Chemistry Webbook \n or if LiquidVaporEq file is there\n"
            T_boil0 = raw_input("Meanwhile, enter a T_boil at 1 atm:")
            T_boil0 = float(T_boil0)
    else:
        T_boil0 = T_boilo

    h_fg, T_boil, rho_l, h_l = _process_liquid_Droplet(gas,liquid,reaction_tuple,T_boil0,P,Hrefchoice=Hrefchoice)


    # Compute (F/O)_{\Phi=1}
    gas.X = phi_to_X(1, reaction_tuple)
    F_over_O_Phi1_mass = gas.Y[Fu_ind]/gas.Y[Ox_ind]


    ## Compute flow rates at inlet
    FLOW_INLET,INITIAL_CONDS=get_initial_flow_vals(gas,reaction_tuple,chamber_params,inlet_conds,rho_l,F_over_O_Phi1_mass)


    #### Call ODE
    cc_ode = ode_creator(gas,reaction_tuple,inlet_conds,chamber_params,FLOW_INLET,h_fg,F_over_O_Phi1_mass, T_boil,rho_l,h_l) 
    
    cc_int_ode = integrate.ode(cc_ode)

    # at this point, there are many ODE integrators that can be used.
    # Matlab's ode45 is scipy's dopri5, a Runge-Kutta method
    # there's also vode, Real-valued Variable-coefficient ODE solver

    cc_int_ode.set_integrator("dopri5")

    # set initial conditions
    T_g0   = INITIAL_CONDS.T_g
    Dsq0   = INITIAL_CONDS.D**2
    v_d0   = INITIAL_CONDS.v_d
    
    cc_int_ode.set_initial_value([T_g0,Dsq0,v_d0],0) # x_0 = 0

    # Additional Python step: create vectors to store data                                   
    x = np.zeros(NUM_STEPS+1)
    Temp_g              = np.zeros(NUM_STEPS+1)
    Dsquared            = np.zeros(NUM_STEPS+1)
    v_droplet           = np.zeros(NUM_STEPS+1)
    
    x[0] = 0.
    Temp_g[0]              = T_g0
    Dsquared[0]            = Dsq0
    v_droplet[0]           = v_d0

    


    ## Integrate the ODE(s) across each delta_x timestep                                           
    delta_x = (chamber_params.l_cc -0.)/(float(NUM_STEPS))

#    return cc_int_ode #, x, , Temp_g, Dsquared, v_droplet
    
    print "Initiating quasi-static inlet flows into combustion chamber, and timing it \n"
    start_time=time.time()
    k = 1

    while cc_int_ode.successful() and k <= NUM_STEPS and cc_int_ode.y[1]>=0.: # cc_int_ode.y[1]>=0. is a check if the droplet diameter isn't negative i.e. droplet has evaporated away and died. # and cc_int_ode.y[0] > 0: # check T_g >0K (no negative temperatures)
        try:
            cc_int_ode.integrate(cc_int_ode.t + delta_x)

        # Store the results to plot later                                                            
            x[k]                   = cc_int_ode.t
            Temp_g[k]              = cc_int_ode.y[0]
            Dsquared[k]            = cc_int_ode.y[1]
            v_droplet[k]           = cc_int_ode.y[2]

            k += 1
        except RuntimeError:
            print "Runtime Error \n"

    # for the "cleaned" data points, we cut off if the droplet dies, i.e. when 
    # droplet diameter is less than 0
    x_clean = x[:k-1]
    Temp_g_clean = Temp_g[:k-1]
    Dsquared_clean = Dsquared[:k-1]
    v_droplet_clean = v_droplet[:k-1]

    print "Completed loop, integrating over entire chamber.\n Took --- %s seconds ---" % (time.time()-start_time)

    #    return cc_int_ode, x,Temp_g, Dsquared, v_droplet
    results = {'ODEobj'  :cc_int_ode,
               'x'       :x,
               'xclean'  :x_clean,
               'rawdat'  :T_gDv_d_Tuple(T_g=Temp_g,D=np.sqrt(Dsquared),v_d=v_droplet),
               'cleandat':T_gDv_d_Tuple(T_g=Temp_g_clean,D=np.sqrt(Dsquared_clean),v_d=v_droplet_clean) #, 
#               'flow_inlet':FLOW_INLET
               }
                   
    return results


###############################################################################################
##### Miscellaneous functions
###############################################################################################

def phi_to_X(phi, reaction_tuple):
    """
    phi_to_X = phi_to_X(phi,reaction_tuple)
    From an equivalence ratio, return X, species mole fraction, specified as a string, 
    that Cantera can use; 
    cf. http://www.cantera.org/docs/sphinx/html/cython/thermo.html#cantera.ThermoPhase.X
 
    INPUTS/Parameters
    -----------------
    phi            - <float>
    equivalence ratio

    reaction_tuple - <Reaction_tuple, a declared named_tuple>

    OUTPUT/Result
    -------------
    X_str         - <string>

    """
    n_F_over_n_Ox = phi*reaction_tuple.n_F/float(reaction_tuple.n_Ox) # \frac{n_F}{n_{Ox}} = \phi*\left( \frac{n_F}{n_{Ox}} \right)_{st}
    return (
        reaction_tuple.Fuelstr
        + r':'
        + str(n_F_over_n_Ox)
        + ' '
        + reaction_tuple.Oxstr
        + ':1'
    )
    

if __name__ == "__main__":
    gri30_gas = ct.Solution("gri30.cti","gri30_mix")
    ctHep = ct.Heptane()

