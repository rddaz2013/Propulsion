## CombustionChamber.py
## Combustion Chamber model
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160225
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
import cantera as ct
import sympy

from collections import namedtuple

import numpy as np
from scipy.optimize import newton

try:
    import LiquidVaporEq
except:
    print "LiquidVaporEq not here!"

try:
    import Liquiddensity
except:
    print "Liquiddensity not here!"

A_CS        = 0.157                    # m^2 combustion chamber cross-sectional area
L_CC        = 0.75                     # m length of combustion chamber

A_TOTINJCS  = 0.0157                   # m^2 total fuel injector cross-sectional area



D_IDROPLET  = (30,50,80,100,200)       # \mu m (microns) initial (injected) droplet diameter 
T_I         = 600                      # K inlet gas temperature 
P_CC        = 3.4474                   # MPa Combustion Chamber Pressure


PHI_OVERALL = 1.139                    # overall equivalence ratio 

# Assume oxidizer injected as gas and fuel partly as gas with initial equivalence ratio 0.45
PHI_I       = 0.45                     # initial equivalence ratio
V_0         = 10                       # m/s initial droplet velocity

# Wrap up these given physical parameters and inlet conditions into namedtuple classes

Chamber_Params = namedtuple('Chamber_parameters',['A_cc','l_cc','A_totinj'])
# 'A_cc' - combustion chamber cross-sectional area
# 'l_cc' - length of combustion chamber
# 
CHAMBER_PARAMS = Chamber_Params(A_cc=A_CS, l_cc=L_CC, A_totinj=A_TOTINJCS)

Inlet_Conds  = namedtuple('Inlet_Conditions',['TP','phi_g','phi_overall','D','v_d'])
INLET_CONDS0 = Inlet_Conds(TP=(T_I,P_CC*1000000),phi_g=PHI_I,phi_overall=PHI_OVERALL,D=D_IDROPLET[0], v_d=V_0) 

###########################################################################
##### LOX-methane rocket
###########################################################################
# Suppose stoichiometric reaction for LOX-methane is
# CH4 + 2O2 -> CO2+2H2O

# Then recall the equivalence ratio relation (cf. https://en.wikipedia.org/wiki/Air-fuel_ratio)
# phi = ( n_Fu/n_Ox )/( n_Fu/n_Ox)_st, so (n_Fu/n_Ox)=phi*(n_Fu/n_Ox)_st = 0.45*0.5

# reaction_dict = dict.fromkeys(['Fuelstr','Oxstr','n_F','n_Ox','Prod_str'])

Reaction_tuple = namedtuple('Reaction_tuple',['Fuelstr','Oxstr','n_F','n_Ox','Prod_str'])

CH4_TUP = Reaction_tuple(Fuelstr='CH4',Oxstr='O2',n_F=1,n_Ox=2,Prod_str='CO2:1 H2O:2')

def ode_rhs(x, DDropsq, rho,u_g,A, gas):
    """
    ode_rhs = ode_rhs(DDrop)
    
    INPUTS/Parameters
    DDrop - <float>
    Diameter of Droplet

    rho   - <float>
    mass density of gas
    
    u_g     - <float>
    velocity of gas
    
    A     - <float>
    cross-sectional area

    gas   - <cantera._cantera.Solution or cantera Phase object>
    """
    assert DDrop > 0.
    if DDrop > 0.:
#        ddot = rho*u*A
        rho_g = gas.density # kg/m^3
        T_g, p_g, X_g = gas.TPX
        MW_g = gas.mean_molecular_weight

        A = A_CS
        # v_g; p_g = P = 3.44474 MPa
        v_g = dotm_g0*ct.gas_constant*T_g/(MW_g*p_g*A)   # ct.gas_constant = 8314 J/kmol-K

        # Re
        mu_g = gas.viscosity # [Pa-s] = [kg/m-s]
        Re_Dg =  DDrop*np.abs(v_d-v_g)*rho_l/mu_g
        
        # C_D
        C_D = 24./Re_Dg + 6./(1 + np.sqrt(Re_Dg)) + 0.4

        # evaporation constant K
        EVAPK = evap_K(gas,T_g,p_g,phi,reaction_tuple,species_name)
        
        # dm_l/dx and dphi/dx
        dm_l_over_dx = -3/2.*dotm_l0/D_0**3*DDrop*EVAPK/v_d

        dm_g_over_dx = - dm_l_over_dx # by continuity

        ## Calculate (F/O)_{phi=1}
        Fuelind = gas.species_index(reaction_tuple.Fuelstr)
        Oxind   = gas.species_index(reaction_tuple.Oxstr)
        F_over_O_phi1_mass = gas.Y[Fuelind]/gas.Y[Oxind]

        dphi_g_over_dx = 1/F_over_O_phi1_mass*dm_g_over_dx
        
        dh_g_over_dphiforRHS = dh_g_over_dphi(gas, T_g,P,phi_g, reaction_tuple)    
        dh_g_over_dT_gforRHS = dh_g_over_dT_g(gas, T_g,P,phi_g, reaction_tuple)    

        
        # RHS 
        dDsq_over_dxRHS = -EVAPK/v_d
        dv_d_over_dxRHS = (3.*C_D*rho_g*(v_g-v_d)*np.abs( v_g-v_d))/(4.*rho_l*DDrop*v_d)
        dT_g_over_dxRHS = ((-1./dotm_g)*dm_g_over_dx*(h_g-h_l)-dh_g_over_dphiforRHS*dphi_g_over_dx)/dh_g_over_dT_gforRHS


        RHS = [dDsq_over_dxRHS, dv_d_over_dxRHS, dT_g_over_dxRHS]
        
    else:
        RHS = [0,0,0]
    return RHS

def evap_K(gas,T_g,P,phi, reaction_tuple, species_name="methane" ):
    """
    evap_K = evap_K(T_g,P,phi)
    
    Calculates evaporation constant K = K(T_g,P,phi)

    INPUTS/Parameters
    -----------------
    T_g          - <float>

    P            - <float>
    pressure P in pascal

    phi          - <float>
    
    species_name - <string>
    """
    try:
        dat = LiquidVaporEq.cleaned_Phase_data(species_name)
    except:
        print "Check if species name (species_name) is in the NIST Chemistry Webbook \n"    
    # Calculate T_b from Clausius-Clapeyron equation
    # Remember, that the Clausius-Clapeyron equation uses, to fix the integration
    # constant, the boiling temperature measured at 1 atm, as a standard.  
    # Clausius-Clapeyon equation tells us how the temperature at which the substances
    # boils or vaporizes changes with pressure
    DeltaH=float(dat.DeltaH_vap[0][2][0])*1000*1000 # kJ/mol -> J/mol -> J/kmol
    CCeqn =LiquidVaporEq.fit_ClausClapEqn(dat.T_boil['Value'],DeltaH)
    T_boil_lambd = sympy.lambdify(LiquidVaporEq.T,CCeqn.rhs-P)
    T_boil = newton( T_boil_lambd, dat.T_boil['Value']) # K

    X_gx = phi_to_X(phi,reaction_tuple) # "X_g(x)" i.e. phi_g(x) but converted to Cantera's .X form for mole fractions
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
    rho_l = Liquiddensity.fluid_density(T_boil,P,species_name)
    
    gas.X = reaction_tuple.Prod_str
    gas.TP = avgT,P
    Deltah_Prc = gas.enthalpy_mass

    # NU i.e. \nu
    gas.X = reaction_tuple.Fuelstr+r":"+str(reaction_tuple.n_F)+r" "+reaction_tuple.Oxstr+r":"+str(reaction_tuple.n_Ox)

    Ox_index = gas.species_index(reaction_tuple.Oxstr)
    F_index  = gas.species_index(reaction_tuple.Fuelstr)

    NU = gas.Y[Ox_index]/gas.Y[F_index]
    
    Deltah_c = Deltah_Fc + NU*Deltah_Oxc - (NU+1)*Deltah_Prc

    T_infty=T_I
    T_s    = T_boil


    # get heat of vaporization of fuel for the droplet
    h_Fg = float(dat.DeltaH_vap[0][2][0])/dat.MW*1000 # kJ/mol/(g/mol) -> kJ/kg

    B_OQ = ((Deltah_c/NU)+C_P*(T_infty-T_s))/(0+h_Fg*1000)

    K = -8.* k_g/ (rho_l*C_P)*np.log( 1. + B_OQ)
    
    #return B_OQ, K 
    return K

def dh_g_over_dphi(gas, T_g,P,phi, reaction_tuple, epsilon = 0.00001 ):
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
    dh_g_over_dphi_result = (h_1 - h_2)/(phi_1-phi_2)
    return dh_g_over_dphi_result

def dh_g_over_dT_g(gas,T_g,P,phi,reaction_tuple, epsilon=0.00001):
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
    dh_g_over_dT_g_result = (h_1-h_2)/(T_g1-T_g2)
    return dh_g_over_dT_g_result


def dh_g_over_dT_g_v2(gas,T_g,P,phi,reaction_tuple, epsilon=0.00001):
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
    dh_g_over_dT_g_result = (h_1-h_2)/(T_g1-T_g2)
    return dh_g_over_dT_g_result

def main(reaction_tuple=CH4_TUP,chamber_params=CHAMBER_PARAMS,inlet_conds=INLET_CONDS0):
    """
    main - main

    INPUTS/Parameters
    -----------------
    reaction_tuple = CH4_TUP - <namedtuple defined as Reaction_tuple>

    chamber_params = CHAMBER_PARAMS - <namedtuple defined as Chamber_Params>
    
    inlet_conds    = INLET_CONDS - <namedtuple defined as Inlet_Conds>

    """

    gas_gri30 = ct.Solution('gri30.cti','gri30_mix')
    Fu_ind    = gas_gri30.species_index(reaction_tuple.Fuelstr)
    Ox_ind    = gas_gri30.species_index(reaction_tuple.Oxstr) 



    # INLET CONDITIONS    
    ## Define constants and conditions at inlet
    DeltaH=float(dat.DeltaH_vap[0][2][0])*1000*1000 # kJ/mol -> J/mol -> J/kmol
    CCeqn =LiquidVaporEq.fit_ClausClapEqn(dat.T_boil['Value'],DeltaH)
    T_boil_lambd = sympy.lambdify(LiquidVaporEq.T,CCeqn.rhs-P)
    T_boil = newton( T_boil_lambd, dat.T_boil['Value']) # K
    rho_l = Liquiddensity.fluid_density(T_boil, inlet_conds.TP[1])

    ## Compute flow rates at inlet
    

    # density of pure methane (but remember methane boiling temperature at chamber pressure)
    gas_gri30.X = Fuelstr+r':1'
    gas_gri30.TP = T_I,P_CC*1000000
    rho_Fu_i    = gas_gri30.density # kg/m^3 mass density of fuel at inlet conditions

    # FLOW RATES at INLET
    # oxidizer injected as gas, fuel partly as gas
    gas_gri30.X = phi_to_X(PHI_I,Fuelstr=Fuelstr,Oxstr=Oxstr,stoich_Xratio=stoich_Xratio)
    gas_gri30.TP = T_I,P_CC*1000000 # (K, Pa)
    rho_vap_i   = gas_gri30.density # kg/m^3 total vapor (gas) mass density at inlet, initially
    
    dotm_VAP = rho_vap_i*U_I*A_TOTINJCS # \dot{m} = rho_vap * u * A
    gas_gri30.equilibrate('HP')
    INLET_T_g = gas_gri30.T
    
    # Equilibrate gas with inlet stoichiometry to get inlet T_g
    gas_gri30.X  = phi_to_X(inlet_conds.phi_g, reaction_tuple)
    gas_gri30.TP = inlet_conds.TP
    gas_gri30.equilibrate('HP') # adiabatic flame temperature
    T_g0 = gas_gri30.T # inlet T_g

    # Call ODE

    return gas_gri30

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
    X_str = reaction_tuple.Fuelstr+r':'+str(n_F_over_n_Ox)+' '+reaction_tuple.Oxstr+':1'
    return X_str




    

    
if __name__ == "__main__":
    gri30_gas = ct.Solution("gri30.cti","gri30_mix")
    

