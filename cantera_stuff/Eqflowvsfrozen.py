## Eqflowvsfrozen.py
## Equilibrium flow vs. Frozen flow
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160205
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
from sets import Set

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

try:
    from sound_speed import equilSoundSpeeds
except:
    print "sound_speed.py is not in the working file directory!  We need function equilSoundSpeeds"

# "load" or "grab" the species, and species only, with thermo and transport model, from h2o2.cti
H2O2spec = ct.Species.listFromFile('h2o2.cti')
# Indeed, the thermo and transport models are there, for each species, if one does
dir(H2O2spec[0])
H2O2spec[0].thermo
H2O2spec[0].transport

################################################################################
##### NASA 9 polynomials
################################################################################
# http://www.cantera.org/docs/sphinx/html/cti/species.html#mcbride2002 
# http://www.cantera.org/docs/sphinx/html/cti/species.html#the-nasa-9-coefficient-polynomial-parameterization
# NASA 9 polynomials are in nasa.cti
# https://github.com/Cantera/cantera/blob/master/data/inputs/nasa.cti
# on your local hard drive, nasa.cti's location (what file directory it's in, 
# depends on your installation method.  
# Using MacPorts on a Mac, I found the .cti files here:
# /opt/local/share/cantera/data
# Using Homebrew to install, I found it in 
# /usr/local/Cellar/cantera/2.2.0/share/cantera/data

####   nasa.cti
# open it up.  It has a lot of species with NASA 9 polynomials up to 6000 K temperatures 
# for many species (5727 C or 10340 F!)

######################################################################
####  The problem with h2o2.cti (and GRI30)
######################################################################
# open up h2o2.cti 
# https://github.com/Cantera/cantera/blob/master/data/inputs/h2o2.cti
# notice that it only goes up to 3500 K for many species (3227 C, 5840 F)

#### Solution; grab the species from nasa.cti, and then 
#### pick out the ones you want from another file, h2o2.cti in this case,
#### and grab the reactions from this other file, h2o2.cti in this case

nasaspec = ct.Species.listFromFile('nasa.cti')
H2O2rxns = ct.Reaction.listFromFile('h2o2.cti')

# But we only want the species involved in H2O2 reactions.  Do this:
nasaH2O2spec = [spec for spec in nasaspec if spec.name in [spec1.name for spec1 in H2O2spec] ]

# Look at the list nasaH2O2spec.  Oh no, Argon wasn't included from nasa.cti!
# This is because the name in h2o2.cti is 'AR'.  So this is either a problem with the naming 
# convention in GRI30, being case-insensitive, because I think it really should be "Ar"
# so add this in manually
nasaH2O2spec = nasaH2O2spec + [spec for spec in nasaspec if spec.name=='Ar']
# generalize this into a function.

def parse_spec(desiredspec, srcspec):
    """
    parse_spec
    parse species for a desired list of species, desiredspec, 
                  from a source list of species, srcspec

    INPUTS/PARAMETERS:
    desiredspec - <list of Cantera Species>
    srcspec     - <list of Cantera Species>


    EXAMPLE of USAGE:
    parsed_H2O2spec = parse_spec( H2O2spec, nasaspec )

    """
    desirednames   = [spec.name for spec in desiredspec]
    n_desirednames = len(desirednames)
    desiredspecfromsrc = [spec for spec in srcspec if spec.name in desirednames]
    if len(desiredspecfromsrc) == n_desirednames:
        return desiredspecfromsrc
    desirednames = Set(desirednames)
    srcnames = Set([spec.name for spec in desiredspecfromsrc])
    missingnames = desirednames - srcnames
    # try this: matching lower cases
    loweredmissingnames = Set([name.lower() for name in missingnames])
    missingspecfromsrc = [spec for spec in srcspec if spec.name.lower() in loweredmissingnames]
    return desiredspecfromsrc + missingspecfromsrc
        
# If one runs this
# H2O2gas = ct.Solution(thermo="IdealGas",kinetics="GasKinetics",species=nasaH2O2spec,reactions=H2O2rxns)
# At this point, I obtained an error, again due to capitalization, uppercase vs. lowercase of Ar or "AR" in the reactions.  So for right now, I manually editted the reactions and placed the edits in 
# h2o2_hiT.cti

def isen_flow(gas,p):
    """
    isen_flow = isenflow(gas)
    """
    gas.SP = None,p
    gas.equilibrate('SP')
    return gas

def equil_flow(TP_0,Y_0str,gas,p_c,NDeltap_c):
    """
    equil_flow = equil_flow(TP_0,Y_0str,gas,p_c,NDeltap_c)
    """
    ps = np.linspace(p_c,0,NDeltap_c)[1:-1]
    p_h_rho_a_X = []
    for p in ps:
        gas.Y  = Y_0str
        gas.TP = TP_0
        gas.equilibrate('HP')
#        h_0 = gas.enthalpy_mass
        gas = isen_flow(gas,p)
        p_h_rho_a_X.append( (p,gas.enthalpy_mass, gas.density, equilSoundSpeeds(gas), gas.X,gas.T) )
    return p_h_rho_a_X

def frozen_flow(TP_0,Y_0str,gas,p_c,NDeltap_c):
    """
    frozen_flow = frozen_flow(TP_0,Y_0str,gas,p_c,NDeltap_c)
    """
    ps = np.linspace(p_c,0,NDeltap_c)[1:-1]
    p_h_rho_a_X = [] # store pressure p, enthalpy h, density rho, speed of sound(s) a, mole fraction X
    for p in ps:
        # combustion (after inlets, in start of combustion chamber)
        gas.Y = Y_0str
        gas.TP = TP_0
        gas.equilibrate('HP')
        gas.SP = None,p
        p_h_rho_a_X.append((p,gas.enthalpy_mass,gas.density,equilSoundSpeeds(gas),gas.X,gas.T))
    return p_h_rho_a_X


try:
    H2O2Arrxns = ct.Reaction.listFromFile('h2o2_hiT.cti')
    H2O2gas = ct.Solution(thermo="IdealGas",kinetics="GasKinetics",species=nasaH2O2spec,reactions=H2O2Arrxns)
    H2O2gas.Y = 'O2:6,H2:1'
    H2O2gas.TP = 200., 2*10**7

    H2O2gas.equilibrate('HP')
    h_0 = H2O2gas.enthalpy_mass
    T_C, p_C = H2O2gas.TP

    H2O2p_h_rho_a_X = equil_flow((200.,2*10**7),'O2:6,H2:1',H2O2gas,2*10**7,50000)
    H2O2u = np.sqrt( 2*(h_0 - np.array([row[1] for row in H2O2p_h_rho_a_X])) )
    H2O2rho = np.array([row[2] for row in H2O2p_h_rho_a_X])
    H2O2rhostar = (H2O2rho*H2O2u)[ np.where( (H2O2rho*H2O2u)==max(H2O2rho*H2O2u))]
    H2O2AstarAratio = H2O2rhostar/(H2O2rho*H2O2u)
    equil_exh_index = np.abs( H2O2AstarAratio-77.5).argmin()
    
    equil_u_exh   = H2O2u[equil_exh_index]
    equil_p_exh   = H2O2p_h_rho_a_X[equil_exh_index][0]
    equil_rho_exh = H2O2p_h_rho_a_X[equil_exh_index][2] 

    equil_Isp_vac = (equil_u_exh + equil_p_exh/(equil_rho_exh*equil_u_exh))/9.8
    equil_thrust_vac = H2O2rhostar*(np.pi*(0.2616**2/4))*9.8*equil_Isp_vac

    equil_Isp_sea = (equil_u_exh + (equil_p_exh-ct.one_atm)/(equil_rho_exh*equil_u_exh))/9.8
    equil_thrust_sea = H2O2rhostar*(np.pi*(0.2616**2/4))*9.8*equil_Isp_sea
    equil_dotm = H2O2rhostar*(np.pi*(0.2616**2/4))

    # frozen flow
    H2O2p_h_rho_a_X_frozen = frozen_flow((200.,2*10**7),'O2:6,H2:1',H2O2gas,2*10**7,50000)
    frozen_results = {}
    frozen_results['u'] = np.sqrt(2*(h_0-np.array([row[1] for row in H2O2p_h_rho_a_X_frozen])))
    frozen_results['rho'] = np.array([row[2] for row in H2O2p_h_rho_a_X_frozen])
    frozenrhoustar = (frozen_results['rho']*frozen_results['u'])[ np.where( (frozen_results['rho']*frozen_results['u'])==max(frozen_results['rho']*frozen_results['u']))]
    frozenAstarAratio = frozenrhoustar/(frozen_results['rho']*frozen_results['u'])
    frozen_exh_index = np.abs( frozenAstarAratio-77.5).argmin()
    
    frozen_u_exh   = frozen_results['u'][frozen_exh_index]
    frozen_p_exh   = H2O2p_h_rho_a_X_frozen[frozen_exh_index][0]
    frozen_rho_exh = H2O2p_h_rho_a_X_frozen[frozen_exh_index][2] 
                                                                  
    frozen_Isp_vac = (frozen_u_exh + frozen_p_exh/(frozen_rho_exh*frozen_u_exh))/9.8
    frozen_thrust_vac = frozenrhoustar*(np.pi*(0.2616**2/4))*9.8*frozen_Isp_vac

    frozen_Isp_sea = (frozen_u_exh + (frozen_p_exh-ct.one_atm)/(frozen_rho_exh*frozen_u_exh))/9.8
    frozen_thrust_sea = frozenrhoustar*(np.pi*(0.2616**2/4))*9.8*frozen_Isp_sea
    frozen_dotm = frozenrhoustar*(np.pi*(0.2616**2/4))



except:
    print "h2o2_hiT.cti not here!"

