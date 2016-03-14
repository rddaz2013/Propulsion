## CC_params.py
## Combustion chamber parameters; the (numerical) physical parameters 
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
from collections import namedtuple

Chamber_Params = namedtuple('Chamber_Parameters',['A_totinj','A_cc','l_cc'])

Inlet_Conds    = namedtuple('Inlet_Conditions',['TP','phi_g','phi_overall','v_d','D'])

Reaction_Tuple = namedtuple('Reaction_Tuple',['Fuelstr','Oxstr','n_F','n_Ox','Prod_str'])


##########################################################################################
##### (Physical) Chamber Parameters for a rocket using CH4 (methane) as a fuel
##########################################################################################
CH4CH_PARAMS = Chamber_Params(A_totinj= 0.0157, # m^2 total fuel injector cross-sectional area  
                              A_cc    = 0.157,  # m^2 combustion chamber cross-sectional area
                              l_cc    = 0.75)   # m length of combustion chamber

##########################################################################################
##### (Physical) Chamber Parameters for a rocket using C7H16 (n-heptane) as a fuel
##########################################################################################
HEPCH_PARAMS = Chamber_Params(A_totinj= 0.0157, # m^2 total fuel injector cross-sectional area  
                              A_cc    = 0.157,  # m^2 combustion chamber cross-sectional area
                              l_cc    = 0.725)  # m length of combustion chamber


##########################################################################################
##### Inlet Conditions for a rocket using CH4 (methane) as a fuel
##########################################################################################
D_IDROPLET  = (30,50,80,100,200)       # \mu m (microns) initial (injected) droplet diameter 

# Inlet Conditions that vary with size of droplet 
CH4INLETS_VARY_D_0 = [Inlet_Conds(
        TP=(600,                       # K inlet gas temperature 
            3.4474*10**6),            # Pa Combustion Chamber Pressure
# Assume oxidizer injected as gas and fuel partly as gas with initial equivalence ratio 0.45
        phi_g      = 0.45,             # initial equivalence ratio
        phi_overall= 1.139,            # overall equivalence ratio 
        v_d        = 10,               # m/s initial droplet velocity
        D=D_i*10**(-6)) for D_i in D_IDROPLET]


##########################################################################################
##### Inlet Conditions for a rocket using C7H16 (n-heptane) as a fuel
##########################################################################################
HEPINLETS_VARY_D_0 = [Inlet_Conds(
        TP=(801,                       # K inlet gas temperature 
            3.4474*10**6),            # Pa Combustion Chamber Pressure
# Assume oxidizer injected as gas and fuel partly as gas with initial equivalence ratio 0.45
        phi_g      = 0.45,             # initial equivalence ratio
        phi_overall= 2.3,              # overall equivalence ratio 
        v_d        = 10,               # m/s initial droplet velocity
        D=D_i*10**(-6)) for D_i in D_IDROPLET]


##########################################################################################
##### Chemistry for CH4 (methane) as a fuel, O2 as an oxidizer
##########################################################################################
# Suppose stoichiometric reaction for LOX-methane is
# CH4 + 2O2 -> CO2+2H2O
CH4_TUP   = Reaction_Tuple(Fuelstr='CH4',Oxstr='O2',n_F=1,n_Ox=2,Prod_str='CO2:1 H2O:2')

##########################################################################################
##### Chemistry for CH4 (methane) as a fuel, O2 as an oxidizer
##########################################################################################
# Suppose stoichiometric reaction for LOX-n-heptane is
# C7H16 + 11O2 -> 7CO2+8H2O

C7H16_TUP = Reaction_Tuple(Fuelstr='NC7H16',Oxstr='O2',n_F=1,n_Ox=11,Prod_str='CO2:7 H2O:8')
