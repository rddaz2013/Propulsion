## tut7.py
## tut7.m implemented in Python for cantera 
## cf. http://www.cantera.org/docs/sphinx/html/matlab/tutorials/tut7.html
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160125
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
#  Tutorial 7: Thermodynamic Properties
#
import cantera as ct

# A variety of thermodynamic property methods are provided.

gas = ct.Solution('air.cti') # or gas = ct.Solution('air.cti','air')

gas.TP = 800, ct.one_atm


# temperature, pressure, density
T   = gas.T
P   = gas.P
rho = gas.density
n   = gas.density_mole


# species non-dimensional properties
hrt = gas.standard_enthalpies_RT            # vector of h_k/RT 


# mixture properties per mole
hmole = gas.enthalpy_mole
umole = gas.int_energy_mole
smole = gas.entropy_mole
gmole = gas.gibbs_mole


# mixture properties per unit mass
hmass = gas.enthalpy_mass
umass = gas.int_energy_mass
smass = gas.entropy_mass
gmass = gas.gibbs_mass

#################################################################
