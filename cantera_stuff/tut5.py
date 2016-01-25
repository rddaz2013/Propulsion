## tut5.py
## tut5.m implemented in Python for cantera 
## cf. http://www.cantera.org/docs/sphinx/html/matlab/tutorials/tut5.html
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
# Tutorial 5:   Reaction information and rates
#
#    Topics:
#       - stoichiometric coefficients
#       - reaction rates of progress
#       - species production rates
#       - reaction equations
#       - equilibrium constants
#       - rate multipliers
#
import cantera as ct
import numpy as np 
import matplotlib.pyplot as plt
g = ct.Solution("gri30")
g.TP = 1500, ct.one_atm
g.X = np.ones(g.n_species)

# Methods are provided that compute many quantities of interest for
# kinetics. Some of these are:


# 1) Stoichiometric coefficients


nu_r   = g.reactant_stoich_coeffs()   # reactant stoichiometric coefficient matrix 
nu_p   = g.product_stoich_coeffs()    # product stoichiometric coefficient matrix 
nu_net = nu_p-nu_r                    # net (product - reactant) stoichiometric
                                      # coefficient matrix 

# cf. https://code.google.com/p/cantera/source/browse/cantera/trunk/interfaces/matlab/toolbox/@Kinetics/stoich_net.m?r=1314&spec=svn1726

# For any of these, the (k,i) matrix element is the stoichiometric
# coefficient of species k in reaction i. These coefficient
# matrices are implemented as Python numpy arrays

# 2) Reaction rates of progress

# Methods forward_rates_of_progress, reverse_rates_of_progress, and 
# net_rates_of_progress return column vectors containing
# the forward, reverse, and net (forward - reverse) rates of
# progress, respectively, for all reactions.

qf = g.forward_rates_of_progress
qr = g.reverse_rates_of_progress
qn = g.net_rates_of_progress
rop = [qf,qr,qn]

# cf. https://code.google.com/p/cantera/source/browse/cantera/trunk/interfaces/matlab/toolbox/@Kinetics/rop_f.m?r=3050&spec=svn3050

# This plots the rates of progress
plt.figure(1)
plt.bar(qf,label='forward')
plt.bar(qr,label='reverse')
plt.bar(qn,label='net')
plt.legend()
plt.show()

# 3) Species production rates

# Methods creationRates, destructionRates, and netProdRates return
# column vectors containing the creation, destruction, and net
# production (creation - destruction) rates, respectively, for all species.

cdot = g.creation_rates
ddot = g.destruction_rates
wdot = g.net_production_rates
rates = [cdot, ddot, wdot]

# This plots the production rates
plt.figure(2)
plt.bar(cdot,label='creation')
plt.bar(ddot,label='destruction')
plt.bar(wdot,label='net')
plt.legend()
plt.show()

# For comparison, the production rates may also be computed
# directly from the rates of progress and stoichiometric
# coefficients. 

cdot2 = np.dot(nu_p,qf) + np.dot(nu_r,qr)
creation = [cdot, cdot2, cdot - cdot2]

ddot2 = np.dot(nu_r,qf) + np.dot(nu_p,qr)
destruction = [ddot, ddot2, ddot - ddot2]

wdot2 = np.dot( nu_net,qn)
net = [wdot, wdot2, wdot - wdot2]

# 4) Reaction equations

e8    = g.reaction_equation(8)          # equation for reaction 8
e1_10 = g.reaction_equations(range(10)) # equation for rxns 1 - 10
eqs   = g.reaction_equations()          # all equations

# 5) Equilibrium constants

# The equilibrium constants are computed in concentration units,
# with concentrations in kmol/m^3.

kc = g.equilibrium_constants
for i in range( g.n_reactions):
    print "%50s  %13.5g" % eqs[i], kc[i]


# 6) Multipliers
 
# For each reaction, a multiplier may be specified that is applied
# to the forward rate coefficient. By default, the multiplier is
# 1.0 for all reactions.

for i in range(g.n_reactions):
    g.set_multiplier(2*i,i)
    m = g.multiplier(i)

# EY : 20160125 Notice that the order of the arguments in the method 
#  set_multiplier in Python is REVERSE of the Matlab setMultiplier
# cf. https://code.google.com/p/cantera/source/browse/cantera/trunk/interfaces/matlab/toolbox/@Kinetics/setMultiplier.m?r=1314&spec=svn1726

#################################################################
#   end of tutorial 5
#################################################################

