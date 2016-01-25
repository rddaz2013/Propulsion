## tut4.py
## tut4.m implemented in Python for cantera 
## cf. http://www.cantera.org/docs/sphinx/html/matlab/tutorials/tut4.html
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
# Tutorial 4:   Chemical Equilibrium
#
#   Topics:
#     - the equilibrate method
#     - specifying fixed TP, HP, UV, SV, or SP
#     - checking reaction rates of progress
#
import cantera as ct
# To set a gas mixture to a state of chemical equilibrium, use the
# 'equilibrate' method.
#
g = ct.Solution('gri30.xml')
g.TP = 1200.0, ct.one_atm
g.X = 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('TP')

# The above statement sets the state of object 'g' to the state of
# chemical equilibrium holding temperature and pressure
# fixed. Alternatively, the specific enthalpy and pressure can be held
# fixed:

print "fixed H and P:"
g.TP = 1200.0, ct.one_atm
g.X = 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('HP')

# Other options are
#     'UV'   fixed specific internal energy and specific volume
#     'SV'   fixed specific entropy and specific volume
#     'SP'   fixed specific entropy and pressure

print "fixed U and V:"
g.TP = 1200.0, ct.one_atm
g.X = 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('UV')

print "fixed S and V:"
g.TP = 1200.0, ct.one_atm
g.X = 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('SV')

print "fixed S and P:"
g.TP = 1200.0, ct.one_atm
g.X = 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('SP')

# How can you tell if 'equilibrate' has correctly found the
# chemical equilibrium state? One way is verify that the net rates of
# progress of all reversible reactions are zero.

# Here is the code to do this:

g.TP = 2000.0, ct.one_atm
g.X = 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('TP')
rf = g.forward_rates_of_progress
rr = g.reverse_rates_of_progress

for i in range(g.n_reactions):
    if g.is_reversible(i):
        print i, rf[i], rr[i], (rf[i]-rr[i])/rf[i]

# EY : 20160125 see the function check_equil in equil.py and LOXmeth_eq.py 
# for its implementation

# You might be wondering how 'equilibrate' works. (Then again, you might
# not, in which case you can go on to the next tutorial now.)  Method
# 'equilibrate' invokes Cantera's chemical equilibrium solver, which
# uses an element potential method. The element potential method is
# one of a class of equivalent 'nonstoichiometric' methods that all
# have the characteristic that the problem reduces to solving a set of
# M nonlinear algebraic equations, where M is the number of elements
# (not species). The so-called 'stoichiometric' methods, on the other
# hand, (including Gibbs minimization), require solving K nonlinear
# equations, where K is the number of species (usually K >> M). See
# Smith and Missen, "Chemical Reaction Equilibrium Analysis" for more
# information on the various algorithms and their characteristics.
#
# Cantera uses a damped Newton method to solve these equations, and
# does a few other things to generate a good starting guess and to
# produce a reasonably robust algorithm. If you want to know more
# about the details, look at the on-line documented source code of
# Cantera C++ class 'ChemEquil' at http://www.cantera.org. 

#################################################################
#   end of tutorial 4
#################################################################


