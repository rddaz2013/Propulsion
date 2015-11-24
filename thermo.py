## thermo.py
## I implement symbolic computation in thermodynamics to solve examples and problems
##   
############################################################################ 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
##                                                            
## 20151019
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
## You can have received a copy of the GNU General Public License             
## along with this program; if not, write to the Free Software Foundation, Inc.,  
## S1 Franklin Street, Fifth Floor, Boston, MA                      
## 02110-1301, USA                                                             
##                                                                  
## Governing the ethics of using this program, I default to the Caltech Honor Code: 
## ``No member of the Caltech community shall take unfair advantage of        
## any other member of the Caltech community.''                               
##                                                         
## Donate, and support my other scientific and engineering endeavors at 
## ernestyalumni.tilt.com                                                      
##                          
## Facebook     : ernestyalumni                                                   
## linkedin     : ernestyalumni                                                    
## Tilt/Open    : ernestyalumni                                                    
## twitter      : ernestyalumni                                                   
## youtube      : ernestyalumni                                                   
## wordpress    : ernestyalumni                                                    
##                                                                                  
############################################################################ 

import sympy
from sympy import Symbol, Eq, N
from sympy.solvers import solve
from sympy.abc import gamma
from sympy import Rational as Rat

import decimal 
from decimal import Decimal

import Physique
from Physique import FCconv, KCconv, FundConst, conv, T_C, T_K, T_F

#########################
##### Physical Constants
#########################

N_A = FundConst[ FundConst["Quantity"].str.contains("Avogadro")].loc[42,:]
k_BOLTZ = FundConst[ FundConst["Quantity"].str.contains("Boltzmann") ].loc[49,:]
P_frmATM = conv[ conv["Toconvertfrom"].str.contains("atm") ].loc[15,:] # 1 atm to Pascal conversion for pressure

###########################################################################
##### Kittel, Kroemer. Thermal Physics
###########################################################################

##################################################
####   Chapter 10: Phase Transformation
##################################################

#############################################
##       2. Calculation of $dT/dp$ for water
#############################################

L = 2260 # J g^{-1}
boilingwatertemp_K = KCconv.subs(T_C,100).rhs # room temperature in Kelvin
tau_boilingwater = boilingwatertemp_K * k_BOLTZ.Value
P = 1 # atm

dpdtau_1002 = L*18.0153/float(N_A.Value)/( tau_boilingwater**2 )
dtaudp_1002 = 1./ dpdtau_1002 # in J/Pascal
dTdp_1002 = dtaudp_1002 / (k_BOLTZ.Value )  # 28.4348535111262 K/atm

#############################################
##       3. Heat of vaporization of ice
#############################################

L1003 = ( 4.58 - 3.88 )/(0 - (-2.) )*((KCconv.subs(T_C,1.).rhs)**2)/( (4.58-3.88)/(0 - (-2)) * ( 1. ) + 3.88 ) * k_BOLTZ.Value*N_A.Value  # 51705.6757640485


###########################################################################
##### Ralph Baierlein, Thermal Physics, Cambridge University Press, 1999
###########################################################################

######################################## 
####   1 Background
###     Problems
########################################

###################################
##       4. Adiabatic compression
###################################

p_i   = Symbol("p_i", positive=True)
V_i   = Symbol("V_i", positive=True)
tau_i = Symbol("tau_i", positive=True)
N_i   = Symbol("N_i", positive=True)
idealgaslaw_i = Eq( p_i*V_i, N_i*tau_i)

p_f   = Symbol("p_f", positive=True)
V_f   = Symbol("V_f", positive=True)
tau_f = Symbol("tau_f", positive=True)
idealgaslaw_f = Eq( p_f*V_f, N_i*tau_f)

adia_tV = Eq( tau_i*V_i**(gamma-1) , tau_f*V_f**(gamma-1) )

##############################
#         (a)
##############################

roomtemp_K = KCconv.subs(T_C,20).rhs # room temperature in Kelvin

Prob0104ans = adia_tV.subs(gamma,1.4).subs(V_f,1).subs(V_i,15).subs(tau_i, roomtemp_K) # answer to Problem 4 of Chapter 1

Prob0104ans = N( Prob0104ans.lhs) # 866.016969686253 K 
Prob0104ansC = solve( KCconv.subs( T_K, Prob0104ans), T_C )[0] # 592.866969686253 C 
solve( FCconv.subs( T_C, Prob0104ansC ), T_F)[0] # 1099.16054543526 F

##############################
#         (b)
##############################

15*( Prob0104ans / roomtemp_K )
