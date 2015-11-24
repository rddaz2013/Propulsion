## NozzleTheory.py
## Nozzle theory for rocket propulsion implemented in sympy
###################################################################################### 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20151119                                                                              
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
## Donate, and support my other scientific and engineering endeavors at 
## ernestyalumni.tilt.com                                                      
###################################################################################### 
import decimal
from decimal import Decimal

import sympy
from sympy import *
from sympy.abc import a, A, p, R, t, u, gamma, psi, rho, theta
from sympy import Rational as Rat

from sympy.utilities.lambdify import lambdify, implemented_function

import numpy as np

import scipy 
import scipy.optimize
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

import pylab

import Physique
from Physique import FCconv, KCconv, FundConst, conv, plnfacts, T_C, T_K, T_F

# Physical Constants
k_Boltz = FundConst[ FundConst["Quantity"].str.contains("Boltzmann") ].loc[49,:]
N_Avog = FundConst[FundConst["Quantity"].str.contains("Avogadro") ].loc[42,:]
g_0std = FundConst[FundConst['Quantity'].str.contains("gravity")].loc[303:,]

# Isentropic relations of ideal gas

p = Symbol('p',real=True)
p_i = Symbol('p_i',real=True)
tau = Symbol('tau',real=True)
tau_i = Symbol('tau_i',real=True)
rho = Symbol('rho',positive=True)
rho_i = Symbol('rho_i',positive=True)
gamma = Symbol('gamma',positive=True)
V = Symbol('V',positive=True)
V_i = Symbol('V_i',positive=True)

isen_pvstau = Eq( p/p_i , ( tau/tau_i)**( gamma / (gamma-1) ) )
isen_tauvsV = Eq( tau/tau_i , (V_i / V)**(gamma -1) )
isen_pvsV   = Eq( p/p_i , (V_i/V)**gamma)
N_0 = Symbol('N_0',positive=True)

# Alternatively, you can make substitutions of the ideal gas law in as you'd like
idealgaslaw   = Eq( p*V, N_0*tau)
idealgaslaw_i = Eq( p_i*V_i, N_0*tau_i)


# exhaust velocity cf. Sutton and Biblarz pp. 52 3.3. Isentropic Flow Through Nozzles, Velocity, (3-15b)
p_0 = Symbol('p_0', positive=True)
T_0 = Symbol('T_0',real=True) # stagnation temperature
v_1 = Symbol('v_1',real=True)
velocity_along_flow = sqrt( 2*gamma*R*T_0/(gamma-1)*(1- (p/p_0)**((gamma-1)/gamma) ) + v_1**2)

# stagnation temperature to temp ratio to Mach number
tau_0 = Symbol('tau_0',real=True)
Mach = Symbol('Mach',positive=True)
tempstoMach = Eq( tau_0/tau,Rat(1) + (gamma-Rat(1))/Rat(2)*Mach**2)

# Area Ratio to Mach numbers
Mach_1 = Symbol('Mach_1',positive=True)
Mach_2 = Symbol('Mach_2',positive=True)
A_1    = Symbol('A_1', positive=True)
A_2    = Symbol('A_2', positive=True)
AreastoMachs = Eq( A_2/A_1 , Mach_1/Mach_2*( ( (Rat(1) + (gamma-Rat(1))/(Rat(2) )*Mach_2**2 )/(Rat(1) + (gamma-Rat(1))/(Rat(2) )*Mach_1**2 ) )**( (gamma +1)/(gamma-1) ) )**0.5 )

massflow = Symbol('massflow',real=True)
Astar = Symbol('A^*',positive=True)
T_0 = Symbol('T_0',real=True) # stagnation temperature

p_0 = Symbol('p_0', positive=True)
massflowrateExp = Eq( massflow/Astar,(p_0*sqrt(gamma))/(sqrt(R*T_0))*( Rat(2)/( gamma +1 ) )**( (gamma +1)/(Rat(2)*(gamma -1) ) ) )

# k_Boltz.Value/ (Decimal(23)/N_Avog.Value )

# Thrust
p_a = Symbol('p_a',real=True)
Thrust = massflowrateExp.rhs*Astar*velocity_along_flow +(p-p_a)*A


#######################################################
##### AE 121
#######################################################

#######################################################
####   PS 5
#######################################################

#######################################################
###     Problem 1: Nozzle flow in liquid rocket engines
#######################################################


Viking5Cnozzle = massflowrateExp.subs(gamma, 1.2).subs(massflow, 275.2).subs(p_0,5800*1000).subs(T_0,3350).subs(R, k_Boltz.Value/ (Decimal(23/1000.)/N_Avog.Value ))
Viking4Bnozzle = massflowrateExp.subs(gamma, 1.2).subs(massflow, 278.0).subs(p_0,5850*1000).subs(T_0,3350).subs(R, k_Boltz.Value/ (Decimal(23/1000.)/N_Avog.Value ))

Viking5CAstar = solve(Viking5Cnozzle, Astar)[0]
Viking4BAstar = solve(Viking4Bnozzle, Astar)[0]

# part(b)

# cf. https://sites.google.com/a/aims-senegal.org/scipy/roots-finding-numerical-integrations-and-differential-equations

Viking5CMachEq = AreastoMachs.subs(Mach_1,Rat(1) ).subs(gamma,1.2).subs(A_2,A_1*10.5) 

Viking5CMach = lambdify(Mach_2, Viking5CMachEq.rhs - Viking5CMachEq.lhs ) # Remember to move all terms to 1 side, and so the other side equals 0
# cf. http://docs.sympy.org/dev/modules/utilities/lambdify.html
# cf. http://stackoverflow.com/questions/22742951/solve-an-equation-using-a-python-numerical-solver-in-numpy

# cf. https://sites.google.com/a/aims-senegal.org/scipy/roots-finding-numerical-integrations-and-differential-equations

scipy.optimize.newton( Viking5CMach, 3) # Newton Raphson method
# 3.3123573073570207

scipy.optimize.bisect( Viking5CMach, 3,4)  # Bisection method
# 3.312357307356251

Viking4BMachEq = AreastoMachs.subs(Mach_1,Rat(1) ).subs(gamma,1.2).subs(A_2,A_1*30.8) 
Viking4BMach   = lambdify(Mach_2, Viking4BMachEq.rhs - Viking4BMachEq.lhs )

plot(Viking4BMachEq.rhs, Mach_2) # plotting helps to visually estimate what the solution is
plot(Viking4BMachEq.rhs, (Mach_2,3.5,4.5)) # 

scipy.optimize.newton( Viking4BMach, 3.5) # 4.0573823996942675
scipy.optimize.bisect( Viking4BMach, 3.5,4.5) # 4.057382399693779

tempstoMach.subs(gamma,1.2).subs(tau_0,3350)

# static temperature at nozzle exit
solve(tempstoMach.subs(gamma,1.2).subs(tau_0,3350).subs(Mach,scipy.optimize.newton(Viking5CMach,3)),tau) # [1597.38993681819]
solve(tempstoMach.subs(gamma,1.2).subs(tau_0,3350).subs(Mach,scipy.optimize.newton(Viking4BMach,3)),tau) # [1265.94945450477]

Viking5CTrat = tempstoMach.subs(gamma,1.2).subs(Mach,scipy.optimize.newton(Viking5CMach,3)).rhs # Viking 5C temperature ratio
Viking4BTrat = tempstoMach.subs(gamma,1.2).subs(Mach,scipy.optimize.newton(Viking4BMach,3.5)).rhs # Viking 5C temperature ratio

# static pressure at nozzle exit
Viking5Cp_e = solve(isen_pvstau.subs(gamma,1.2).subs(tau_i,Viking5CTrat*tau).subs(p_i,5800*10**3),p)[0] # [68174.9481403000]
Viking4Bp_e = solve(isen_pvstau.subs(gamma,1.2).subs(tau_i,Viking4BTrat*tau).subs(p_i,5850*10**3),p)[0] # [17036.6899478252]

Viking5Crho_0 = Decimal(5800*10**3)/(k_Boltz.Value*N_Avog.Value/Decimal(23*10**(-3))*Decimal(3350))
Viking4Brho_0 = Decimal(5850*10**3)/(k_Boltz.Value*N_Avog.Value/Decimal(23*10**(-3))*Decimal(3350))

# static density at nozzle exit

Viking5CrhoEq = isen_pvsV.subs(gamma,1.2).subs( p,p_i*isen_pvstau.subs(gamma,1.2).subs(tau_i,Viking5CTrat*tau).rhs).subs(V,rho_i).subs(V_i,rho).subs(rho_i,Viking5Crho_0)
# 0.0117543014035 == rho**1.2*rho_i**(-1.2)

Viking4BrhoEq = isen_pvsV.subs(gamma,1.2).subs( p,p_i*isen_pvstau.subs(gamma,1.2).subs(tau_i,Viking4BTrat*tau).rhs).subs(V,rho_i).subs(V_i,rho).subs(rho_i,Viking4Brho_0)
# 0.00291225469193595 == rho**1.2*rho_i**(-1.2)

# plot(Viking5Crho.rhs,(rho,0,0.2))
# plot(Viking4Brho.rhs,(rho,0,0.1))

Viking5Crho = lambdify(rho,Viking5CrhoEq.rhs-Viking5CrhoEq.lhs)
Viking4Brho = lambdify(rho,Viking4BrhoEq.rhs-Viking4BrhoEq.lhs)

scipy.optimize.newton(Viking5Crho,0) # Newton Raphson method
scipy.optimize.bisect(Viking5Crho,0,0.2)  # Bisection method

scipy.optimize.newton(Viking4Brho,0) # Newton Raphson method
scipy.optimize.bisect(Viking4Brho,0,0.1)  # Bisection method

# part (c)

Viking5Cthrust = Thrust.subs(gamma,1.2).subs(T_0,3350.).subs(v_1,0).subs(p_0,5800.*10**3).subs(R,k_Boltz.Value*N_Avog.Value/Decimal(23*10**(-3))).subs(p,Viking5Cp_e).subs(Astar,Viking5CAstar).subs(A,Viking5CAstar*10.5)
Viking4Bthrust = Thrust.subs(gamma,1.2).subs(T_0,3350.).subs(v_1,0).subs(p_0,5850.*10**3).subs(R,k_Boltz.Value*N_Avog.Value/Decimal(23*10**(-3))).subs(p,Viking4Bp_e).subs(Astar,Viking4BAstar).subs(A,Viking4BAstar*30.8)

#Viking5Cthrust = lambdify(p_a, Viking5Cthrust)
#Viking4Bthrust = lambdify(p_a, Viking4Bthrust)

# U.S. Standard Atmosphere data 

# 1. Manually inputting in data
Manual_US_std_atm = {"h_alt":np.array([0.,1.,2.,5.,10.,12.,15.,20.,25.,30.,40.,50.,60.,70.,80.,90.,100.]),
                     "p_mbar":np.array([1013,899,795,540,265,194,121,55,25,12,3,0.8,0.2,0.06,0.01,0.002,0.0003])}

# cf. http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.curve_fit.html
# cf. http://www.pdas.com/hydro.pdf
# cf. https://engineering.purdue.edu/~andrisan/Courses/AAE490A_S2002/Atmosphere.pdf
def US_std_atm_pvsH_exp_func(h,p_b,h_0,h_1):
    return p_b*np.exp(-h_0*h+h_1)

def US_std_atm_lnpvsH_poly_func(h,a,b,c,d,e):
    return (a*h+b*h**2+c*h**3+d*h**4+e)

def US_std_atm_lnpvsH_ln_func(h,lnp_b,h_0,g_0):
    return ( lnp_b + g_0*np.log(1 - h_0*h ) )

lnpvsHpopt,lnpvsHpcov = curve_fit(US_std_atm_lnpvsH_poly_func, Manual_US_std_atm['h_alt']*1000, np.log(Manual_US_std_atm['p_mbar']*100 ))

Height = Symbol('Height',real=True)

US_std_atm_pvsH_polyExp = exp(lnpvsHpopt[0]*Height+lnpvsHpopt[1]*Height**2+lnpvsHpopt[2]*Height**3+lnpvsHpopt[3]*Height**4+lnpvsHpopt[4])

Viking5Cthrust_lambda = lambdify(Height,Viking5Cthrust.subs(p_a,US_std_atm_pvsH_polyExp))
Viking4Bthrust_lambda = lambdify(Height,Viking4Bthrust.subs(p_a,US_std_atm_pvsH_polyExp))

Heights = pylab.linspace(0.,100000.,10000)

# cf. http://stackoverflow.com/questions/8079061/function-application-over-numpys-matrix-row-column
plt.figure(1)
plt.subplot(211)
plt.plot(Heights, np.vectorize(Viking5Cthrust_lambda)(Heights) )
plt.subplot(212)
plt.plot(Heights, np.vectorize(Viking5Cthrust_lambda)(Heights)/(float(g_0std.Value)*275.2))

plt.figure(2)
plt.subplot(211)
plt.plot(Heights, np.vectorize(Viking4Bthrust_lambda)(Heights) )
plt.subplot(212)
plt.plot(Heights, np.vectorize(Viking4Bthrust_lambda)(Heights)/(float(g_0std.Value)*278.0))

# part (d)
Viking5Cthrustideal = 275.2*velocity_along_flow.subs(v_1,0).subs(gamma,1.2).subs(T_0,3350.).subs(R,k_Boltz.Value/(Decimal(23/1000.)/N_Avog.Value)).subs(p_0,5800.*10**3).subs(p,US_std_atm_pvsH_polyExp)
Viking4Bthrustideal = 278.0*velocity_along_flow.subs(v_1,0).subs(gamma,1.2).subs(T_0,3350.).subs(R,k_Boltz.Value/(Decimal(23/1000.)/N_Avog.Value)).subs(p_0,5850.*10**3).subs(p,US_std_atm_pvsH_polyExp)

Viking5Cthrustideal_lambda = lambdify(Height,Viking5Cthrustideal)
Viking4Bthrustideal_lambda = lambdify(Height,Viking4Bthrustideal)

plt.figure(3)
plt.subplot(211)
plt.plot(Heights, np.vectorize(Viking5Cthrustideal_lambda)(Heights) )
plt.subplot(212)
plt.plot(Heights, np.vectorize(Viking5Cthrustideal_lambda)(Heights)/(float(g_0std.Value)*275.2))

plt.figure(4)
plt.subplot(211)
plt.plot(Heights, np.vectorize(Viking4Bthrustideal_lambda)(Heights) )
plt.subplot(212)
plt.plot(Heights, np.vectorize(Viking4Bthrustideal_lambda)(Heights)/(float(g_0std.Value)*278.0))

plt.figure(5)
plt.subplot(211)
plt.plot(Heights, np.vectorize(Viking5Cthrust_lambda)(Heights) )
plt.plot(Heights, np.vectorize(Viking5Cthrustideal_lambda)(Heights) )
plt.subplot(212)
plt.plot(Heights, np.vectorize(Viking5Cthrust_lambda)(Heights)/(float(g_0std.Value)*275.2))
plt.plot(Heights, np.vectorize(Viking5Cthrustideal_lambda)(Heights)/(float(g_0std.Value)*275.2))

plt.figure(6)
plt.subplot(211)
plt.plot(Heights, np.vectorize(Viking4Bthrust_lambda)(Heights) )
plt.plot(Heights, np.vectorize(Viking4Bthrustideal_lambda)(Heights) )
plt.subplot(212)
plt.plot(Heights, np.vectorize(Viking4Bthrust_lambda)(Heights)/(float(g_0std.Value)*278.0))
plt.plot(Heights, np.vectorize(Viking4Bthrustideal_lambda)(Heights)/(float(g_0std.Value)*278.0))


#######################################################
###     Problem 3: Duct flow with heating
#######################################################
# part a
(1.*10**6)/(40.*15.6) # 1602.5641025641025

# part b
massconsEq = Eq(massflow,p_0*sqrt(gamma/(R*T_0))*A*Mach*(1+(gamma-Rat(1))/Rat(2)*Mach**2)**((gamma+1)/(-Rat(2)*(gamma-1))))
massconsProb0503 = massconsEq.subs(massflow,40.*10**(-3)).subs(gamma,1.4).subs(p_0,6.8*10**(6)).subs(T_0,673.).subs(A,N(pi)*(10**(-2)/2.)**2).subs(R,k_Boltz.Value/(Decimal(2.0159*10**(-3))/N_Avog.Value))

MachProb0503lamb = lambdify(Mach,massconsProb0503.rhs-massconsProb0503.lhs) # Remember to move all terms to 1 side, and so the other side equals 0
# cf. http://docs.sympy.org/dev/modules/utilities/lambdify.html
# cf. http://stackoverflow.com/questions/22742951/solve-an-equation-using-a-python-numerical-solver-in-numpy
# cf. https://sites.google.com/a/aims-senegal.org/scipy/roots-finding-numerical-integrations-and-differential-equations
plot(massconsProb0503.rhs,(Mach,0,5))
MachProb0503 = scipy.optimize.newton(MachProb0503lamb,0)

T_01=Symbol('T_01',real=True)
heataddTvsMachEq=Eq(T_0/T_01,((Rat(1)+gamma*Mach_1**2)/(Rat(1)+gamma*Mach**2)*Mach/Mach_1)**2*(Rat(1)+(gamma-Rat(1))/Rat(2)*Mach**2)/(1+(gamma-Rat(1))/Rat(2)*Mach_1**2) )

heataddTvsMachProb0503=heataddTvsMachEq.subs(gamma,1.4).subs(T_0,675).subs(T_01,675+(1.*10**6)/(40.*15.6)).subs(Mach,MachProb0503)
Mach1Prob0503lamb=lambdify(Mach_1,heataddTvsMachProb0503.rhs-heataddTvsMachProb0503.lhs)
plot(heataddTvsMachProb0503.rhs,(Mach_1,0,10))
Mach1Prob0503=scipy.optimize.newton(Mach1Prob0503lamb,0.1)

# part (c)

AtoMProb0503=AreastoMachs.subs(A_2,A_1*100).subs(Mach_1,1).subs(gamma,1.4)
# plot(AtoMProb0503.rhs,(Mach_2,0,7))
AtoMProb0503lamb=lambdify(Mach_2,AtoMProb0503.rhs-AtoMProb0503.lhs)
MachProb0503c=scipy.optimize.newton(AtoMProb0503lamb,3)  # 6.936

u_eProb0503cEq=Mach*sqrt(gamma*R*T_0/(Rat(1) + (gamma-Rat(1))/Rat(2)*Mach**2))
u_eProb0503c=u_eProb0503cEq.subs(gamma,1.4).subs(T_0,675+(1.*10**6)/(40.*15.6)).subs(Mach,MachProb0503c).subs(R,k_Boltz.Value/(Decimal(2.0159*10**(-3))/N_Avog.Value))
