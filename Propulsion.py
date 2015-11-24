## Propulsion.py
## Problems and solutions in Propulsion
###################################################################################### 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20151112                                                                              
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
from sympy.abc import a,A, t, u, psi, rho, theta
from sympy import Rational as Rat

from sympy.utilities.lambdify import lambdify, implemented_function

import Physique
from Physique import FCconv, KCconv, FundConst, conv, plnfacts, T_C, T_K, T_F

t_p = Symbol('t_p', real=True) # burn rate
g_0 = Symbol('g_0', positive=True) # standard gravity
F_thrust = Function('F_thrust')(t)
I_t = integrate(F_thrust,(t,0,t_p))  # total impulse
m = Function('m')(t) # mass of propellant flowing out
W_p = Symbol('W_p',positive=True) # weight of propellant
I_sp = I_t/W_p
I_sp.subs( W_p, (g_0*integrate(m.diff(t),(t,0,t_p))) ) # specific impulse

M0 = Symbol('M0',positive=True)
m_p = Symbol('m_p',positive=True)
M = Function('M')(t) # mass of rocket+propellant system

massflow = Symbol('massflow',real=True)
u_e = Symbol('u_e',real=True) # effective exhaust velocity, $c$ for Bibliarz and Sutton

u=Function('u')(t)

M_constantflow = M0 - t*m_p/t_p


# assume constant mass flow
I_t.subs(F_thrust, massflow*u_e).doit() 
I_sp.subs(F_thrust,massflow*u_e).subs(W_p, g_0*massflow*t_p).doit()  # u_e/g_0

# cf. 4.1 Gravity-Free, Drag-Free Space Flight
# Biblarz, Sutton, Rocket Propulsion Elements (2001)

gravityfreedragfreespaceflight = Eq(  u.diff(t), massflow*u_e/M )
gravityfreedragfreespaceflight.subs(M,M_constantflow)
Deltau_g0D0 = integrate( gravityfreedragfreespaceflight.subs(M,M_constantflow).rhs , (t,0,t_p)).simplify() # \Delta u for g_0 = 0, D = 0(gravity-free, drag-free)


# cf. 4.2 Forces acting on a Vehicle in the Atmosphere
# Biblarz, Sutton, Rocket Propulsion Elements (2001)

C_L = Symbol('C_L', positive=True)
C_D = Symbol('C_D', positive=True)
Lift = C_L*(Rat(1)/Rat(2))*rho*A*u**2
Drag = C_D*(Rat(1)/Rat(2))*rho*A*u**2


theta = Function('theta')(t)

flightpathdirection = Eq( u.diff(t), F_thrust/M*cos(psi-theta) - Drag/M - g_0*sin(theta) ) 
tangentialflightdirection = Eq( u*theta.diff(t) , F_thrust/M*sin(psi-theta)+ Lift/M - g_0*cos(theta) )


# Example 4-1 Biblarz, Sutton, Rocket Propulsion Elements (2001)
# assume constant thrust
F_thrust0 = Symbol('F_thrust0',real=True)
I_sp = Symbol('I_sp',positive=True)
I_spEq = Eq(I_sp, I_t/W_p) # specific impulse equation

# Given
# Launch weight                        4.0 lbf
# Useful propellant mass               0.4 lbm
# Effective specific impulse           120 sec
# Launch angle (relative to horizontal 80 degrees
# Burn time (with constant thrust)     1.0 sec
theta_0 = Symbol("theta_0",real=True)

I_spEq.subs(I_t,F_thrust0*t_p) # I_sp == F_thrust0*t_p/W_p
solve( I_spEq.subs(I_t,F_thrust0*t_p).subs(I_sp,120.).subs(t_p,1.0).subs(W_p, 0.4), F_thrust0) # [48.0000000000000] lbf

# "The direction of thrust and the flight path are the same
udot = Matrix([ [flightpathdirection.rhs],[tangentialflightdirection.rhs]])
Rot  = Matrix([[ cos(theta), -sin(theta)],[sin(theta),cos(theta)]])

# assume negligible Drag (low velocity), no lift (wingless)
udot.subs(Lift,0).subs(Drag,0).subs(psi,theta)

( Rot * udot.subs(Lift,0).subs(Drag,0).subs(psi,theta)).expand() # This reproduces the acceleration in x and y components of the powered flight stage

( Rot * udot.subs(Lift,0).subs(Drag,0).subs(psi,theta)).expand().subs(F_thrust, 48.0).subs(M,4.0/32.2).subs(g_0,32.0).subs( theta, 80./180.*N(pi) ) # 67.1 ft/sec^2 in x direction, 348.5 ft/sec^2 in y direction

uxuydot = (Rot*udot.subs(Lift,0).subs(Drag,0).subs(psi,theta)).expand().subs(F_thrust,48.0).subs(g_0,32.0).subs(theta,80./180.*N(pi)).subs(M,M_constantflow).subs(M0,4.0/32.2).subs(m_p,0.4/32.2).subs(t_p,1.0)

u_p = integrate(uxuydot,(t,0,1.0) ) # Matrix([
# [70.6944361984026], 70.7 ft/sec
# [368.928070760125]]) 375 ft/sec
# EY : 20151113 Launch weight is in 4.0 lbf, useful propellant mass was in 0.4 lbm, and yet Biblarz and Sutton divides by 32.2 ft/sec^2 for both, and lbf and lbm are along the same footing as the units for initial weight and final weight on pp. 114; is this wrong?  Answer is no.  lbm is lbf, but in different contexts, see this clear explanation: https://youtu.be/4ePaKh9QyC8

atan( u_p[1]/u_p[0])*180/N(pi) # 79.1524086456152

# Problems Ch. 4 Flight Performance pp. 154
# 3. 
Problem0403 = flightpathdirection.subs(Drag,0).subs(psi,theta).subs(theta,pi/2).subs(F_thrust, m_p/t_p*u_e).subs(M,M_constantflow).subs(m_p,M0*0.57).subs(t_p,5.).subs(u_e,2209.).subs(g_0,9.8).factor(M0).rhs  # Chapter 4 Flight Performance, Problem 3
integrate( Problem0403, (t,0,5.0) ) # 1815.32988528061

integrate( integrate(Problem0403,(t,0,t) ),(t,0,5.0) ) # 3890.37850288891

# Problem 6, Ch. 4 Flight Performance pp. 155

M_earth = plnfacts.loc[plnfacts['Planet']=="EARTH","Mass (1024kg)"].values[0]*10**(24) # in kg
R_earth = plnfacts.loc[plnfacts['Planet']=="EARTH","Diameter (km)"].values[0]/Decimal(2)

Gconst = FundConst[ FundConst["Quantity"].str.contains("gravitation") ].loc[243,"Value"]
v0406 = sqrt( Gconst*M_earth/((R_earth + Decimal(500))*10**3) )  # velocity of satellite v of Chapter 4, Problem 6 of Biblarz and Sutton
T0406 = (2.*N(pi)*float((R_earth + Decimal(500))*10**3 )**(3./2))/float(sqrt( Gconst*M_earth)) 

Eperm0406 = Gconst*M_earth*(-1/(2*((R_earth+Decimal(500))*10**3)) + 1/(R_earth*10**3)) # Energy per mass
Eperm0406.quantize(Decimal('100000.'))
# cf. https://gist.github.com/jackiekazil/6201722
# cf. http://stackoverflow.com/questions/6913532/display-a-decimal-in-scientific-notation
'%.6E' % Eperm0406


##############################
##### AE 121
##############################

#########################
####   PS 2
#########################

####################
###     Problem 1
####################

gstd = FundConst[ FundConst["Quantity"].str.contains("gravity") ].loc[303,:].Value
M_0 = Symbol('M_0',positive=True)
Deltau = -I_sp*g_0*ln( (M_0 -m_p)/M_0) 
# part (a)
Deltau.subs(I_sp,268.8).subs(g_0,gstd).subs(M_0,805309.).subs(m_p, (1-0.1396)*586344) # 2595.74521034101 m/s

# part (b)
Deltau.subs(I_sp,452.1).subs(g_0,gstd).subs(M_0,183952+35013.).subs(m_p, (1-0.1110)*183952) # 6090.68716730318 m/s

# part (c)
1.5*805309./268.8 # 4493.911830357143

#################### 
###     Problem 3 
####################

import scipy
from scipy import exp, array
from scipy.integrate import ode

import matplotlib.pyplot as plt

M_cannonball = (7.8*(10**2)**3/(10**3))*4./3.*N(pi)*(15./2./100.)**3
(1.225)*(0.1)/(2.*M_cannonball)*(N(pi)*(15./2./100.)**2)  # 7.85256410256411e-5

def deriv(t,u): # return derivatives of the array u
    """
    cf. http://bulldog2.redlands.edu/facultyfolder/deweerd/tutorials/Tutorial-ODEs.pdf

    """
    uxdot = (7.853*10**(-5))*exp( -u[3]/(10000.))*(u[0]**2 + u[1]**2)**(0.5)*(-u[0])
    uydot = -9.8 + (7.853*10**(-5))*exp(-u[3]/(10000.))*(u[0]**2 + u[1]**2)**(0.5)*(-u[1])
    return array([ uxdot,uydot,u[0],u[1] ])

u0 = [300.*cos(50./180.*N(pi)), 300.*sin(50./180.*N(pi)),0,0]

Prob0203 = ode(deriv).set_integrator('dopri5')  # Problem 3 from Problem Set 2 for AE121 Fall 2015
# cf. http://stackoverflow.com/questions/26738676/does-scipy-integrate-ode-set-solout-work
Prob0203.set_initial_value(u0)

t1 = 41.575
dt = 0.005
while Prob0203.successful() and Prob0203.t < t1:
    Prob0203.integrate(Prob0203.t+dt)
    print(" %g " % Prob0203.t )
    print Prob0203.y

Prob0203.set_initial_value(u0)
Prob0203_solution = []
while Prob0203.successful() and Prob0203.t < t1:
    Prob0203_solution.append( [Prob0203.t+dt,] + list( Prob0203.integrate(Prob0203.t+dt) ) )
# take the transpose of a list of lists
Prob0203_solution = map(list, zip(*Prob0203_solution))

plt.figure(1)
plt.plot( Prob0203_solution[3],Prob0203_solution[4])
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Cannonball trajectory with Drag: Variable density')

# part (b)
def deriv_b(t,u): # return derivatives of the array u
    """
    cf. http://bulldog2.redlands.edu/facultyfolder/deweerd/tutorials/Tutorial-ODEs.pdf

    """
    uxdot = (7.853*10**(-5)) *(u[0]**2 + u[1]**2)**(0.5)*(-u[0])
    uydot = -9.8 + (7.853*10**(-5)) *(u[0]**2 + u[1]**2)**(0.5)*(-u[1])
    return array([ uxdot,uydot,u[0],u[1] ])

Prob0203b = ode(deriv_b).set_integrator('dopri5')
Prob0203b.set_initial_value(u0)
Prob0203b.integrate(41.23)

t1b = 41.225
Prob0203b.set_initial_value(u0)
Prob0203b_solution = []
while Prob0203b.successful() and Prob0203b.t < t1b:
    Prob0203b_solution.append( [Prob0203b.t+dt,] + list( Prob0203b.integrate(Prob0203b.t+dt) ) )
Prob0203b_solution = map(list, zip(*Prob0203b_solution))

plt.figure(2)
plt.plot( Prob0203b_solution[3],Prob0203b_solution[4])
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Cannonball trajectory with Drag: Constant density')


# part (c)

300.**2/9.8*sin(2.*50./180.*N(pi) ) # 9044.15283378558

#parabola trajectory data
Prob0203c_x = [i*10 for i in range(905)]
Prob0203c_y = [ tan(50./180.*N(pi))*x - (9.8/2.)*x**2/(300.*cos(50./180.*N(pi)))**2 for x in Prob0203c_x]

plt.figure(3)
plt.plot( Prob0203_solution[3],Prob0203_solution[4], label="Drag: Variable density")
plt.plot( Prob0203b_solution[3],Prob0203b_solution[4], label="Drag: Constant density")
plt.plot( Prob0203c_x,Prob0203c_y, label="No Drag")
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Trajectories of cannonball with Drag of variable density, Drag of constant density, and no drag')
plt.legend()


#########################
####   PS 4
#########################

####################
###     Problem 1
####################

# (b)

k_Boltz = FundConst[ FundConst["Quantity"].str.contains("Boltzmann") ].loc[49,:]
k_Boltz.Value
k_Boltz.Unit
N_Avog = FundConst[FundConst["Quantity"].str.contains("Avogadro") ]

c_V = float( Decimal(1.5)*(N_Avog.Value)*(k_Boltz.Value))/M_0
c_P = float( Decimal(2.5)*(N_Avog.Value)*(k_Boltz.Value))/M_0
c_V.subs(M_0, 39.948/1000.) # 312.198102337360
c_V.subs(M_0, 131.293/1000.) # 94.9912774647001
c_P.subs(M_0, 39.948/1000.) # 520.330170562267
c_P.subs(M_0, 131.293/1000.) # 158.318795774500

tau = Symbol("tau",real=True)
tau_0 = Symbol("tau_0",real=True)
GAMMA = Symbol("GAMMA",positive=True)
MachNo = Symbol("MachNo",positive=True)

TratiovsMachNo = Eq( tau_0/tau, Rat(1) + (GAMMA - Rat(1))/Rat(2)*MachNo )
