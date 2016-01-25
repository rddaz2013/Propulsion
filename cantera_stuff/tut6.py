## tut6.py
## tut6.m implemented in Python for cantera 
## cf. http://www.cantera.org/docs/sphinx/html/matlab/tutorials/tut6.html
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
# Tutorial 6:   Transport properties
#
#    Topics:
#       - mixture-averaged and multicomponent models
#       - viscosity
#       - thermal conductivity
#       - binary diffusion coefficients
#       - mixture-averaged diffusion coefficients
#       - multicomponent diffusion coefficients
#       - thermal diffusion coefficients
#

#################################################################
import cantera as ct
import numpy as np

# Methods are provided to compute transport properties. By
# default, calculation of transport properties is not enabled. If
# transport properties are required, the transport model must be
# specified when the gas mixture object is constructed.

# Currently, two models are implemented. Both are based on kinetic
# theory expressions, and follow the approach described in Dixon-Lewis
# (1968) and Kee, Coltrin, and Glarborg (2003). The first is a full
# multicomponent formulation, and the second is a simplification that
# uses expressions derived for mixtures with a small number of species
# (1 to 3), using approximate mixture rules to average over
# composition.

# To use the multicomponent model with GRI-Mech 3.0, call function
# GRI30 as follows:

g1 = ct.Solution('gri30.xml')
g1.transport_model='Multi'

# To use the mixture-averaged model:

g2 = ct.Solution('gri30.xml')
g2.transport_model='Mix'


# Both models use a mixture-averaged formulation for the viscosity.
visc = [g1.viscosity, g2.viscosity]


# The thermal conductivity differs, however.
lambdas = [ g1.thermal_conductivity, g2.thermal_conductivity ] # lambda is in namespace of Python of names you CANNOT use; one will obtain SyntaxError: invalid syntax if you try to assign anything to lambda


# Binary diffusion coefficients
bdiff1 = g1.binary_diff_coeffs
bdiff2 = g2.binary_diff_coeffs


# Mixture-averaged diffusion coefficients. For convenience, the
# multicomponent model implements mixture-averaged diffusion
# coefficients too.
dmix1 = g1.mix_diff_coeffs
dmix2 = g2.mix_diff_coeffs

# EY : 20160125 cf. https://code.google.com/p/cantera/source/browse/cantera/trunk/interfaces/matlab/toolbox/@Transport/mixDiffCoeffs.m?r=1126


# Multicomponent diffusion coefficients. These are only implemented
# if the multicomponent model is used. 
dmulti = g1.multi_diff_coeffs


# Thermal diffusion coefficients. These are only implemented with the
# multicomponent model.  These will be very close to zero, since
# the composition is pure H2.
dt= g1.thermal_diff_coeffs


# Now change the composition and re-evaluate
g1.X = np.ones(g1.n_species)
dt = g1.thermal_diff_coeffs 

# Note that there are no singularities for pure gases. This is
# because a very small positive value is added to all mole
# fractions for the purpose of computing transport properties. 

#################################################################


#################################################################
#   end of tutorial 6
#################################################################
