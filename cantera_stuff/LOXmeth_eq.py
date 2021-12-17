## LOXmeth_eq.py
## equil.m implemented in Python for cantera 
## cf. http://www.cantera.org/docs/sphinx/html/matlab/examples/equil.html
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160122
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
import numpy as np
import matplotlib.pyplot as plt

def equil(g=None,p=ct.one_atm,solver_choice="gibbs"):
    gas = ct.Solution("gri30.xml") if g is None else g
    nsp = gas.n_species

    # find methane, nitrogen, and oxygen indices
    ich4 = gas.species_index('CH4')
    io2  = gas.species_index('O2')

    phi = np.linspace(0.75,5,50)  # O/F mass ratio (ranging from 0.75-5

    tad = []
    xeq = []
    for phi_i in phi:
        y = np.zeros(nsp)
        y[ich4] = 1.0  # we want mass fractions
        y[io2] = phi_i

        gas.TP = 300.0, p
        gas.Y  = y
        gas.equilibrate('HP',solver=solver_choice)
        tad.append( gas.T)
        xeq.append( gas.X )

    # preprocess xeq as a numpy array for easier indexing
    xeq = np.array(xeq)
    return phi, tad, xeq, gas

# make plots
def make_equil_plts(phi, tad, xeq, gas):
    """
    make_equil_plts
    make_equil_plts makes plots from Python function equil
    """
    plt.figure(1)
    plt.subplot(211)
    phitadplt = plt.plot(phi,tad)
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Temperature (K)')
    plt.title('Adiabatic Flame Temperature')

    plt.subplot(212)
    phixeqsemilogyplt = plt.semilogy(phi,xeq) # matplotlib semilogy takes in the correct dimensions! What I mean is that say phi is a list of length 50.  Suppose xeq is a numpy array of shape (50,53).  This will result in 53 different (line) graphs on the same plot, the correct number of line graphs for (50) (x,y) data points/dots!
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Mole Fraction')
    plt.title('Equilibrium Composition')

    j = 10
    for k in range(gas.n_species):
        plt.text(phi[j],1.5*xeq[j,k],gas.species_name(k))
        j += 2
        if j > 46:
            j = 10
    plt.show()
    return phitadplt, phixeqsemilogyplt

def make_equil_plts_2(phi, tad, xeq, gas):
    """
    make_equil_plts_2
    make_equil_plts_2 makes plots from Python function equil; it's the same as make_equil_plts, but with 2 separate figures
    """
    plt.figure(1)
    phitadplt = plt.plot(phi,tad)
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Temperature (K)')
    plt.title('Adiabatic Flame Temperature')

    plt.figure(2)
    plt.ylim(10.**(-14),1)
    plt.xlim(phi[0],phi[49])
    phixeqsemilogyplt = plt.semilogy(phi,xeq) # matplotlib semilogy takes in the correct dimensions! What I mean is that say phi is a list of length 50.  Suppose xeq is a numpy array of shape (50,53).  This will result in 53 different (line) graphs on the same plot, the correct number of line graphs for (50) (x,y) data points/dots!
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Mole Fraction')
    plt.title('Equilibrium Composition')

    j = 10
    for k in range(gas.n_species):
        plt.text(phi[j],1.5*xeq[j,k],gas.species_name(k))
        j += 2
        if j > 46:
            j = 10
    plt.show()
    return phitadplt, phixeqsemilogyplt


# for the mean molecular weight, recall that you do
# gas.mean_molecular_weight
# for all the molecular weights, it's an attribute of the (instantiated) class and its (type is) a numpy array
# gas.molecular_weights

# recall that the ratio of specific heat is also called the heat capacity ratio or specific heat ratio and it's specific heat at constant pressure C_p to the specific heat at constant volume C_V and you can easy get both and divide by doing the following:
# gas.cp
# gas.cv
# gas.cp/gas.cv

def char_velocity( MW, T_0, gamma=1.6 ):
    """
    char_velocity
    Note: char_velocity uses cantera's gas constant, which is in J/(K kmol) not mol

    ARGUMENTS/INPUTS:
    MW <number> is the molecular weight in amu's or (equivalently) g/mol =kg/kmol
    
    """
    ceestar = (1./gamma)*((gamma+1.)/2)**((gamma+1.)/(gamma-1) )*ct.gas_constant/MW*T_0
    ceestar = np.sqrt( ceestar )
    return ceestar

# for (rocket) Thrust Coefficient is C_F, if p_ambient = p_e (not likely, but in the ratio of (p_e - p_ambient)/p_c, if p_e is close to p_c, then take C_F = 1), and p_e=p_c,
# I_sp = C_f c*/g_0 = c*/g_0


# Let's "generalize" equil function to any choice of oxidizer and fuel
# this also includes saving data (as a function of O/F mass ratio) for 
#  molecular weight MWeq
#  ratio of specific heat gammas
#  characteristic velocity char_v

def equil_general(Ox='O2',Fuel='CH3',g=None,p=ct.one_atm,solver_choice="gibbs"):
    """
    equil_general
    equil_general is a "superset" of equil, but "generalizes" to any choice of fuel and oxidizer (specify oxidizer and fuel as a string for the argument Ox, Fuel; e.g. Ox='O2', Fuel='CH3' (which is the default) or Ox='N2O4',Fuel='N2H4')
    equil_general will also save and output values (beyond equil) for 
     * molecular weight MWeq
     * ratio of specific heat gammas
     * characteristic velocity char_v

    OUTPUT/RESULTS:
    phi,tad,xeq,gas,MWeq,gammas,char_v

    EXAMPLES of USAGE:
    phi,tad,xeq,gas,MWeq,gammas,char_v = equil_general() # for methane
    
    spec = ct.Species.listFromFile('Ae121.cti')
    rxns = ct.Reaction.listFromFile('gri30.cti')
    gas = ct.Solution(thermo='IdealGas',kinetics='GasKinetics',species=spec,reactions=rxns)
    phi,tad,xeq,gas,MWeq,gammas,char_v = equil_general(Ox='N2O4',Fuel='N2H4',g=gas,p=ct.one_atm*80.,solver_choice="auto")

    """
    gas = ct.Solution("gri30.xml") if g is None else g
    nsp = gas.n_species

    iFu = gas.species_index(Fuel)  # index of Fuel
    iOx  = gas.species_index(Ox)  # index of Oxidizer

    phi = np.linspace(0.75,5,50)  # O/F mass ratio (ranging from 0.75-5

    tad  = []
    xeq  = []
    MWeq = []
    gammas = []
    char_v = []
    for phi_i in phi:
        y = np.zeros(nsp)
        y[iFu] = 1.0  # we want mass fractions
        y[iOx] = phi_i

        gas.TP = 300.0, p
        gas.Y  = y
        gas.equilibrate('HP',solver=solver_choice)
        tad.append( gas.T)
        xeq.append( gas.X )
        MWeq.append( gas.mean_molecular_weight )
        gammas.append( gas.cp/gas.cv )
        char_v.append( char_velocity( gas.mean_molecular_weight, gas.T, gas.cp/gas.cv) )

    # preprocess xeq as a numpy array for easier indexing
    xeq = np.array(xeq)
    return phi, tad, xeq, gas, MWeq, gammas, char_v

def check_equil(gas):
    """
    check_equil
    check_equil does a baseline check if "equilibrate" has correctly found the 
     chemical equilibrium state by showing that the net rates of progress of 
     all reversible reactions are 0, by outputting the max. net rate of progress
    
    See tut4.py or tut4.m for the rationale for check_equil, to 
    tell if "equilibrate" has correctly found the chemical equilibrium state
    """
    rf = gas.forward_rates_of_progress
    rr = gas.reverse_rates_of_progress
    results = [(i,rf[i],rr[i],(rf[i]-rr[i])/rf[i]) for i in range(gas.n_reactions) if gas.is_reversible(i)]
    results = np.array(results)
    maxnetrate = results[:,3]
    maxnetrate = [i for i in maxnetrate if not np.isnan(i)]
    maxnetrate = np.array(maxnetrate)
    maxnetrate = maxnetrate.max()
    return results, maxnetrate


###################################
##### Special case of Hydrazine
###################################

# spec = ct.Species.listFromFile('Ae121.cti')
# rxns = ct.Reaction.listFromFile('gri30.cti')
# gas = ct.Solution(thermo='IdealGas',kinetics='GasKinetics',species=spec,reactions=rxns)
# Explanation of the above 3 lines 
# cf. http://www.cantera.org/docs/sphinx/html/cython/importing.html
# We want to build the Solution objects using Species objects and Reaction objects separately.  
# 'IdealGas' is a built into Cantera thermodynamics i.e. thermo model  
# 'GasKinetics' will calculate the kinetics; it's a built in kinetics model
# cf. https://github.com/Cantera/cantera/blob/5926d2db7c0d4919b75ee50828b0adab4e691a51/include/cantera/kinetics.h
# it appears these models have to be specified if one sought to create or instantiate a Solution class with the species, spec, separately; 
