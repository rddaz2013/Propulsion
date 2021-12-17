## equil.py
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

def equil(g=None):
    gas = ct.Solution("gri30.cti") if g is None else g
    nsp = gas.n_species

    # find methane, nitrogen, and oxygen indices
    ich4 = gas.species_index('CH4')
    io2  = gas.species_index('O2')
    in2  = gas.species_index('N2')

    phi = [0.2 + 0.05*(i+1) for i in range(50)]

    tad = []
    xeq = []
    for phi_i in phi:
        x = np.zeros(nsp)
        x[ich4] = phi_i
        x[io2] = 2.0
        x[in2] = 7.52
        gas.TP = 300.0, ct.one_atm
        gas.X  = x
        gas.equilibrate('HP')
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
    
    
