## CC_out.py
## Combustion Chamber Droplet output
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
import numpy as np

import cantera as ct

import matplotlib.pyplot as plt

import CCDroplet
import CC_params

from CC_params import CH4_TUP, C7H16_TUP,CH4CH_PARAMS,HEPCH_PARAMS,CH4INLETS_VARY_D_0,HEPINLETS_VARY_D_0

gri30_gas = ct.Solution("gri30.cti","gri30_mix")
prf_gas   = ct.Solution("PRF_171.cti")
ctCH4     = ct.Methane()
ctHep     = ct.Heptane()

def CH4result(key):
    """
    calculations for a single initial droplet size only
    """
    result = CCDroplet.main(gri30_gas,ctCH4,"methane",CH4_TUP,CH4CH_PARAMS,CH4INLETS_VARY_D_0[key],Hrefchoice='y', Tboilo=111.5)
    return result

def CH4results():
    choicesaveresults = raw_input("Save results? (y/n)")
    results = []
    for D_0 in CH4INLETS_VARY_D_0:
        result = CCDroplet.main(gri30_gas,ctCH4,"methane",CH4_TUP,CH4CH_PARAMS,D_0,Hrefchoice='y', T_boilo=111.5)
        results.append(result)

    if choicesaveresults == 'y':
        for result in results:
            # sanity check
            try:
                Dsizestr = str((result['rawdat'].D)[0,0])
            except KeyError:
                print "KeyError still! I'll print the type:"
                print type(result)
            filename = "datCCDroplet"+str(CCDroplet.NUM_STEPS)+CH4_TUP.Fuelstr+"D"+str((result['rawdat'].D)[0,0]*10**6).split('.')[0]
            x = np.concatenate([array for array in result['cleandat'] ])
            np.save(filename,x)
        
    return results

"""
EXAMPLES of USAGE
CH4RES = CH4results()



"""

