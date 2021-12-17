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
from collections import namedtuple


import numpy as np

import cantera as ct

import matplotlib.pyplot as plt

import pandas as pd

import CCDroplet
import CC_params

from CC_params import CH4_TUP, C7H16_TUP,CH4CH_PARAMS,HEPCH_PARAMS,CH4INLETS_VARY_D_0,HEPINLETS_VARY_D_0

Output_Vary_D = namedtuple('Output_Vary_D',['D_0','df']) # 'df' for (pandas) dataframe

gri30_gas = ct.Solution("gri30.cti","gri30_mix")
prf_gas   = ct.Solution("PRF_171.cti")
ctCH4     = ct.Methane()
ctHep     = ct.Heptane()

def CH4result(key):
    """
    calculations for a single initial droplet size only
    """
    return CCDroplet.main(
        gri30_gas,
        ctCH4,
        "methane",
        CH4_TUP,
        CH4CH_PARAMS,
        CH4INLETS_VARY_D_0[key],
        Hrefchoice='y',
        Tboilo=111.5,
    )

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
                Dsizestr = str(result['cleandat'].D[0]*10**6).split('.')[0]
            except KeyError:
                print "KeyError still! I'll print the type:"
                print type(result)
            filename = CH4_TUP.Fuelstr+str(CCDroplet.NUM_STEPS)+"steps"+"D"+Dsizestr+'v_d'+str(result['cleandat'].v_d[0]).split('.')[0]
            tosave = pd.DataFrame(result['cleandat']._asdict(),index=result['xclean'])

# EY : 20160311 I wanted to implement hdf5 requiring PyTables but pip install
# of PyTables was nontrivial; received warning
# Warning! ***HDF5 library version mismatched error*** 
# Abort trap: 6
# would like to know best practices for serialization of pandas DataFrame
#            try:
#                tosave.to_hdf("./data/"+filename,filename)
#            except:
#                print "We're not in a directory where there's a subdirectory for data.  Saving into working directory."
#                try:
#                    tosave.to_hdf(filename,filename)
#                except:
#                    print "PyTables may not be installed. Try pickle.  Be careful of reading unknown pickles (see pandas doc for warning)"
            try:
                tosave.to_pickle("./data/"+filename+'.pkl')
            except:
                print "We're not in a directory where there's a subdirectory for data.  Saving into working directory."
                tosave.to_pickle(filename+'.pkl')

    
    # (surjective) transforming from list of dicts, that includes raw data, to list of namedtuples of pandas dataframe, which only includes clean data
    results = [pd.DataFrame(result['cleandat']._asdict(),index=result['xclean']) for result in results]
#    D_0s = [int(str(result.D[0]*10**6).split('.')[0]) for result in results]
    D_0s = [str(result.D[0]*10**6).split('.')[0] for result in results]

    results = pd.concat(results,axis=1, keys=D_0s)
    filename = CH4_TUP.Fuelstr+str(CCDroplet.NUM_STEPS)+"steps"+ str(results[D_0s[0]].v_d[0.0]).split('.')[0]
    try:
        results.to_pickle("./data/"+filename+'.pkl')
    except:
        print "We're not in a directory where there's a subdirectory for data.  Saving into working directory."
        results.to_pickle(filename+'.pkl')
    

    return results



def output_from_main(gas,liquid,speciesname,reaction_tuple,chamber_params,inlet_conds,Hrefchoice,T_boilo):
    """
    output_from_main(gas,liquid,speciesname,reaction_tuple,chamber_params,inlet_conds,T_boilo)
    """
    choicesaveresults = raw_input("Save results? (y/n)")

    results = []
    for D_0 in inlet_conds:
        result = CCDroplet.main(gas,liquid,speciesname,reaction_tuple,chamber_params,D_0,Hrefchoice, T_boilo)
        results.append(result)


    if choicesaveresults == 'y':
        for result in results:
            try:
                Dsizestr = str(result['cleandat'].D[0]*10**6).split('.')[0]
            except KeyError:
                print "KeyError still! I'll print the type:"
                print type(result)
            filename = speciesname+"steps"+str(CCDroplet.NUM_STEPS)+"D"+Dsizestr+'v_d'+str(result['cleandat'].v_d[0]).split('.')[0]

            
            tosave = pd.DataFrame(result['cleandat']._asdict(),index=result['xclean'])

            try:
                tosave.to_pickle("./data/"+filename)
            except:
                print "We're not in a directory where there's a subdirectory for data.  Saving into working directory."
                tosave.to_pickle(filename)
    
    # (surjective) transforming from list of dicts, that includes raw data, to list of namedtuples of pandas dataframe, which only includes clean data
    results = [pd.DataFrame(result['cleandat']._asdict(),index=result['xclean']) for result in results]
#    D_0s = [int(str(result.D[0]*10**6).split(',')[0]) for result in results]
    D_0s = [str(result.D[0]*10**6).split(',')[0] for result in results]    

    results = pd.concat(results,axis=1, keys=D_0s)
    
    filename = speciesname+str(CCDroplet.NUM_STEPS)+"steps"+ str(results[D_0s[0]].v_d[0.0]).split('.')[0] +'v_d'
    try:
        results.to_pickle("./data/"+filename)
    except:
        print "We're not in a directory where there's a subdirectory for data.  Saving into working directory."
        results.to_pickle(filename)
    
    return results



def compute_flows_phi_g_v_g(gas,liquid,reaction_tuple,chamber_params,inlet_conds,T_boil,oderesults,Hrefchoice='y'):
    """
    compute_flows_phi_g_v_g = compute_flows_phi_g_v_g(gas,liquid,reaction_tuple,T_boil,P,Hrefchoice)
    
    oderesults - <CCDroplet.T_gDv_d_Tuple> 
    """
    P = inlet_conds.TP[1]
    h_fg,T_boil,rho_l,h_l=CCDroplet._process_liquid_Droplet(gas,liquid,reaction_tuple,T_boil,P,Hrefchoice)

    Fu_ind = gas.species_index(reaction_tuple.Fuelstr)
    Ox_ind = gas.species_index(reaction_tuple.Oxstr)

    gas.X = CCDroplet.phi_to_X(1, reaction_tuple)
    F_over_O = gas.Y[Fu_ind]/gas.Y[Ox_ind]

    flow_inlet,initial_conds = CCDroplet.get_initial_flow_vals(gas,reaction_tuple,chamber_params,inlet_conds,rho_l,F_over_O)

    dotm_Ox0 = flow_inlet.dotm_Ox # kg/s
    
    # Additional Python step: create vectors to store data                             

    # try if it's a named_tuple; then try if it's a panda DataFrame
    try:
        NUM_STEPS = len(oderesults[0])
    except:
        NUM_STEPS = oderesults.count()[0]

    dotm_l = np.zeros(NUM_STEPS)
    dotm_g              = np.zeros(NUM_STEPS)
    phi_g               = np.zeros(NUM_STEPS)
    v_g                 = np.zeros(NUM_STEPS)
    
    dotm_l[0] = flow_inlet.dotm_l
    dotm_g[0]           = flow_inlet.dotm_g
    phi_g[0]            = inlet_conds.phi_g
    v_g[0]              = flow_inlet.v_g

    # first try if oderesults is a namedtuple; then try if it's a panda DataFrame
    try:
        D_0 = oderesults.D[0,0]
    except:
        D_0 = oderesults.D.loc[0]

    for k in range(1,NUM_STEPS):
        # first try if oderesults is a namedtuple; then try if it's a panda DataFrame
        try:
            dotm_lx = flow_inlet.dotm_l* (oderesults.D[k]/D_0)**3
        except:
            dotm_lx = flow_inlet.dotm_l*(oderesults.D.iloc[k]/D_0)**3

        dotm_gx = flow_inlet.dotm_g + flow_inlet.dotm_l - dotm_lx
        phi_gx  = ((dotm_gx/dotm_Ox0)-1)/F_over_O
        
        gas.X = CCDroplet.phi_to_X(phi_gx,reaction_tuple)
        MW_gphi_g = gas.mean_molecular_weight

        try:
            v_gx    = dotm_gx*ct.gas_constant*oderesults.T_g[k]/(MW_gphi_g*P*chamber_params.A_cc)
        except:
            v_gx    = dotm_gx*ct.gas_constant*oderesults.T_g.iloc[k]/(MW_gphi_g*P*chamber_params.A_cc)

        dotm_l[k]                   = dotm_lx
        dotm_g[k]                   = dotm_gx
        phi_g[k]                    = phi_gx
        v_g[k]                      = v_gx

    flows = CCDroplet.Flow_Tuple(dotm_l=dotm_l,dotm_Ox=dotm_Ox0,dotm_g=dotm_g,v_g=v_g)
    
    return flows, phi_g
#    return flow_inlet,dotm_l,dotm_g, phi_g

def compute_flows_phi_g_v_g_vary_D_0(gas,liquid,reaction_tuple,chamber_params,inlet_conds_vary_D_0,T_boil,oderesultsdf,Hrefchoice='y'):
    """
    compute_flows_phi_g_v_g_vary_D_0
    """
    assert len(inlet_conds_vary_D_0) == len(oderesultsdf.columns.levels[0])

    number_D_0s = len(inlet_conds_vary_D_0) # number of D_0s 
     
    choicesaveresults = raw_input("Save results? (y/n)")    

    D_0s = oderesultsdf.columns.levels[0]
    results = []
    for i in range(number_D_0s):
        inlet_cond = inlet_conds_vary_D_0[i]
        D_0label   = D_0s[i]

        flow, phi_g = compute_flows_phi_g_v_g(gas,liquid,reaction_tuple,chamber_params,inlet_cond,T_boil,oderesultsdf[D_0label],Hrefchoice)

        number_flowpts = len(flow[0]) # number of flow pts.
        flowdf = pd.DataFrame(flow._asdict(),index=oderesultsdf.index.values[:number_flowpts])

        number_phi_g = len(phi_g) # number of phi_g pts.
        phi_gdf = pd.DataFrame(phi_g,index=oderesultsdf.index.values[:number_phi_g],columns=['phi_g'])
        flow_phi_gdf = pd.concat([flowdf,phi_gdf],axis=1)
        results.append( flow_phi_gdf )
    
    results = pd.concat(results,axis=1,keys=D_0s)

    if choicesaveresults=='y':
    
        filename = "flows"+reaction_tuple.Fuelstr+str(CCDroplet.NUM_STEPS)+"steps"+ str(oderesultsdf[(D_0s[0],'v_d')].iloc[0]).split('.')[0] + 'v_d'
        try:
            results.to_pickle("./data/"+filename+'.pkl')
        except:
            print "We're not in a directory where there's a subdirectory for data.  Saving into working directory."
            results.to_pickle(filename)
    return results

        

def load_droplet_dat(filename):
    """
    load_droplet_dat = load_droplet_dat(filename)
    """
    dat = np.load(filename)

    return {'x'       : dat[0], 
               'xclean'  : dat[4],
               'rawdat'  : CCDroplet.T_gDv_d_Tuple(T_g=dat[1],D=dat[2],v_d=dat[3]),
               'cleandat': CCDroplet.T_gDv_d_Tuple(T_g=dat[5],D=dat[6],v_d=dat[7])
               }

             



"""
EXAMPLES of USAGE
CH4RES = CH4results()
# or if you're importing from outside, i.e. import ccdroplet
# from ccdroplet import CC_out
CH4RES = CC_out.CH4results()

# use pandas to deal with our data!
import pandas as pd





# say the data from solving the ODE system for a droplet of initial size 100 microns for methane, for 500 pts. was solved for in datCCDroplet500CH4D100.npy
dat500100 = load_droplet_dat("datCCDroplet500CH4D100.npy") 
flowphi500100 = compute_flows_phi_g_v_g(gri30_gas,ctCH4,CH4_TUP,CH4CH_PARAMS,CH4INLETS_VARY_D_0[3],111.,dat500100['cleandat'])

# If you already have the results from solving the ode's saved as a panda DataFrame
import pandas as pd
import ccdroplet
from ccdroplet import CC_out
CH4RES = pd.read_pickle("./data/CH41000steps10")
test0 =CC_out.compute_flows_phi_g_v_g(CC_out.gri30_gas,CC_out.ctCH4,CC_out.CH4_TUP,CC_out.CH4CH_PARAMS,CC_out.CH4INLETS_VARY_D_0[0],111.5,CH4RES['30'])

CH4flows = CC_out.compute_flows_phi_g_v_g_vary_D_0(CC_out.gri30_gas,CC_out.ctCH4,CC_out.CH4_TUP,CC_out.CH4CH_PARAMS,CC_out.CH4INLETS_VARY_D_0,111.5,CH4RES)

D_0s = CH4flows.columns.levels[0]

# Examples of plotting now
import matplotlib.pyplot as plt
from matplotlib import cm

CH4flows[[(D_0,'phi_g') for D_0 in D_0s]].plot()

CH4RES[[(D_0,'T_g') for D_0 in D_0s]].plot.area(stacked=False)

T_g_D_v_dstr = [str(col) for col in CH4RES.columns.levels[1]]
T_g_D_v_dstr[0] = "$"+T_g_D_v_dstr[0]+"$ (K)"
T_g_D_v_dstr[1] = "$"+T_g_D_v_dstr[1]+"$ (m)"
T_g_D_v_dstr[2] = "$"+T_g_D_v_dstr[2]+"$ (m/s)"
plttitles = [y+" vs. $x$ (m) for methane ("+CC_out.CH4_TUP.Fuelstr+") with "+CC_out.CH4_TUP.Oxstr+" as oxidizer for various $D_0 (\mu m)$" for y in T_g_D_v_dstr]

# cf. http://stackoverflow.com/questions/21487329/add-x-and-y-labels-to-a-pandas-plot
CH4T_g_ax = CH4RES[[(D_0,'T_g') for D_0 in D_0s]].plot(title= plttitles[0],grid=True,colormap=cm.Blues, style=['*',':','-.','--','-'],linewidth=7 ) # cm.hot # cm.cool # cm.winter # cm.Spectral
CH4T_g_ax.set_xlabel("Combustion chamber length $x$ (m)")
CH4T_g_ax.set_ylabel("gas temp. $T_g$ (K)")
CH4T_g_ax.legend(loc=4)

CH4D_ax = CH4RES[[(D_0,'D') for D_0 in D_0s]].plot(title= plttitles[1],grid=True ,colormap=cm.Blues,style=['*',':','-.','--','-'] ,linewidth=7)
CH4D_ax.set_xlabel("Combustion chamber length $x$ (m)")
CH4D_ax.set_ylabel("droplet diameter $D$ (m)")
CH4D_ax.legend()

# cf.http://stackoverflow.com/questions/26721353/normalize-each-column-of-the-pandas-dataframe
CH4Dnorm = CH4RES[[(D_0,'D') for D_0 in D_0s]]
for col in CH4Dnorm:
    CH4Dnorm[col] /= CH4Dnorm[col].iloc[0]
Dnormtitle = "$D/D_0$ vs. $x$ (m) for methane ("+CC_out.CH4_TUP.Fuelstr+") with "+CC_out.CH4_TUP.Oxstr+" as oxidizer for various $D_0 (\mu m)$"
CH4Dnorm_ax = CH4Dnorm.plot(title=Dnormtitle,grid=True,colormap=cm.Blues,style=['*',':','-.','--','-'] ,linewidth=7)
CH4Dnorm_ax.set_xlabel("Combustion chamber length $x$ (m)")
CH4Dnorm_ax.set_ylabel("$D/D_0$")
CH4Dnorm_ax.legend()


CH4v_d_ax = CH4RES[[(D_0,'v_d') for D_0 in D_0s]].plot(title= plttitles[2],grid=True,colormap=cm.Blues,style=['*',':','-.','--','-'] ,linewidth=7)
CH4v_d_ax.set_xlabel("Combustion chamber length $x$ (m)")
CH4v_d_ax.set_ylabel("droplet speed $v_d$ (m/s)")
CH4v_d_ax.legend(loc=4)

phi_gtitle = "$\phi_g$ vs. $x$ (m) for methane ("+CC_out.CH4_TUP.Fuelstr+") with "+CC_out.CH4_TUP.Oxstr+" as oxidizer for various $D_0 (\mu m)$"
CH4phi_g_ax = CH4flows[[(D_0,'phi_g') for D_0 in D_0s]].plot(title=phi_gtitle,grid=True,colormap=cm.Blues,style=['*',':','-.','--','-'] ,linewidth=7)
CH4phi_g_ax.set_xlabel("Combustion chamber length $x$ (m)")
CH4phi_g_ax.set_ylabel("$\phi_g$ equivalence ratio of fuel to oxidizer for gases only")
CH4phi_g_ax.legend()

v_gtitle = "$v_g (m/s)$ vs. $x$ (m) for methane ("+CC_out.CH4_TUP.Fuelstr+") with "+CC_out.CH4_TUP.Oxstr+" as oxidizer for various $D_0 (\mu m)$"
CH4v_g_ax = CH4flows[[(D_0,'v_g') for D_0 in D_0s]].plot(title=v_gtitle,grid=True,colormap=cm.Blues,style=['*',':','-.','--','-'] ,linewidth=7)
CH4v_g_ax.set_xlabel("Combustion chamber length $x$ (m)")
CH4v_g_ax.set_ylabel("$v_g$ speed of gas (m/s)")
CH4v_g_ax.legend()

CH4flowstr = [(D_0,flow) for D_0 in D_0s for flow in CH4flows.columns.levels[1][:3]]

CH4flowstr = [(D_0,flow) for D_0 in D_0s for flow in CH4flows.columns.levels[1][1:3]]

dotm_gtitle = "$\dot{m}_g (kg/s)$ with $\dot{m}_{Ox} (kg/s)$ vs. $x$ (m) for methane ("+CC_out.CH4_TUP.Fuelstr+") with "+CC_out.CH4_TUP.Oxstr+" as oxidizer for various $D_0 (\mu m)$"
CH4dotm_gOx_ax = CH4flows[CH4flowstr].plot(title=dotm_gtitle,grid=True,colormap=cm.Blues,style=['*','*',':',':','-.','-.','--','--','-','-'] ,linewidth=5)
CH4dotm_gOx_ax.set_xlabel("Combustion chamber length $x$ (m)")
CH4dotm_gOx_ax.set_ylabel("$\dot{m}_g$ and $\dot{m}_{Ox} (kg/s)$")
CH4dotm_gOx_ax.legend(loc=4)



dotm_ltitle = "$\dot{m}_l (kg/s)$ with $\dot{m}_{Ox} (kg/s)$ vs. $x$ (m) for methane ("+CC_out.CH4_TUP.Fuelstr+") with "+CC_out.CH4_TUP.Oxstr+" as oxidizer for various $D_0 (\mu m)$"

CH4flowstr = [(D_0,flow) for D_0 in D_0s for flow in CH4flows.columns.levels[1][:2]]
CH4dotm_lOx_ax = CH4flows[CH4flowstr].plot(title=dotm_ltitle,grid=True,colormap=cm.Blues,style=['*','*',':',':','-.','-.','--','--','-','-'] ,linewidth=5)
CH4dotm_lOx_ax.set_xlabel("Combustion chamber length $x$ (m)")
CH4dotm_lOx_ax.set_ylabel("$\dot{m}_l$ and $\dot{m}_{Ox} (kg/s)$")
CH4dotm_lOx_ax.legend(loc=2)

# overlay plots
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
CH4RES[[(D_0,'T_g') for D_0 in D_0s]].plot(title= plttitles[0],grid=True,colormap=cm.Blues, style=['*',':','-.','--','-'],linewidth=10, ax=ax)
CH4flows[[(D_0,'phi_g') for D_0 in D_0s]].plot(grid=True,colormap=cm.Purples, style=['*',':','-.','--','-'],linewidth=4, ax=ax2)
ax.set_ylabel("gas temp. $T_g$ (K)")
ax2.set_ylabel("$\phi_g$ equivalence ratio of fuel to oxidizer for gases only")   

CH4v_d_v_g = pd.concat([CH4RES[[(D_0,'v_d') for D_0 in D_0s]],CH4flows[[(D_0,'v_g') for D_0 in D_0s]]],axis=1)

CH4v_d_v_g.plot(title="CH4 $v_g,v_d$ (m/s) vs. x (m)",grid=True,colormap=cm.Blues,style=['*',':','-.','--','-','*',':','-.','--','-'],linewidth=5)


# Now let's do it for n-Heptane
import ccdroplet
from ccdroplet import CC_out
HEPRES = CC_out.output_from_main(CC_out.prf_gas,CC_out.ctHep,"heptane",CC_out.C7H16_TUP,CC_out.HEPCH_PARAMS,CC_out.HEPINLETS_VARY_D_0,'n',371.5)

HEPflows = CC_out.compute_flows_phi_g_v_g_vary_D_0(CC_out.prf_gas,CC_out.ctHep,CC_out.C7H16_TUP,CC_out.HEPCH_PARAMS,CC_out.HEPINLETS_VARY_D_0,371.5,HEPRES,'n')


D_0s = HEPflows.columns.levels[0]

T_g_D_v_dstr = [str(col) for col in HEPRES.columns.levels[1]]
T_g_D_v_dstr[0] = "$"+T_g_D_v_dstr[0]+"$ (K)"
T_g_D_v_dstr[1] = "$"+T_g_D_v_dstr[1]+"$ (m)"
T_g_D_v_dstr[2] = "$"+T_g_D_v_dstr[2]+"$ (m/s)"
plttitles = [y+" vs. $x$ (m) for n-heptane ("+CC_out.C7H16_TUP.Fuelstr+") with "+CC_out.C7H16_TUP.Oxstr+" as oxidizer for various $D_0 (\mu m)$" for y in T_g_D_v_dstr]

HEPT_g_ax = HEPRES[[(D_0,'T_g') for D_0 in D_0s]].plot(title= plttitles[0],grid=True,colormap=cm.Blues, style=['*',':','-.','--','-'],linewidth=7 ) # cm.hot # cm.cool # cm.winter # cm.Spectral
HEPT_g_ax.set_xlabel("Combustion chamber length $x$ (m)")
HEPT_g_ax.set_ylabel("gas temp. $T_g$ (K)")
HEPT_g_ax.legend()

HEPD_ax = HEPRES[[(D_0,'D') for D_0 in D_0s]].plot(title= plttitles[1],grid=True ,colormap=cm.Blues,style=['*',':','-.','--','-'] ,linewidth=7)
HEPD_ax.set_xlabel("Combustion chamber length $x$ (m)")
HEPD_ax.set_ylabel("droplet diameter $D$ (m)")
HEPD_ax.legend()


HEPv_d_ax = HEPRES[[(D_0,'v_d') for D_0 in D_0s]].plot(title= plttitles[2],grid=True,colormap=cm.Blues,style=['*',':','-.','--','-'] ,linewidth=7)
HEPv_d_ax.set_xlabel("Combustion chamber length $x$ (m)")
HEPv_d_ax.set_ylabel("droplet speed $v_d$ (m/s)")
HEPv_d_ax.legend(loc=4)

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
HEPRES[[(D_0,'T_g') for D_0 in D_0s]].plot(title= plttitles[0],grid=True,colormap=cm.Blues, style=['*',':','-.','--','-'],linewidth=10, ax=ax)

HEPflows[[(D_0,'phi_g') for D_0 in D_0s]].plot(grid=True,colormap=cm.Purples, style=['*',':','-.','--','-'],linewidth=4, ax=ax2)
ax.set_ylabel("gas temp. $T_g$ (K)")
ax2.set_ylabel("$\phi_g$ equivalence ratio of fuel to oxidizer for gases only")   

HEPv_d_v_g = pd.concat([HEPRES[[(D_0,'v_d') for D_0 in D_0s]],HEPflows[[(D_0,'v_g') for D_0 in D_0s]]],axis=1)

HEPv_d_v_g.plot(title="n-heptane $v_g,v_d$ (m/s) vs. x (m)",grid=True,colormap=cm.Blues,style=['*',':','-.','--','-','*',':','-.','--','-'],linewidth=5)




"""

