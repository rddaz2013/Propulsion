## LiquidVaporEq.py
## Obtaining Liquid Vapor Equilibrium data from NIST Chemistry WebBook
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160205
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
import requests 
from bs4 import BeautifulSoup

import urlparse
from urllib import urlencode

from collections import namedtuple, OrderedDict

import sympy
from sympy import Eq, exp, log, Symbol
from sympy.abc import A,B,C,p,T
from sympy.plotting import plot

WEBBK_URL = "http://webbook.nist.gov"
WEBBK_NAME_URL = "http://webbook.nist.gov/cgi/cbook.cgi?Name=methane&Units=SI"

# PHYSICAL CONSTANTS
GASCONSTANT = 8314.4621 # J/(kmol * K )
ONEATM      = 101325.0 # pascal


AntoineEqn = Eq(p, 10**A*exp(-B*log(10)/(T+C) ) )

PhaseData = namedtuple('PhaseData',['MW','T_boil','AntEqhdrs','AntEqParams','DeltaH_vap'])

def _get_Phase_data(name="methane"):
    """
    _get_Phase_data = _get_Phase_data(name="methane")
    Returns Phase change data from the NIST Chemistry Webbook as a BeautifulSoup object

    INPUTS/Parameters
    -----------------
    name="methane" - <string>
    Name of the (chemical) species that someone wants to obtain Phase data on

    OUTPUTS/Results
    ---------------
    Phasesoup - <BeautifulSoup object>
    BeautifulSoup object from the url contents of the Phase change data webpage of the 
    specified substance

    EXAMPLES of USAGE:
    ------------------
    propane_dat = _get_Phase_data("propane")
    """
    Name_urlparse = urlparse.urlparse( WEBBK_NAME_URL )
    Name_url_parts = list(Name_urlparse)
    namequery = urlparse.parse_qs( Name_urlparse.query )

    namequery['Name'] = name
    Name_url_parts[4] = urlencode(namequery)
    target_url = urlparse.urlunparse(Name_url_parts)

    url = requests.get( target_url )
    soup = BeautifulSoup( url.content, "html.parser")

    PhaseURL = WEBBK_URL + soup.find('a',string="Phase change data")['href']

    Phaseurl = requests.get( PhaseURL )
    Phasesoup = BeautifulSoup( Phaseurl.content, "html.parser")

    return Phasesoup
    
def cleaned_Phase_data(name="methane"):
    """
    cleaned_Phase_data = cleaned_Phase_data(name="methane")
    Returns Phase change data from the NIST Chemistry Webbook that is already "cleaned" or 
    formatted for use in Python, in the form of a named tuple, PhaseData
    
    INPUTS/Parameters
    -----------------
    name="methane" - <string>
    Name of the (chemical) species that someone wants to obtain Phase data on

    OUTPUTS/Results
    ---------------
    phasedat_result - <namedtuple object>
    Note that the instantiation or "declaration" of this namedtuple is done globally in this file. 
    phasedat_result, and PhaseData namedtuple in general, has the following entry categories, 
    or "keys": "MW", "T_boil", "AntEqhdrs", "AntEqParams", "DeltaH_vap", 
    where "AntEq" is Antoine Equation, "hdrs" is headers, "Params" is Parameters, 
    "vap" is vaporization    

    EXAMPLES of USAGE:
    ------------------
    propane_dat = cleaned_Phase_data("propane")
    """

    Phasesoup = _get_Phase_data(name)
    # Phase change data headers
    Phasechangedatahdrs = [hdr.string for hdr in Phasesoup.find_all('table')[2].find_all('tr')[0].find_all('th')]
    
    # Some substances will have a boiling point; some won't
                           
    # Temperatures table for Phase change                              
    tempstbl = Phasesoup.find_all('table')[2].find_all('tr')[1].find_all('td') 
    if tempstbl[0].text == u'Tboil':
        Tboiltext = [td.text for td in tempstbl]
        Tboildict = OrderedDict( zip( Phasechangedatahdrs , Tboiltext) )
        Tboildict['Valuestr'] = Tboildict['Value']
        # just the value of T_boil as a float with no uncertainties, no plus-minus values
        Tboildict['Value'] = float(Tboildict['Valuestr'].split()[0])
    else:
        Tboildict = None

    # obtain table for Antoine Equation Parameters
    AntEqsoup = Phasesoup.find('table',summary="Antoine Equation Parameters")
    AntEqhdrs = [hdr.string for hdr in AntEqsoup.find('tr').find_all('th')]

    AntEqtbl = []
    for row in AntEqsoup.find_all('tr')[1:]:
        AntEqtbl.append( [td.text for td in row.find_all('td')] )

    # split up the temperature range
    AntEqtbl_split = []
    for row in AntEqtbl:
        temp = row[0].replace(' ','').split('-')
        AntEqtbl_split.append( temp + row[1:] )

    # split up the temperature slot for the headers
    AntEqhdrs_split = [ AntEqhdrs[0] +u' lo', AntEqhdrs[0]+u' hi']+AntEqhdrs[1:]

    # scrape Molecular Weight 
    MW = Phasesoup.find('a',string="Molecular weight").find_previous('li').text
    MW = float(MW.split(':')[1].strip())
    
    # scrape Enthalpy of vaporization 
    vaptbls = Phasesoup.find_all('table',summary="Enthalpy of vaporization")
    DeltaHtbls = []

    for tbl in vaptbls:
        cleantbl = []
        hdrs = [hdr.text for hdr in tbl.find('tr').find_all('th')]
        cleantbl.append(hdrs)
        for row in tbl.find_all('tr')[1:]:
            cleantbl.append( [td.text for td in row.find_all('td')] )
        DeltaHtbls.append(cleantbl)
    
    # some substances have only 1 Enthalpy of vaporization table; some have more than 2 and
    # are of different format
    # So, unfortunately, the tables on this page isn't "standard" in that the headers are 
    # in a single row: for the second table, it's a single COLUMN.  
    # Instead of trying to mind guess the rationale for this html page layout, 
    # I will hack away expediently.  

    if len(vaptbls) > 1:
        DeltaHvaptbl = []
        for row in vaptbls[1].find_all('tr'): 
            DeltaHvaptbl.append( (row.find('th').text , row.find('td').text) )
        DeltaHtbls[-1] = DeltaHvaptbl

    
    phasedat_result = PhaseData(MW=MW,T_boil=Tboildict,AntEqhdrs=AntEqhdrs_split,AntEqParams=AntEqtbl_split,DeltaH_vap=DeltaHtbls)    

    return phasedat_result


# This is the Clausius-Clapeyron relation for phase coexistence, see my derivation in thermo.pdf or thermo.tex, rewritten in terms of enthalpy of vaporization
h_fg = Symbol('h_fg',real=True)
p_0  = Symbol('p_0',positive=True)
T_b  = Symbol('T_b',positive=True)
ClausClapEqn = Eq(p, p_0*exp(-h_fg*(1/T - 1/T_b) ) )

def fit_ClausClapEqn(T_boil, DeltaH ):
    """
    fit_ClausClapEqn = fit_ClausClapEqn( MW, DeltaHtbl )

    EXAMPLES of USAGE:
    methanedat = cleaned_Phase_data()
    methaneDeltaH = float(methanedat[4][0][2][0]) * 1000 * 1000 # kJ/mol -> J/mol -> J/kmol
    methaneCCeq = fit_ClausClapEqn( methanedat[0], methaneDeltaH )
    plot( methaneCCeq.rhs, (T,50,150) )
    """
    CCEqn = ClausClapEqn

    CCEqn = CCEqn.subs( T_b, T_boil )
    CCEqn = CCEqn.subs( h_fg,DeltaH/GASCONSTANT)
    CCEqn = CCEqn.subs( p_0, ONEATM )
    return CCEqn

def plot_all_Ant_fits( AntEqtbl_split ):
    """
    plot_all_Ant_fits = plot_all_Ant_fits( AntEqtbl_split )

    EXAMPLES of USAGE:
    propane_dat = cleaned_Phase_data("propane")
    propane_plts = plot_all_Ant_fits( propane_dat.AntEqParams )
    """
    fits = []
    for row in AntEqtbl_split:
        fit = AntoineEqn.subs(dict(zip([A,B,C], [float(no) for no in row[2:5]])))
        fits.append( fit )
    to_plot = []
    for fit, row in zip(fits, AntEqtbl_split):
        range = (T, float(row[0]), float(row[1]))
        to_plot.append( (fit.rhs, range ) )
    plot( *to_plot )
    return to_plot

# Not essential function; purpose of these functions are to demonstrate or output data
def comparing_AntEq_vs_ClausClap():
    """
    comparing_AntEq_vs_ClausClap
    """
    CH4dat = cleaned_Phase_data("methane")
    CH4DeltaH = float(CH4dat.DeltaH_vap[0][2][0])*1000*1000 # kJ/mol -> J/mol -> J/kmol
    CH4T_boil = CH4dat.T_boil['Value']
    CH4CCEqn = fit_ClausClapEqn(CH4T_boil, CH4DeltaH)
    CH4fits = []
    for row in CH4dat.AntEqParams:
        fit = AntoineEqn.subs(dict(zip([A,B,C],[float(no) for no in row[2:5]])))
        CH4fits.append(fit)
    return CH4dat, CH4CCEqn, CH4fits

if __name__ == "__main__":
    print "The functions to use are \n \t cleaned_Phase_data \n \t fit_ClausClapEqn \n \t plot_all_Ant_fits \n"
    
    CH4dat, CH4CCEqn, CH4fits = comparing_AntEq_vs_ClausClap()
    
    CH4Antplts = plot_all_Ant_fits(CH4dat.AntEqParams)

    plot( *CH4Antplts,title="Methane phase p vs. T", ylabel="p (bar)", xlabel="T (K)", legend=True)
    plot( CH4CCEqn.rhs,(T,95,180),title="Clausius-Clapeyron equation for Methane, p vs. T",ylabel="p (Pa)",xlabel="T(K)")

    plot(*(CH4Antplts+[(CH4CCEqn.rhs/100000,(T,95,180)),]),title="Compare Clausius-Clapeyron equation and Antoine Equation for Methane, p vs. T",ylabel="p (bar)",xlabel="T(K)")




    




