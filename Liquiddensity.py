## Liquiddensity.py
## Obtaining Liquid density data from NIST Chemistry WebBook
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160227
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

WEBBK_URL = "http://webbook.nist.gov"
WEBBK_NAME_URL = "http://webbook.nist.gov/cgi/cbook.cgi?Name=methane&Units=SI"

FLUID_URLQ = "http://webbook.nist.gov/cgi/fluid.cgi?TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoBar&RefState=DEF&Action=Page&ID=C74828"

FLUID_URLQ_HTML = "http://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C74828&Type=IsoBar&Digits=5&P=3.4&THigh=174.1&TLow=174.1&TInc=&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm"


def _get_Fluid_data(T,P,name="methane"):
    """
    _get_Fluid_data = _get_Fluid_data(name="methane")
    Returns Fluid change data from the NIST Chemistry Webbook as a BeautifulSoup object

    INPUTS/Parameters
    -----------------
    T              - <float>
    Temperature T (K) 
    
    P              - <float>
    pressure P (Pa)

    name="methane" - <string>
    Name of the (chemical) species that someone wants to obtain Phase data on
    
    """
    P_MPa = P/1000000. # Pascal -> MPa

    Name_urlparse = urlparse.urlparse( WEBBK_NAME_URL )
    Name_url_parts = list(Name_urlparse)
    namequery = urlparse.parse_qs( Name_urlparse.query )

    namequery['Name'] = name
    Name_url_parts[4] = urlencode(namequery)
    target_url = urlparse.urlunparse(Name_url_parts)

    url = requests.get( target_url )
    soup = BeautifulSoup( url.content, "html.parser")

    FluidURL = WEBBK_URL + soup.find('a',text="Fluid Properties")['href']

    Fluid_urlparse = urlparse.urlparse(FluidURL)
    Fluid_ID = urlparse.parse_qs( Fluid_urlparse.query )['ID'][0] # something like u'C74828'
    
    # Isobaric Properties for Fluid, directly obtained
    Isobar_urlparse = urlparse.urlparse( FLUID_URLQ_HTML )
    Isobar_url_parts = list(Isobar_urlparse)
    Isobarquery = urlparse.parse_qs( Isobar_urlparse.query )

    for key in Isobarquery:
        Isobarquery[key] = Isobarquery[key][0]
    Isobarquery['ID'] = Fluid_ID
    Isobarquery['THigh'] = str(T)
    Isobarquery['TLow'] = str(T)
    Isobarquery['P']    = str(P_MPa)
    Isobar_url_parts[4] = urlencode( Isobarquery )
# EY : 20160227 Something's wrong with my urlencode where I obtain urls where the special characters are turned into escaped html codes such as this:
# Digits=%5B%275%27%5D&RefState=%5B%27DEF%27%5D&DUnit=%5B%27kg%2Fm3%27%5D&THigh=%5B%27174.1%27%5D&HUnit=%5B%27kJ%2Fkg%27%5D&WUnit=%5B%27m%2Fs%27%5D&VisUnit=%5B%27Pa%2As%27%5D&TLow=%5B%27174.1%27%5D&PUnit=%5B%27MPa%27%5D&P=%5B%273.4%27%5D&TUnit=%5B%27K%27%5D&STUnit=%5B%27N%2Fm%27%5D&Action=%5B%27Load%27%5D&Type=%5B%27IsoBar%27%5D&ID=%5B%27C74828%27%5D

    Isobar_target_url = urlparse.urlunparse(Isobar_url_parts)
    
    Isobarurl = requests.get( Isobar_target_url)
    Isobarsoup = BeautifulSoup( Isobarurl.content, "html.parser")

    Isobartbl = Isobarsoup.find("h2",text="Fluid Data").find_next("table")

    Isobarhdrs = [hdr.text for hdr in Isobartbl.find('tr').find_all('th')]
    Isobarvalues = []
    for row in Isobartbl.find_all('tr')[1:]:
        Isobarvalues.append( [td.text for td in row.find_all('td')] )
        
    Isobarvaluesdict = [dict(zip(Isobarhdrs,row)) for row in Isobarvalues]
    
    # Sanity check
    if len(Isobarvaluesdict)==1: # which should be the case, since only 1 value of T was specified 
        Isobarvaluesdict = Isobarvaluesdict[0]
    else: # something went wrong with url response back
        print "Something went wrong with the url response back!"
        
    return Isobarvaluesdict

def fluid_density(T,P,name="methane"):
    """
    fluid_density = fluid_density(T,P,name="methane")
    """
    Isobardict = _get_Fluid_data(T,P,name)
    return float( Isobardict['Density (kg/m3)'] )


if __name__ == "__main__":
    print "The functions to use are \n \t fluid_density = fluid_density(T,P,name)\n"
    

