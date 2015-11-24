## toPd.py
## Using Pandas (Panda dataframe) for NIST Physical Constants
## NIST National Institute of Standards and Technology 
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

import re
import decimal
from decimal import Decimal, InvalidOperation

import pandas as pd
import numpy as np

from scrape_BS import scraping_allascii, make_conv_lst, make_NASAPLNFACTS, make_Braeunig_ATMOS

###################################
#### Preprocessing for Pandas
###   From a list
###################################
lines,title,src,header,rawtbl,tbl=scraping_allascii()
FundConst = pd.DataFrame(tbl,columns=header)

#########################
#### Examples of usage
#########################

# FundConst.columns ## Index([u'Quantity', u'Value', u'Uncertainty', u'Unit'], dtype='object')

# ks = FundConst[ FundConst["Quantity"].str.contains("Boltzmann") ]
## look at what you want and see the index; it happens to be 49 in this example.
# ks.loc[49,:]
# ks.loc[49,:].Quantity
## 'Boltzmann constant'
# ks.loc[49,:].Value
## Decimal('1.38064852E-23')
# ks.loc[49,:].Unit
## 'J K^-1'
# ks.loc[49,:].Uncertainty
## Decimal('7.9E-30')

###############################################################################################
# NIST Official conversions
# from http://physics.nist.gov/cuu/Reference/unitconversions.html
# Appendix Bof NIST Special Publication 811
# are the official values of conversion factors
###############################################################################################
##################################################
### NIST Guide to SI
## B.8 Factors for Units Listed Alphabetically
##################################################

def make_pd_conv_1st(filename="DF_conv"):
    """
    make_pd_conv_lst
    
    DESCRIPTION
    run make_pd_conv_1st first to put the panda DataFrame saved locally
    """
    headers,convdata,convdata2 = make_conv_lst()
    DF_conv = pd.DataFrame(convdata2,columns=headers)
    DF_conv.to_pickle('./rawdata/'+filename) 
    return DF_conv

conv = pd.read_pickle('./rawdata/DF_conv')

##################################################
##### NASA Planetary Data
##################################################

def make_pd_NASAPLNFACTS(filename="DF_plnfacts"):
    """
    make_pd_NASAPLNFACTS

    DESCRIPTION
    Make a pandas (panda? singular or plural?) DataFrame from html table of 
    NASA Planetary Fact Sheet
    """
    tbl = make_NASAPLNFACTS()
    # preprocess, data wrangle planets header row
    tbl[0][0] = 'Planet'

    # clean the unicode for each planet cf. http://stackoverflow.com/questions/15321138/removing-unicode-u2026-like-characters-in-a-string-in-python2-7
    for i in range(1,len(tbl[0])):
        tbl[0][i] = tbl[0][i].encode('ascii','ignore')
    
    # take the "transpose" of the planetary data
    # cf. http://stackoverflow.com/questions/6473679/python-list-of-lists-transpose-without-zipm-thing
    tbl = map(list,zip(*tbl))
    
    # strip the footnotes * character
    tbl=[[entry.replace("*","") for entry in row] for row in tbl]
    
    # detach header from data
    header = tbl[0]
    tbldata = tbl[1:]
    
    # preprocess data from strings to specific types (decimals and integers)
    tbldata2 = []
    for row in tbldata:
        rowdecis = []
        for j in range(1,18):
            try:
                rowdecis.append( Decimal(row[j].replace(",","")))
            except decimal.InvalidOperation:
                rowdecis.append(None)        
        tbldata2.append( [row[0],]+rowdecis + [int(row[18]),]+row[19:] )
    
    DF_plnfacts = pd.DataFrame(tbldata2,columns=header) # pandas DataFrame of Planetary Fact Sheet
#    DF_conv.to_pickle('./rawdata/'+filename)    
    DF_plnfacts.to_pickle('./rawdata/'+filename)

    return DF_plnfacts

try:
    plnfacts = pd.read_pickle('./rawdata/DF_plnfacts')
except IOError:
    pass

######################################################################
##### Rocket & Space Technology
##### Robert A. Braeunig
####   http://www.braeunig.us/space/index.htm
######################################################################

def make_pd_Braeunig_ATMOS(filename="DF_Braeunig_ATMOS"):
    standatmdict = make_Braeunig_ATMOS()
    DF_standatm  = pd.DataFrame(standatmdict["data"],columns=standatmdict["header"])
    DF_standatm.to_pickle('./rawdata/'+filename)
    return DF_standatm

try:
    Braeunig_standatm = pd.read_pickle('./rawdata/DF_Braeunig_ATMOS')
except IOError:
    pass
    

#########################
#### Examples of usage
#########################

# If you are running this for the FIRST TIME, be sure to do the following commands to add data locally to the subdirectory 'rawdata':
## make_pd_conv_1st()
## make_pd_NASAPLNFACTS()


###############
##### main
###############

if __name__ == "__main__":
    print "FundConst - NIST Fundamental Constants as a Panda DataFrame "
    print "conv      - NIST SI Conversions as a Panda DataFrame"


