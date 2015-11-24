## scrape_BS.py
## Screen scraping i.e. web scraping with Beautiful Soup (BS) and requests
##  
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

import urllib2
import re
import decimal
from decimal import Decimal

import requests
from bs4 import BeautifulSoup as BS

############################################################
## NIST - National Institute of Standards and Technology
############################################################

class scraped_BS(object):
    """
    scraped_BS
    scraped_BS is a class so we can have, as an object, all that results from scraping 
    with requests and BeautifulSoup (BS)
    
    Example of Usage:
    python -i scrape_BS.py
    >>> mainNISTcuu = scraped_BS(NISTCUU_URL)
    >>> NISTconvBS  = scraped_BS(NISTCONValpha)
    """
    def __init__(self,url):
        self.url  = url
        self.req  = requests.get(url)
        self.soup = BS(self.req.content) 
        self.req.close()


#################################################################
## NIST Reference on Constants, Units, and Uncertainty mainpage
#################################################################

FPCascii    = "http://physics.nist.gov/cuu/Constants/Table/allascii.txt" # Fundamental Physical Constants

#####
# On looking at mainNISTcuu.soup.find_all("a") and your favorite Web Inspector of the 
# NIST cuu webpage, NISTCUU_URL, one realizes that to get to Frequently Used Constants 
# or All values, one does a query.  Let's grab all values in ascii.  
# Do this with FPCasciitbl = scraped_BS(FPCascii) if you want to use requests
# Do this with urllib2 as follows, using urllib
#####

def retrieve_file(url=FPCascii,filename="allascii.txt"):
    with open('./rawdata/'+filename,'wb') as targetfile:
        u = urllib2.urlopen(url) 
        targetfile.write(u.read())
        targetfile.close()
    return u

####################
# At this point, running retrieve_file() should put allascii.txt into the rawdata subdirectory
# You can open up allascii.txt and see you have all the Fundamental Physical Constants (!!!)
# Time to parse this table:
####################

def line_reading(filename='./rawdata/allascii.txt'):
    openedfile = open(filename,'rb')
    lines = openedfile.read().splitlines()
    lines = [line for line in lines if line != '']
    return lines

def scraping_allascii(filename='./rawdata/allascii.txt'):
    lines  = line_reading(filename)
    title  = lines[0].strip()
    src    = lines[1].strip()
    header = lines[2].split()

    # cf. http://stackoverflow.com/questions/12866631/python-split-a-string-with-at-least-2-whitespaces for re
    rawtbl = []
    for line in lines[4:]:
        rawtbl.append( re.split(r'\s{2,}', line) )
    tbl = []
    for rawline in rawtbl:
        line = []
        line.append(rawline[0])
        try:
            line.append(Decimal(rawline[1].replace(" ","")))
        except decimal.InvalidOperation:
#            line.append(Decimal(rawline[1].replace(" ","").replace(".","") ))
            value = rawline[1].replace(" ","")
            value = ''.join( re.split(r'[.]{2,}',value))
            line.append( Decimal(value) )
        try:
            line.append(Decimal(rawline[2].replace(" ","")))
        except decimal.InvalidOperation:
#            line.append(rawline[2].replace(" ","")
            # EY : 20150823 for SQLAlchemy, instead of "(exact)" use None type None for 
            # SQLAlchemy to denote "(exact)" in the Uncertainty column
            line.append(None)
        line.append(rawline[3])
        tbl.append(line)

    return lines, title, src, header, rawtbl, tbl



###########################################################################
# Now you should be able to put this into your favorite database of choice
###########################################################################


###############################################################################################
# NIST Official conversions
# from http://physics.nist.gov/cuu/Reference/unitconversions.html
# Appendix Bof NIST Special Publication 811
# are the official values of conversion factors
###############################################################################################
##################################################
### NIST Guide to SI
## B.8 Factors for Units Listed Alphabetically
NISTCONValpha =  "http://physics.nist.gov/Pubs/SP811/appenB8.html"

def make_conv_lst(url=NISTCONValpha,tablecls='texttable'):
    convBS = scraped_BS(url)
    convBS.convtbls = convBS.soup.find_all("table",{"class": tablecls})
    convdata = []
    convdata2 = []

    headers = convBS.convtbls[0].find_all('tr')[1].find_all('th')
    headers = [ele.text.replace(' ','') for ele in headers]

    for tbl in convBS.convtbls:
        for row in tbl.find_all('tr'):
            if row.find_all('td') != []:
                if row.text != '':
                    rowsplit = row.text.replace("\n",'',1).split('\n')
                    try:
                        rowsplit = [pt.replace(u'\xa0',u' ').strip() for pt in rowsplit]
                    except UnicodeDecodeError as err:
                        print rowsplit
                        Break
                        raise err
                    convdata.append( rowsplit )
                    if len(row.find_all('td')) == (len(headers)+1):
                        convdata2.append( row.find_all('td'))
    convdata3 = []
    for row in convdata2:
        rowout = []
        rowout.append( row[0].text.strip())
        rowout.append( row[1].text.strip())
        value = (row[2].text+row[3].text).strip().replace(u'\xa0',' ').replace(u'\n',' ').replace(' ','')
        
        rowout.append(Decimal( value ))
        convdata3.append(rowout)
    convdata2 = convdata3
                        
    return headers, convdata, convdata2

######################################################################
##### NASA Planetary Fact Sheet (Metric)
# Author/Curator:
# Dr. David R. Williams, dave.williams@nasa.gov
# NSSDCA, Mail Code 690.1
# NASA Goddard Space Flight Center
# Greenbelt, MD 20771
# +1-301-286-1258
# NASA Official: Ed Grayzeck, edwin.j.grayzeck@nasa.gov
# Last Updated: 17 July 2015, DRW
######################################################################

NASAPLNFACTSURL = "http://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html" # NASA Planetary Data URL

def make_NASAPLNFACTS(url=NASAPLNFACTSURL):
    NASAPLNsoup = scraped_BS(url)
    tbl = NASAPLNsoup.soup.find("table")

    tbldata = []
    for row in tbl.find_all('tr'):
        if row.find_all('td') != []:
            tbldata.append( [ele.text for ele in row.find_all('td')] )
    tbldata.pop()
    return tbldata
    
# now, go to toPd.py

######################################################################
##### Rocket & Space Technology
##### Robert A. Braeunig
####   http://www.braeunig.us/space/index.htm
######################################################################

Braeunig_ATMOS_URL="http://www.braeunig.us/space/atmos.htm"

def make_tbls_from_soup(url=Braeunig_ATMOS_URL):
    soup = scraped_BS(url)
    tbls = soup.soup.find_all("table")

    hdrsdata = []
    tblsdata = []
    for tbl in tbls:
        hdrdata = []
        tbldata = []
        for row in tbl.find_all('tr'):
            if row.find_all('td') != []:
                tbldata.append([ele.text for ele in row.find_all('td')])
            elif row.find_all('th') != []:
                hdrdata.append([ele.text for ele in row.find_all('th')])
        tblsdata.append(tbldata)
        hdrsdata.append(hdrdata)
    return tbls, tblsdata, hdrsdata

def make_Braeunig_ATMOS(url=Braeunig_ATMOS_URL):
    tbls, tblsdata, hdrsdata = make_tbls_from_soup(url)
    
    # Physical Properties of U.S. Standard Atmosphere, 1976 in SI Units
    standatm_FOOTNOTE = tblsdata[1].pop()
    standatm_DATA  = [[Decimal(ele.replace(',','')) for ele in row] for row in tblsdata[1]]
    standatm_TITLE = hdrsdata[1][0][0]
    standatm_HDR   = hdrsdata[1][1]
    standatm_dict = {"footnote":standatm_FOOTNOTE,"data":standatm_DATA,"title":standatm_TITLE,"header":standatm_HDR}
    return standatm_dict

# now, go to toPd.py
