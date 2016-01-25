## tut2.py
## tut2.m implemented in Python for cantera 
## cf. http://www.cantera.org/docs/sphinx/html/matlab/tutorials/tut2.html
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160124
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
# Tutorial 2: Working with input files
#
#   Topics:
#     - using functions 'importPhase' and 'importInterface'
#     - input files distributed with Cantera
#     - the Cantera search path
#     - CTML files
#     - converting from CK format
#
import cantera as ct
import timeit

# In the last tutorial, we used function GRI30 to create an object
# that models an ideal gas mixture with the species and reactions of
# GRI-Mech 3.0. Another way to do this is shown here, with statements
# added to measure how long this takes:

gas1 = ct.Solution('gri30.cti','gri30')
msg = timeit.timeit(stmt="gas1 = ct.Solution('gri30.cti','gri30')", setup="import cantera as ct", number=1)
print "time to create gas1: %f:" % msg

# class Solution constructs an object representing a phase of
# matter by reading in attributes of the phase from a file, which in
# this case is 'gri30.cti'. This file contains several phase
# specifications; the one we want here is 'gri30', which is specified
# by the second argument.  This file contains a complete specification
# of the GRI-Mech 3.0 reaction mechanism, including element data
# (name, atomic weight), species data (name, elemental composition,
# coefficients to compute thermodynamic and transport properties), and
# reaction data (stoichiometry, rate coefficient parameters). The file
# is written in a format understood by Cantera, which is described in
# the document "Defining Phases and Interfaces."

# EY : 20160124 according to 
# https://github.com/Cantera/cantera/blob/5926d2db7c0d4919b75ee50828b0adab4e691a51/interfaces/cython/cantera/base.pyx
# the second argument 'gri30' is the argument for phaseid, 
#  with the default being phaseid=''
# the first name, 'gri30.cti' is the argument for infile
#  with the default being infile=''

# On some systems, processing long CTI files like gri30.cti can be a
# little slow. For example, using a typical laptop computer running
# Windows 2000, the statement above takes more than 4 s, while on a
# Mac Powerbook G4 of similar CPU speed it takes only 0.3 s. In any
# case, running it again takes much less time, because Cantera
# 'remembers' files it has already processed and doesn't need to read
# them in again:

msg = timeit.timeit(stmt="gas1b = ct.Solution('gri30.cti','gri30')", setup="import cantera as ct", number=1)
print "time to create gas1: %f:" % msg

# CTI files distributed with Cantera
#-----------------------------------

# Several reaction mechanism files in this format are included in the
# Cantera distribution, including ones that model high-temperature
# air, a hydrogen/oxygen reaction mechanism, and a few surface
# reaction mechanisms. Under Windows, these files may be located in
# 'C:\Program Files\Common Files\Cantera', or in 'C:\cantera\data',
# depending on how you installed Cantera and the options you
# specified.  On a unix/linux/Mac OSX machine, they are usually kept
# in the 'data' subdirectory within the Cantera installation
# directory.

# If for some reason Cantera has difficulty finding where these files
# are on your system, set environment variable CANTERA_DATA to the
# directory where they are located. Alternatively, you can call function
# add_directory to add a directory to the Cantera search path:

# ct.add_directory('/usr/local/cantera/my_data_files')

# EY : 20160125 or just work in the directory that you run python in, 
#  in your Terminal Command Prompt

# Cantera input files are plain text files, and can be created with
# any text editor. See the document 'Defining Phases and Interfaces'
# for more information.


# Importing multiple phases or interfaces
# ---------------------------------------

# A Cantera input file may contain more than one phase specification,
# or may contain specifications of interfaces (surfaces). Here we
# import definitions of two bulk phases and the interface between them
# from file diamond.cti:

gas2 = ct.Solution('diamond.cti','gas')          # a gas
diamond = ct.Solution('diamond.cti','diamond')   # bulk diamond
diamond_surf = ct.Interface('diamond.cti','diamond_100',[gas2,diamond])

# CTML files
# ----------

# Note that when Cantera reads a .cti input file, wherever it is
# located, it always writes a file of the same name but with extension
# .xml *in the local directory*. If you happen to have some other file
# by that name, it will be overwritten. Once the XML file is created,
# you can use it instead of the .cti file, which will result in
# somewhat faster startup.

gas4 = ct.Solution('gri30.xml','gri30')

# Interfaces can be imported from XML files too.
diamond_surf2 = ct.Interface('diamond.xml','diamond_100',[gas2, diamond])

# Converting CK-format files
# --------------------------

# Many existing reaction mechanism files are in "CK format," by which
# we mean the input file format developed for use with the Chemkin-II
# software package. [See R. J. Kee, F. M. Rupley, and J. A. Miller,
# Sandia National Laboratories Report SAND89-8009 (1989).]

# Cantera comes with a converter utility program 'ck2cti' (or
# 'ck2cti.exe') that converts CK format into Cantera format. This
# program should be run from the command line first to convert any CK
# files you plan to use into Cantera format. This utility program can
# also be downloaded from the Cantera User's Group web site.
#
# Here's an example of how to use it:
#
# ck2cti -i mech.inp -t therm.dat -tr tran.dat -id mymech 
#
# EY : 20160125 I want to investigate further into which file is for which input, what the flags mean, and if mymech is the output or not

######################################################################
#   end of tutorial 2
######################################################################
