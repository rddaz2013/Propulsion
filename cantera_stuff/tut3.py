## tut3.py
## tut3.m implemented in Python for cantera 
## cf. http://www.cantera.org/docs/sphinx/html/matlab/tutorials/tut3.html
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
# Tutorial 3:   Getting Help
#
import cantera as ct

# Suppose you have created a Cantera object and want to know what
# methods are available for it, and get help on using the methods. 

g = ct.Solution("gri30.xml")

# The first thing you need to know is the class object g
# belongs to. Type:

type(g)

# This tells you that g belongs to a class called 'Solution'. To find
# the methods for this class, type

dir(ct.Solution)

# This command returns all method names as a Python list. 

# A long list is printed.  Some methods are 
# inherited from other classes. For example, variable P and 
# method set_unnormalized_mass_fractions are
# inherited from a class 'ThermoPhase'. Don't be concerned at this
# point about what these base classes are - we'll come back to them
# later.

# Now that you see what methods are available, you can type 
# 'help(<method_name>)' to print help text for any method. For example,

help(ct.Solution.P)
help(ct.Solution.set_unnormalized_mass_fractions)
help(ct.Solution.net_rates_of_progress)

# For help on how to construct objects of a given class, type 
# 'help(<classname>)'

help(ct.Solution)

# Now that you know how to get help when you need it, you can
# explore using the Cantera Toolbox on your own. But there are a
# few more useful things to know, which are described in the next
# few tutorials.

#################################################################
#   end of tutorial 3
#################################################################


