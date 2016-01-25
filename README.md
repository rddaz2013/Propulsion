# Propulsion
Propulsion - Numerical recipes for (rocket) propulsion, including notes and solutions in LaTeX

## Cantera

### `cantera_stuff` - examples (of usage) and implementations of Cantera
First, the folder [cantera_stuff](https://github.com/ernestyalumni/Propulsion/tree/master/cantera_stuff) contains an implementation in Python of the (useful) Matlab tutorial .m files/examples/ and (some of the) Matlab examples on the [Index of Examples](http://www.cantera.org/docs/sphinx/html/matlab/examples.html) of *Cantera Matlab Toolbox* (but it's now in **Python**). *See* the Cantera Matlab Toolbox examples page http://www.cantera.org/docs/sphinx/html/matlab/examples.html and compare it with the files in [cantera_stuff](https://github.com/ernestyalumni/Propulsion/tree/master/cantera_stuff).

#### (cantera) Tutorials
Of note are the tutorial files for *Cantera* (which I recommend that one works through):
- `tut1.py`
- `tut2.py`
- `tut3.py`
- `tut4.py`
- `tut5.py`
- `tut6.py`
- `tut7.py`
and 

#### (cantera) Examples
- `equil.py`

#### My own (EY's) Examples
- `LOXmeth_eq.py` 

`LOXmeth_eq.py` calculates, as a function of oxidizer/fuel O/F *mass* ratio, the adiabatic flame temperature, equilibrium molecular composition, mean molecular weight, ratio of specific heats, and characteristic velocity, of the combustion of oxidizer oxygen (O2) and methane (CH3).  This function is further generalized (called `equil_general`) to make other species be the oxidizer and fuel, such as dinitrogen tetraoxide (N2O4) and **hydrazine** (N2H4).  

## Physique
Physique is a (small) Python package for web scraping physical constants data and "data wrangling" or "cleaning" the data to use for Python pandas, as a pandas DataFrame object.
