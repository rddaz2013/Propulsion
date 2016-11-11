# Cantera Installation Tips

Installing Cantera on Fedora Linux, straight, directly from the github repository, all the way to being compiled with `scons`, was nontrivial, mostly because of the installation prerequisites, which, in retrospect, can be easily installed *if* one knows what they are with respect to what it is in terms of Fedora/CentOS/RedHat `dnf`.  

| codename                  | directory      | reference webpage (if any) | Description  |
| ------------------------- | :------------- | :------------------------- | :----------: | 
| `cantera_install_success` | `./`           | *None*                     | A *verbose*, but complete Terminal *log* of cantera installation on *Fedora Workstation 23 Linux*, from `git clone`, cloning the githb repository for cantera, directly, all the way to a successful `scons install`. |
| `ClassThermoPhaseExam.cpp` | `./`          | [Computing Thermodynamic Properties, Class ThermoPhase, Cantera C++ Interface User's Guide](http://www.cantera.org/docs/sphinx/html/cxx-guide/thermo.html#example-program) | Simple, complete program creates object representing gas mixture and prints its temperature |
| `chemeqex.cpp`            | `./`           | [Chemical Equilibrium Example Program, Cantera C++ INterface User's Guide](http://www.cantera.org/docs/sphinx/html/cxx-guide/equil-example.html) | `equilibrate` method called to set gas to state of chemical equilibrium, holding temperature and pressure fixed. |

## Installation Prerequisites, ala Fedora Linux, Fedora/CentOS/RedHat `dnf`

While [Cantera mainpage's Cantera Compilation Guide](http://www.cantera.org/docs/sphinx/html/compiling.html) gave the packages in terms of Ubuntu/Debian's package manager:
```
g++ python scons libboost-all-dev libsundials-serial-dev
```
and for the python module
```
cython python-dev python-numpy python-numpy-dev
```



## Troubleshooting installation/(installation) errors that pop up

- `collect2: error: ld returned 1 exit status`
From this page that I found on Google search:
cf. [What does “collect2: error: ld returned 1 exit status” mean? stackoverflow](http://stackoverflow.com/questions/27272525/what-does-collect2-error-ld-returned-1-exit-status-mean)

I realized that I needed to include the Cantera library in this way:
```
-lcantera
```
when compiling with g++.  


- `/usr/lib64/libpthread.so.0: error adding symbols: DSO missing from command line`

I Google searched for this webpage:
cf. [“error adding symbols: DSO missing from command line” while compiling g13-driver, ask ubuntu](http://askubuntu.com/questions/521706/error-adding-symbols-dso-missing-from-command-line-while-compiling-g13-driver)

From this page, I saw the use of the line `LIBS = -lusb-1.0 -l pthread`, and the idea of using the flag `-l pthread` ended up being the solution.  





