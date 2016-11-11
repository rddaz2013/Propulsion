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
for other Linux distributions/flavors, the same libraries have different names for different package managers and some libraries were already installed with the "stock" OS and some aren't (as I found in my situation.  For example, Cantera's mainpage, for Ubuntu/Debian installation (compilation), it's neglected that `boost` is already installed (which I found wasn't for Fedora 23 Workstation Linux).

### Installation Prerequisites for **Fedora 23 Workstation Linux** (make sure to do these `dnf install`s and installation with `scons` will go more smoothly).

I found that you can't get away from `dnf install` on an *administrator* account - be sure to be on a `sudo` or `admin` account to be able to do `dnf install`s.  Also, I found that compiling Cantera had to be done on a `sudo`-enabled or administrator account, in particular, access is needed to be granted to accessing root directories such as `/opt/`, etc. (more on that later).

Also, in general, you'd want to **install the developer version** of the libraries as well, usually suffixed with `-devel`, mostly because the header files will be placed in the right `/usr/*` subdirectory so to be included in the system (when compiling C++ files or installing).  

- **`scons`** - be sure to install `scons` - it seems like there is a push to use scons, a Python program, for installation and (package) compilation, as opposed to (old-school) CMake, or Make.
- **`boost`** - *Boost* is free peer-reviewed portable C++ source libraries.
```
sudo dnf install boost.x86_64
sudo dnf install boost-devel.x86_64
```
- **`lapack`** - *`lapack`, Linear Algebra PACkage*.  Don't take it for granted that `lapack` is already installed (I had to troubleshoot this myself, beyond the Cantera main page documentation, and find where it is).  I had to install it because I found it was missing through the Cantera `scons build`  
```
dnf list lapack*  # find lapack in dnf
sudo dnf install lapack.x86_64
sudo dnf install lapack-devel.x86_64
```
- **`blas`** - *`blas`, Basic Linear Algebra Subprograms*.  Don't take it for granted that `blas` is already installed (I had to troubleshoot this myself, beyond the Cantera main page documentation, and find where it is).  I had to install it because I found it was missing through the Cantera `scons build` 
```
dnf list blas*  # find blas in dnf
sudo dnf install blas.x86_64
sudo dnf install blas-devel.x86_64
```
- **`python-devel`** - Following the spirit of how you'd want to install the developer's version of the library concurrent with the library itself, in that you'd want the headers and symbolic links to be installed and saved onto the respective root `/usr/*` subdirectories (so that your system will know how to include the files), you'd want to install the Python developer's libraries.
```
sudo dnf install python-devel
```
On this note, for Fedora Linux, I *did not* find with `dnf list` `python-numpy` nor `python-numpy-dev` which, supposedly, is found in Ubuntu/Debian - this is an example of how Fedora/CentOS/RedHat package manager is different from Ubuntu/Debian.




## Clean install, from `git clone` to `scons install`

- `git clone https://github.com/Cantera/cantera.git`
```
git clone https://github.com/Cantera/cantera.git
```
-`scons build -j12`
```
scons build -j12
```
`scons build` by itself is ok; I added the flag `-j12` (correct me if I'm wrong) to optimize the compilation on **12** cores.  So if you're on a quad-core CPU processor, then you'd do `-j4`.  
-`sudo scons install`
```
sudo scons install
```
There's no getting around not using sudo for scons install.  

## Troubleshooting installation/(installation) errors that pop up

- `fatal error: Python.h: No such file or directory`
```
fatal error: Python.h: No such file or directory
scons: *** [build/temp-py/_cantera2.os] Error 1
```
I found that I had to `dnf install` `python-devel` to get the header files installed onto the appropriate `/usr/*` root subdirectories.  
- `scons: *** [/usr/local/include/cantera/Edge.h] /usr/local/include/cantera/Edge.h: Permission denied`  
Do `sudo scons install`
- `error: could not create `/usr/local/lib64/python2.7': Permission denied`  
Do `sudo scons install`

## Troubleshooting C++ compilation/(C++ compilation) errors that pop up

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



## Images gallery (that may help you with your installation process; it can be daunting)

```
dnf list boost-*
```
![dnf list boost](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/boostdevel01Screenshot%20from%202016-11-11%2000-26-42.png)

```
sudo dnf install boost-devel.x86_64

```

![sudo dnf install boost-devel.x86_64](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/boostdevel01Screenshot%20from%202016-11-11%2000-27-14.png)


```
dnf list lapack*  # find lapack in dnf
sudo dnf install lapack-devel.x86_64
```

![sudo dnf install lapack-devel.x86_64](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/lapackblasScreenshot%20from%202016-11-11%2001-06-12.png)

```
sudo dnf install python-devel
```

![sudo dnf install python-devel](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/python-develinstallsconsbuildScreenshot%20from%202016-11-11%2002-09-23.png)



```
git clone https://github.com/Cantera/cantera.git
```

![git clone](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/gitclonecanteraScreenshot%20from%202016-11-11%2002-02-34.png)


```
fatal error: Python.h: No such file or directory
scons: *** [build/temp-py/_cantera2.os] Error 1
```
![fatal error: Python.h: No such file or directory](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/scaninstallfailpythonherrorScreenshot%20from%202016-11-11%2002-07-47.png)

### Successful installation/compilation (what we want, what it should look like)

```
scons build
```

![sconsbuildsuccess](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/sconsbuildsuccessfulScreenshot%20from%202016-11-11%2002-10-46.png)

```
sudo scons install
```
There's no way, I found, of getting away from having to use `sudo` for scons install - you'll have to be on a sudo enabled or administrator account logged in.  

It troubleshoots
```
scons: *** [/usr/local/include/cantera/Edge.h] /usr/local/include/cantera/Edge.h: Permission denied
error: could not create `/usr/local/lib64/python2.7': Permission denied
```

![sudo scons install](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/sconsinstallScreenshot%20from%202016-11-11%2002-05-54.png)

![sudo scons install](https://raw.githubusercontent.com/ernestyalumni/Propulsion/master/cantera_stuff/cantera_install_tips/images/sconsinstallScreenshot%20from%202016-11-11%2002-14-43.png)