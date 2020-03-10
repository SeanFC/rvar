# RVAR 

Version: 0.2.0 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3445165.svg)](https://doi.org/10.5281/zenodo.3445165)

Make spatially coherent gridded maps of the palaeoclimate by combining site-based pollen reconstructions and climate model outputs using a conditioned 3D variational data assimilation method. Implements the method shown in [this paper](https://doi.org/10.1029/2019MS001630).

## What this code can do 
Puts into a SQL database, cleaned outputs from:

* Pollen reconstructions from Bartlein et al. (2011)
* The 3rd round of the Palaeoclimate Modelling Intercomparison Project (PMIP3; Braconnot et al. 2012) 
* Climate Research Unit’s CL v2.0 (CRU CL 2.0; New et al. 2002) 

This data is run through a 3D-Variational (3D-Var) data assimilator to find the best a posteriori estimate of the palaeoclimate.
The estimate and its first order variance is saved to the SQL database and other diagnostic quantities (such as condition number) can be output to text files. 
The method also corrects for changes in atmospheric CO<sub>2</sub> concentration between palaeo and reference time periods.

## Setup to make reconstructions
Most of the setup must be user configured and was tested on a machine running archlinux. 

### Compilation 
First install the prerequisite packages:

* [mysql-connector-c++](https://aur.archlinux.org/packages/mysql-connector-c%2B%2B/)
* [netcdf_cxx](https://www.archlinux.org/packages/community/x86_64/netcdf-cxx/)
* [liblbfgs](https://aur.archlinux.org/packages/liblbfgs/)
* [gcc-fortran](https://www.archlinux.org/packages/core/x86_64/gcc-fortran/)
* [openblas](https://www.archlinux.org/packages/community-testing/x86_64/openblas/)
* [cblas](https://www.archlinux.org/packages/extra/x86_64/cblas/)
* [lapacke](https://www.archlinux.org/packages/extra/x86_64/lapacke/)

Run `make configure` to create machine specific files.

Run `make` from the root directory to compile the code, this creates an executable named `rmin`.

### Install 

An SQL database named `clim` must be available with the pollen reconstructions, PMIP3 and CRU CL 2.0 data; running the setup functions from `src/main.cpp` can help with this.
A username and password supplied for the database must be added to the configuration file in the `config` directory.
The configuration file is named with the host name of the machine which is running the code (found from running `hostname`) with the postfix `.cfg`.

### Run

The assimilator can be run with the executable `rmin`. 
Alternatively, the assimilator can be run in the background and tunnelled to an exterior database by running the `bac_rmin` script.

## Further Information
Details of the method that this code implements have been published in the [Journal of Advances in Modeling Earth Systems](https://doi.org/10.1029/2019MS001630). 
This code has been used to generate the data for a paper at [Climate of the Past](https://www.clim-past-discuss.net/cp-2019-55/); the data can be found at the [University of Reading Data Archive](http://dx.doi.org/10.17864/1947.206).

## Contact
Sean F. Cleator

Email: s.f.cleator@surrey.ac.uk

ORCid: [0000-0001-9602-1989](https://orcid.org/0000-0001-9602-1989)

## References

* P. J. Bartlein, S. P. Harrison, S. Brewer, S. Connor, B. A. S. Davis, K. Gajewski, J. Guiot, T. I. Harrison-Prentice, A. Henderson, O. Peyron, I. C. Prentice, M. Scholze, H. Seppä, B. Shuman, S. Sugita, R. S. Thompson, A. E. Viau, J. Williams, and H. Wu. Pollen-based continental climate reconstructions at 6 and 21 ka: a global synthesis. Climate Dynamics, 37:775–802, 2011.
* P. Braconnot, S. P. Harrison, M. Kageyama, P. J. Bartlein, V. Masson-Delmotte, A. Abe-Ouchi, B. Otto-Bliesner, and Y. Zhao. Evaluation of climate models using palaeoclimatic data. Nature Climate Change, 2:417–424, 2012.
* M. New, D. Lister, M. Hulme, and I. Makin. A high-resolution data set for surface climate over
global land areas. Climate Research, 21:1–25, 2002.
