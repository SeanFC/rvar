# RVAR 

Make spatially coherent palaeoclimate gridded maps from site-based pollen reconstructions and climate model outputs using a conditioned 3D variational data assimilation method.

## What this code can do 
Puts into a SQL database, clean outputs from:

* Pollen reconstructions from Bartlein et al. (2011)
* The 3rd round of the Palaeoclimate Modelling Intercomparison Project (PMIP3; Braconnot et al. 2012) 
* Climate Research Unit’s CL v2.0 (CRU CL 2.0; New et al. 2002) to use as a modern climate base line

This data is run through a 3D-Variational (3D-Var) data assimilator to find the best a posteriori estimate of the palaeoclimate.
The estimate and its uncertainty is saved to the SQL database and other diagnostic quantities (such as condition number) can be output to text files. 
The method also corrects for changes in atmospheric CO<sub>2</sub> concentration between the palaeo and reference time periods.

## Setup
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

Run `make` from the root directory to compile the code, this should created an executable named `rmin`.


### Install 

An SQL database named `clim` must be available with the pollen reconstructions, PMIP3 and CRU CL 2.0 data; running the setup functions from `main` can help with this.
Also a username and password supplied for the database must be added to the configuration file in the `config` directory.
The configuration file is named the username of the machine which is running the code, found from running `hostname`, with the postfix `.cfg`.

### Run

The assimilator can be run with executable `rmin`. 
Alternatively, the assimilator can be run in the background, tunnelled to an exterior database, by adding the correct configuration options and running the `bac_rmin` script.

## Further Information
Details of the method that this code implements can be found at the [arXiv](https://arxiv.org/abs/1902.04973). 
This code has been used to generate the data for a paper at [Climate of the Past](https://www.clim-past-discuss.net/cp-2019-55/); the data can be found at the [University of Reading Data Archive](http://dx.doi.org/10.17864/1947.206).


## Contact
Sean F. Cleator

Email: s.cleator@surrey.ac.uk

ORCid: [0000-0001-9602-1989](https://orcid.org/0000-0001-9602-1989)


## References

* P. J. Bartlein, S. P. Harrison, S. Brewer, S. Connor, B. A. S. Davis, K. Gajewski, J. Guiot, T. I. Harrison-Prentice, A. Henderson, O. Peyron, I. C. Prentice, M. Scholze, H. Seppä, B. Shuman, S. Sugita, R. S. Thompson, A. E. Viau, J. Williams, and H. Wu. Pollen-based continental climate reconstructions at 6 and 21 ka: a global synthesis. Climate Dynamics, 37:775–802, 2011.
* P. Braconnot, S. P. Harrison, M. Kageyama, P. J. Bartlein, V. Masson-Delmotte, A. Abe-Ouchi, B. Otto-Bliesner, and Y. Zhao. Evaluation of climate models using palaeoclimatic data. Nature Climate Change, 2:417–424, 2012.
* M. New, D. Lister, M. Hulme, and I. Makin. A high-resolution data set for surface climate over
global land areas. Climate Research, 21:1–25, 2002.
