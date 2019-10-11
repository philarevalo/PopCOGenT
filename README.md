# PopCOGenT

# Purpose
Identifying microbial populations using networks of horizontal gene transfer. For more detailed information on the method, you can read our paper here: 

Arevalo P., VanInsberghe D., Elsherbini J., Gore J., Polz M.F. (2019). [A reverse ecology approachbased on a biological definition of microbial populations](https://linkinghub.elsevier.com/retrieve/pii/S0092867419307366). Cell, 178(4).doi:10.1016/j.cell.2019.06.03

# Dependencies

## All modules
* A linux-based system (N.B., we are working on a way to get all the dependencies working properly on OSX, but as of now OSX is not a supported operating system).
* [Miniconda with python 3.7](https://docs.conda.io/en/latest/miniconda.html)

The required python and (most) R packages can be installed by creating a conda environment with the included `PopCOGenT.yml` file as follows:

`conda config --set restore_free_channel true`

`conda env create -f PopCOGenT.yml`

The first command is necessary as we use some versions of packages that are only availble on the `free` conda channel. The `PopCOGenT.yml` file reflects the versions we used for our publication and at this time we cannot guarantee forward compatability with new versions of the dependencies. 

## PopCOGenT

* [mugsy version 1.2.3](http://mugsy.sourceforge.net/)
* [Infomap version 0.18.3, will likely work with 0.20](https://www.mapequation.org/code.html#Installation)

## Flexible genome sweep identification
* [mmseqs2 version vNatBiotech](https://github.com/soedinglab/MMseqs2)
* [muscle 3.8.31](https://www.drive5.com/muscle/)

## Core genome sweep identification
* [phyml version 3.1](http://www.atgc-montpellier.fr/phyml/versions.php)
* [mugsy version 1.2.3](http://mugsy.sourceforge.net/)
* The `ape` R package. To install, please follow the instructions under "All modules." Then, activate the environment (`source activate PopCOGenT`). Finally, run the `Rscript install_ape.R` from the Core genome sweep identification source directory.

# Usage

Instructions for the usage of each module are provided in each module's source code directory. Populations are identified with PopCOGenT and then core and flexible genome sweeps differentiating the most closely related populations (i.e., populations that are still connected by some gene flow but are still characterized by significantly more gene flow within the population than between populations) can be identified with the core and flexible genome sweep modules. We do not recommend running the sweep modules to identify sweeps between populations that are completely disconnected by gene flow (i.e., species). 
