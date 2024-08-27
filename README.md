# VPRetina
# Linux installation
This code has been tested on Linux (Ubuntu 22.04) with Python 3.10.<br>
Clone the repository:
```
   git clone git@github.com:remihndz/VPRetina.git
   cd VPRetina
```

## Installation of C++ libraries

Build and install the VItA libray (more details in [this README](https://github.com/remihndz/VItA/blob/040e2d0267054fa9b4442b891862555986d99e3c/README.md)):
```
cd src/VItA/
mkdir vita_build
cd vita_build
ccmake ../vita_source
```
In the ccmake interface, press "c" to configure, make sure is set to current directory (vita_build) by setting "CMAKE_INSTALL_PREFIX = .", press "c" again, then "g" to generate makefiles. Then type:
```
make
make install
```
Build the main c++ files:
```
cd ../../SVC/CCO/cpp/
mkdir build
cd build
ccmake ../
```
In the ccmake interface, press "c" to configure, make sure is set to current directory (vita_build) by setting "CMAKE_INSTALL_PREFIX = .", press "c" again, then "g" to generate makefiles. Then type:
```
make
mv svc.py ../../
mv _svc.so ../../
```

## Python virtual environment
Install the necessary python packages in a new conda environment:
```
conda create --name <env_name> python=3.10
conda activate <env_name>
python -m pip install -r requirements.txt
```
or using pip (assuming you have python version >=3.10 installed):
```
python3.1x -m pip install -r requirements.txt
```

# Usage
Virtual vasculatures can be generated using the [scripts](./scripts) provided.<br>
Use `python3 RunSims.py` to create `n` vasculatures.<br>
This will create several folders and files:
  - Input parameters and hyperparameters are contained in `PopulationParameters.csv`. Each row corresponds to one simulation (i.e., one virtual retina) indexed by the column `sim`.
  - [In SSM](./SSM/), you can find the graphs resulting from the Statistical Shape Model. There are *two* files for each virtual retina: `sim_{i}_artery.cco` and `sim_{i}_vein.cco`.
  - [In Coarse](./Coarse/), you can find the graphs resulting from the first stage of CCO growth. There are *six* files for each virtual retina: `sim_{i}_{artery or vein}.cco` which contain the graphs resulting from the CCO, `sim_{i}_{artery or vein}.root` which is a copy of the file found in [SSM](SSM/) and `sim_{i}_{artery or vein}.conf` which contains the hyperparameters for this stage.
  - [In Macula](./Macula/), you can find the results from the second stage of the CCO. The naming convention is the same as in `Coarse`. The `.root` files are copies of the files in `Coarse`.
  - [In AVTrees](./AVTrees/), you can find the merged arteries-veins graph resulting from the second stage of CCO growth (following the same naming convention, i.e., `sim_{i}_AV.cco`).

Use `python3 GlobalMetrics.py AVTrees/*_AV.cco` to compute the morphological metrics for a list of networks. The results can be found in `OCTAMetrics.csv`, where<br>

Use `python3 haemodynamics.py /path/to/PopulationParameters.csv` to run haemodynamics simulations for all the simulations in `PopulationParameters.csv`.<br>

There are also parallel versions of these scripts for faster execution.

# Publications

Rémi J. Hernandez, Savita Madhusudhan, Yalin Zheng, Wahbi K. El-Bouri; Linking Vascular Structure and Function: Image-Based Virtual Populations of the Retina. Invest. Ophthalmol. Vis. Sci. 2024;65(4):40. https://doi.org/10.1167/iovs.65.4.40.

