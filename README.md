# VPRetina

# Linux installation

Clone the repository:
```
   git clone https://github.com/remihndz/VPRetina.git
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
```

## Python virtual environment

TODO:
	- Figure out what are the necessary packages
	- Make sure code is commented and as clean as possible
	- Create the requirement.txt


# Publications

Rémi J. Hernandez, Savita Madhusudhan, Yalin Zheng, Wahbi K. El-Bouri; Linking Vascular Structure and Function: Image-Based Virtual Populations of the Retina. Invest. Ophthalmol. Vis. Sci. 2024;65(4):40. https://doi.org/10.1167/iovs.65.4.40.

