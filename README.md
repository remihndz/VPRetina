# VPRetina

# Linux installation

Clone the repository:

   git clone https://github.com/remihndz/VPRetina.git
   cd VPRetina

Build and install the VItA libray (more details in [this README](https://github.com/remihndz/VItA/blob/040e2d0267054fa9b4442b891862555986d99e3c/README.md)):

   cd src/VItA/
   mkdir vita_build
   cd vita_build
   ccmake ../vita_source

In the ccmake interface, press "c" to configure, make sure is set to current directory (vita_build) by setting "CMAKE_INSTALL_PREFIX = .", press "c" again, then "g" to generate makefiles. Then type:

   make
   make install

Build the main c++ files:

   cd ../../SVC/CCO/cpp/
   mkdir build
   cd build
   ccmake ../

In the ccmake interface, press "c" to configure, make sure is set to current directory (vita_build) by setting "CMAKE_INSTALL_PREFIX = .", press "c" again, then "g" to generate makefiles. Then type:

   make


