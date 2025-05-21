# MPSCode

The current version of the code is compiled using g++.

The libraries Boost and Ula are needed, as well as the installation for the following:

https://www.boost.org/users/download/
sudo apt install gfortran
sudo apt install liblapack-dev
sudo apt install libblas-dev

It's run using a Make file where the root directory must be changed 

IDFLAGS     = -I. -I/home/seb/Desktop/MPO/libs/boost_1_88_0 -I/home/seb/Desktop/MPO/libs/ula/include -I/home/seb/Desktop/MPO/libs/ula/include/ula 

where -I/home/sebs/ as well as where the files (and standalone libraries) is up to the user. 

