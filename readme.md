# PFP-Thickening

This project is a demo to achieve to algorithm of PFP-Thickening.

[![c++20](https://img.shields.io/badge/standard-C++20-blue.svg?style=flat&logo=c%2B%2B)](https://isocpp.org)

![BuildExamplesLinux](https://github.com/rainbowwing/Thickening2/workflows/CMake/badge.svg)

[![License](https://img.shields.io/badge/License-BSD_3--Clause-orange.svg)](https://github.com/rainbowwing/Thickening2/blob/main/LICENSE)
## 📚 Documentation

This code is a demo of the paper [PFP-Thickening: a Parallel Feature-preserving Mesh Offsetting Approach
to Create Variable-thickness Solid for Additive manufacturing]().

It can be use to calculate the offset result of a mesh.

This demo only have the code which working in manifold meshes.

Use the mesh kernel which built by [intelligent visual modeling & simulation lab](https://igame.hdu.edu.cn) in Hangzhou
Dianzi University

This code can be built in Linux and MacOS.

## 📄 Dependencies

[Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page)

[Osqp-eigen](https://github.com/robotology/osqp-eigen#osqp-eigen)

[Cgal](https://www.cgal.org)

[Gflags](https://github.com/gflags/gflags)

[Vcglib](https://github.com/cnr-isti-vclab/vcglib)

## 🛠️ Build

### Linux
The version of g++ must later than 12.0. We test the g++-7.2, it can't compile successfully.

cgal version must be larger than 5.5.
```bash
sudo su root
apt update
apt upgrade
apt install libcgal-dev 
apt install libeigen3-dev 
apt install libgflags-dev 
git clone --recursive https://github.com/osqp/osqp
cd osqp
mkdir build
cd build
cmake .. -DBUILD_SHARED_LIBS=ON
make
sudo make install
cd ../..
git clone https://github.com/robotology/osqp-eigen.git
cd osqp-eigen
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt ../
make
make install
git clone git@github.com:TsukiMiyabiLake/Thickening2.git
cd Thickening2
#git clone https://github.com/cnr-isti-vclab/vcglib # add the submodule vcglib
mkdir build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make

```

### MacOS
The version of macOS must later than 13.0.

cgal version must be larger than 5.5
```bash
sudo su root
brew install cgal
brew install eigen
brew install gflags
git clone --recursive https://github.com/osqp/osqp
cd osqp
mkdir build
cd build
cmake .. -DBUILD_SHARED_LIBS=ON
make
sudo make install
cd ../..
git clone https://github.com/robotology/osqp-eigen.git
cd osqp-eigen
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt ../
make
make install
cd ../..
git clone git@github.com:TsukiMiyabiLake/Thickening2.git
cd Thickening2
cd TetWild
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
cd ../..
cd fTetWild
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
mkdir build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make

```

## Usage:

```bash

USAGE: Thicken -f file_name [options]              file_name is *.obj or *.obj2      

options:
  -m={1|2}        representing result mode which value can be chose in 1 and 2. 
                    mode 1 get the result with do remeshing;
                    mode 2 get the result without  do remeshing;

  -d=<num>          this number is a double represent the expect length of each facet in running invariable thickening.
                    please do not set is very large, because the program will cost large time. 
                    
  -l=<num>          This number is a double.
                    which value indicates how many times the maximum offset distance is the ideal offset distance.
                    You can set it is 1.5 or do not set it.
                    
  -s=<num>          This number is a double.
                    which value indicates how many times the minimum offset distance is the ideal offset distance.
                    You can set it is 1.0 or do not set it. You need to ensure that the value of -s is less than -l.
                    
  -t=<num>          This number is a integer represent thread num please set this value depend the number of cpu cores of you device.
  
  -i={1|2}          representing running mode which value can be chose in 1,2. 
                    1 is offsetting to the outside of the mesh;
                    2 is offsetting to the inside of the mesh.
                    
  -e=<num>          This number is a double means the eps.
                    When the distance of two points is smaller than eps, we will regard these two point as coinciding.
                    You can set it is 0.0001 or do not set it.

                    

example:
./Thicken  -f ../data/mechanical04.obj -m=1 -i=2 -t=8 -d=0.6        
./Thicken  -f ../data/tet.obj2 -m=1 -i=1 -t=8 
```


### The input file format

If you want to do invariable thicken. You can use the obj file than set the expect distance by setting the argument -d.

If you want to do variable thicken, you should use the obj2 format.

obj2 format which is same as *.obj. But in the end line of each face, it needs to add a double representing the ideal offsetting value.


example: tet.obj2
```text
v 0 0 0
v 1 0 0
v 0 1 0
v -1 -1 1
f 1 3 2 0.05
f 1 4 3 0.04
f 2 4 1 0.03
f 2 3 4 0.02

```

### NOTE:

cgal version must be larger than 5.5

If you want to do not use this program for thickening with long offset distance, this program will not have a fast speed.
You need to wait a lot of time before the program finishes running.

### KNOWN ISSUES:

In extremely rare cases, there will be problems with cgal allocating memory under multiple threads. On devices with sufficient memory, this phenomenon does not occur frequently.
So it can be solved by re-running at present.


## 📝 License
Materials in this repository are distributed under the following license:

> All software is licensed under the BSD 3-Clause License. See [LICENSE](https://github.com/rainbowwing/Thickening2/blob/main/LICENSE) file for details.