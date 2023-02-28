# PFP-Thickening

This project is a demo to achieve to algorithm of PFP-Thickening.

[//]: # (|  General  | [![c++14]&#40;https://img.shields.io/badge/standard-C++14-blue.svg?style=flat&logo=c%2B%2B&#41;]&#40;https://isocpp.org&#41; [![License]&#40;https://img.shields.io/badge/License-BSD_3--Clause-orange.svg&#41;]&#40;https://github.com/robotology/osqp-eigen/blob/master/LICENSE&#41; |)

[//]: # (| :-------: | :----------------------------------------------------------: |)

[//]: # (| **CI/CD** | [![Codacy Badge]&#40;https://app.codacy.com/project/badge/Grade/a18710c10f1c4df19bc2759fd50e9cf5&#41;]&#40;https://www.codacy.com/gh/robotology/osqp-eigen/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=robotology/osqp-eigen&amp;utm_campaign=Badge_Grade&#41; [![CI]&#40;https://github.com/robotology/osqp-eigen/workflows/C++%20CI%20Workflow/badge.svg&#41;]&#40;https://github.com/robotology/osqp-eigen/workflows/C++%20CI%20Workflow/badge.svg&#41; [![Azure]&#40;https://dev.azure.com/conda-forge/feedstock-builds/_apis/build/status/osqp-eigen-feedstock?branchName=master&#41;]&#40;https://dev.azure.com/conda-forge/feedstock-builds/_build/results?buildId=341091&view=results&#41; |)

[//]: # (| **conda** | [![Conda Recipe]&#40;https://img.shields.io/badge/recipe-osqp--eigen-green.svg&#41;]&#40;https://anaconda.org/conda-forge/osqp-eigen&#41;  [![Conda Downloads]&#40;https://img.shields.io/conda/dn/conda-forge/osqp-eigen.svg&#41;]&#40;https://anaconda.org/conda-forge/osqp-eigen&#41;  [![Conda Version]&#40;https://img.shields.io/conda/vn/conda-forge/osqp-eigen.svg&#41;]&#40;https://anaconda.org/conda-forge/osqp-eigen&#41;  [![Conda Platforms]&#40;https://img.shields.io/conda/pn/conda-forge/osqp-eigen.svg&#41;]&#40;https://anaconda.org/conda-forge/osqp-eigen&#41; |)

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
git clone https://github.com/cnr-isti-vclab/vcglib # add the submodule vcglib
mkdir build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make

```

### MacOS
The version of macOS must later than 13.0.

cgal version must be larger than 5.5
```bash
sudo su root
apt update
apt upgrade
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
git clone git@github.com:TsukiMiyabiLake/Thickening2.git
cd Thickening2
git clone https://github.com/cnr-isti-vclab/vcglib # add the submodule vcglib 
mkdir build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make

```

## Usage:

```bash

USAGE: Thicken2 -f file_name [options]              file_name is *.obj or *.obj2      

options:
  -m={1|2|3}        representing result mode which value can be chose in 1,2 and 3. 
                    mode 1 get the result without mesh and remeshing;
                    mode 2 get the result with building mesh;
                    mode 3 get the result with building mesh and remeshing.
  -d=<num>          this number is a double represent the expect length of each facet in running invariable thickening.
                    Limited in 0.1~1.5. This number does not represent an absolute distance.
                    Example -d=0.5, indicating that the offset distance is 0.5 times the average mesh edge length.
  -l=<num>          This number is a double.
                    which value indicates how many times the maximum offset distance is the ideal offset distance limited in 1.5~2.7.
                    You can set it is 2.0 or do not set it.
  -g=<num>          This number is a double which means the length of edge length of each cell in grid.
                    If you can not calculate a length with better performance, it can be passed.
                    Then it will use the default value.
  -t=<num>          This number is a integer represent thread num please set this value depend the cpu of you device.
  -i={1|2}          representing running mode which value can be chose in 1,2. 
                    1 is offsetting to the outside of the mesh;
                    2 is offsetting to the inside of the mesh.
  -e=<num>          This number is a double means the eps.
                    When the distance of two points is smaller than eps, we will regard these two point as coinciding.
                    You can set it is 0.00001 or do not set it.
  -s={0|1}          This value is a boolean, when it is 1, the program will skip some cell which is most likely useless.
                    It can improve performance, but may cause holes in the result. We suggest not use this function.
                    

example:
./Thickening2  -f ../data/mechanical04.obj -m=3 -i=2 -t=8 -d=0.6        
./Thickening2  -f ../data/tet.obj2 -m=3 -i=1 -t=8 
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

Please do not use this program for thickening with long offset distance. 
because in this case the speed and the result isn't better than the traditional method.
It can only use in short offset distance.

In the case of some extreme situations, the result will occur some very sharp parts. 
The argument -l can be appropriately reduced to alleviate these situations.

When the mesh complexity is high and the operation is slow, try adding -s=1 to argument.
This action may cause holes in the grid.

It is necessary to ensure that the mesh is closed, and the normal vectors of all faces are outward and correct.

If you find bug in the result, you can use result mode 1 and then use meshlab, tetwild or other tools to build topology.


## 📝 License
Materials in this repository are distributed under the following license:

> All software is licensed under the BSD 3-Clause License. See [LICENSE](https://github.com/rainbowwing/Thickening2/blob/main/LICENSE) file for details.