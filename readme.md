# PFP-Thickening

This project is a demo to achieve to algorithm of PFP-Thickening.

[//]: # (|  General  | [![c++14]&#40;https://img.shields.io/badge/standard-C++14-blue.svg?style=flat&logo=c%2B%2B&#41;]&#40;https://isocpp.org&#41; [![License]&#40;https://img.shields.io/badge/License-BSD_3--Clause-orange.svg&#41;]&#40;https://github.com/robotology/osqp-eigen/blob/master/LICENSE&#41; |)

[//]: # (| :-------: | :----------------------------------------------------------: |)

[//]: # (| **CI/CD** | [![Codacy Badge]&#40;https://app.codacy.com/project/badge/Grade/a18710c10f1c4df19bc2759fd50e9cf5&#41;]&#40;https://www.codacy.com/gh/robotology/osqp-eigen/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=robotology/osqp-eigen&amp;utm_campaign=Badge_Grade&#41; [![CI]&#40;https://github.com/robotology/osqp-eigen/workflows/C++%20CI%20Workflow/badge.svg&#41;]&#40;https://github.com/robotology/osqp-eigen/workflows/C++%20CI%20Workflow/badge.svg&#41; [![Azure]&#40;https://dev.azure.com/conda-forge/feedstock-builds/_apis/build/status/osqp-eigen-feedstock?branchName=master&#41;]&#40;https://dev.azure.com/conda-forge/feedstock-builds/_build/results?buildId=341091&view=results&#41; |)

[//]: # (| **conda** | [![Conda Recipe]&#40;https://img.shields.io/badge/recipe-osqp--eigen-green.svg&#41;]&#40;https://anaconda.org/conda-forge/osqp-eigen&#41;  [![Conda Downloads]&#40;https://img.shields.io/conda/dn/conda-forge/osqp-eigen.svg&#41;]&#40;https://anaconda.org/conda-forge/osqp-eigen&#41;  [![Conda Version]&#40;https://img.shields.io/conda/vn/conda-forge/osqp-eigen.svg&#41;]&#40;https://anaconda.org/conda-forge/osqp-eigen&#41;  [![Conda Platforms]&#40;https://img.shields.io/conda/pn/conda-forge/osqp-eigen.svg&#41;]&#40;https://anaconda.org/conda-forge/osqp-eigen&#41; |)

[![c++20](https://img.shields.io/badge/standard-C++20-blue.svg?style=flat&logo=c%2B%2B)](https://isocpp.org)

![BuildExamplesLinux](https://github.com/rainbowwing/Thickening2/workflows/CMake/badge.svg)

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
```bash
sudo su root
apt update
apt upgrade
brew install cgal
brew install eigen
brew install gflags
g++ -v
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

USAGE: Thicken2 -f file_name [options]    file_name must *.obj or *.obj2      

options:
  -m={1|2|3}                               means result mode which value can be chose in 1,2 and 3. 
                                           mode 1 get the result without mesh and remeshing;
                                           mode 2 get the result with building mesh;
                                           mode 3 get the result with building mesh and remeshing.
  -d=<num>                                 this number is a double means the expect length of each facet in running invariable thickening.
                                           Limited in 0.1~1.5. This number does not represent an absolute distance.
                                           Example -d=0.5, indicating that the offset distance is 0.5 times the average mesh edge length.
  -l=<num>                                 This number is a double.
                                           which value indicates how many times the maximum offset distance is the ideal offset distance limited in 1.5~2.7.
                                           You can set it is 2.0 .
  -r <value>                                   value is (0|1|2) 
  -d <value>                                   value is (0|1|2)
  -l <value>                                   value is (0|1|2) 
  -g <value>                                   value is (0|1|2) 
  -s <value>                                   value is (0|1|2) 
  -t <value>                                   value is (0|1|2) 
```

DEFINE_int32(m, 3,
"This arg is a integer means result mode which value can be chose in 1,2,3. \n mode 1  without step 6;mode 2 is with building mesh;mode 3 is with building mesh and remeshing.");
DEFINE_double(d, 0.5,
"This arg is a double means the length running invariable thickening limited in 0.1~1.5. Example -d 0.5, Indicating that the offset distance is 0.75 times the average mesh edge length.");
DEFINE_double(l, 2.0,
"This arg is a double which value indicates how many times the maximum offset distance is the ideal offset distance limited in 1.5~2.7. .You can set it is 2.0");
DEFINE_double(g, -1,
"This arg is a double means the length of edge length of each cell in grid length. If you can't calculate a length with better performance, it can be passed. Then it will use the default value.");
DEFINE_int32(t, 12, "thread num please set this value depend the cpu of you device.");
DEFINE_string(f, "", "file name, which can choose *.obj2 or *.obj.");
DEFINE_int32(i, 1,
"This arg is a integer means running mode which value can be chose in 1,2. \n 1 is offsetting to the outside of the mesh; 2 is offsetting to the inside of the mesh.");
DEFINE_double(e, 1e-5,
"This arg is a double means the eps. When the distance of two points is smaller than eps, we will regard these two point as coinciding. ");
DEFINE_bool(s, false,
"This arg is a bool value which means the program will skip some cell which is most likely useless. It can improve performance, but may cause holes in the result. We suggest not use this function.");


### The input file *.obj2 format

which is same as *.obj.

v x y z

### NOTE:

cgal version must be larger than 5.5

在一些极限网格的情况中会出现非常尖锐的偏移，可以适当的减小l，来缓解该情况

当网格复杂都高，运行慢的时候可尝试性在arg中加入 -s 1。此举可能造成网格破洞。

args:

-m 0带remeshing 1带重建 2 都不带

-r 带remeshing

-d 1.5 must in (0.1~1.5) 表示进行等距偏移

-l 最大长度(1~2.5)过长在极限情况造成边角

-g grid len 修改大小可能可以提高性能，如果不确定多大可以使用默认的

-s skipmode

-t 线程数目

obj2 格式：

需要保证网格封闭，并且所有面的法向量向外且正确，不封闭的话，可以使用工具进行修复。本demo暂未实现。

如果遇到结果的bug，可以使用mode 1 然后再借助meshlab，tetwild等工具建立拓扑。


