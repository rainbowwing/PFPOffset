## install bash
```bash
sudo su root
apt update
apt upgrade
apt install libcgal-dev # brew install cgal
apt install libeigen3-dev # brew install eigen
apt install libgflags-dev # brew install gflags
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
mkdir build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make


```

NOTE:
cgal version must be larger than 5.5

args:

-m 0带remeshing 1带重建 2 都不带 

-r 带remeshing 

-d 1.5 must in (0.1~1.5) 表示进行等距偏移

-l 最大长度(1~2.5)过长在极限情况造成边角

-g grid len 修改大小可能可以提高性能，如果不确定多大可以使用默认的

-s 当超过多少时使用小格子 

-t 线程数目


obj2 格式：

需要保证网格封闭，并且所有面的法向量向外且正确，不封闭的话，可以使用工具进行修复。本demo暂未实现。

如果遇到结果的bug，可以使用mode 1 然后再借助meshlab，tetwild等工具建立拓扑。


