# HyperPlonk零知识证明实现

HyperPlonk通过将Plonk移植到布尔超立方体上，从而避免了对FFT的依赖。本项目旨在提高HyperPlonk协议的运行效率。目前项目主要针对多项式承诺KZG中的**多标量乘法**进行优化，后续会对协议本身依赖组件进行优化等工作。（注意项目中transcript目前还是简易实现版本，需要进一步的实现）

## 环境配置

本项目基于libff库使用c++实现HyperPlonk，即HyperPlonk_libff。后基于relic进行了更优的多标量乘法实现，由于研究成果未发表原因暂不能开源，只将[SimdMSM](https://github.com/JR-account/SimdMSM)方案整合到了HyperPlonk中，即HyperPlonk_relic。本节主要介绍如何配置运行该项目。

**获取源码**

```bash
git clone git@github.com:See-You-Soon-111/zkp-temp.git
```

**安装libff**

```bash
git clone https://github.com/scipr-lab/libff.git
cd libff
git submodule init && git submodule update

sudo apt-get install build-essential git libboost-all-dev cmake libgmp3-dev libssl-dev libprocps3-dev pkg-config libsodium-dev

make build && cd build
cmake ..

make
make install
```

**运行HyperPlonk_libff**

```bash
cd zkp-temp/HyperPlonk_libff

mkdir build && cd build
cmake ..
make
```

**安装relic库和SimdMSM库**

```bash
cd zkp-temp/HyperPlonk_relic/SimdMSM/AVX-MSM/demo/381
make lib

# 出现权限错误则设置权限
chmod +x ../../preset/x64-pbc-bls12-381.sh
sed -i 's/\r$//' ../../preset/x64-pbc-bls12-381.sh

# 生成静态链接库libmsm.a
make msm
```

在`HyperPlonk_relic/CMakeLists.txt`中添加`libmsm.a`和`librelic_s.a`的绝对路径

最后在`zkp/HyperPlonk_relic`目录下

```bash
mkdir build && cd build 
cmake ..
make
```



在对应的`test.h`文件中可设置参数（为了避免取整原因引发的报错，证据数量默认为2的幂次）。注：目前只支持代码模拟生成的电路，还无法直接连通业界主流的电路前端工具。





