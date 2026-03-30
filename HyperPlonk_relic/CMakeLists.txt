cmake_minimum_required(VERSION 3.16)
project(HyperPlonk_relic)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ==========================================
# 0. 寻找 OpenMP 多线程库 (新增)
# ==========================================
find_package(OpenMP REQUIRED)

# ==========================================
# 1. 寻找外部依赖库 (AVX-MSM)
# ==========================================

#set(AVXMSM_LOCAL_ROOT "/home/lgq/sys/SimdMSM/AVX-MSM") 

set(MSM_LIBRARY "/home/lgq/sys/zkp-temp/HyperPlonk_relic/SimdMSM/AVX-MSM/demo/381/build/libmsm.a")
set(RELIC_LIBRARY "/home/lgq/sys/zkp-temp/HyperPlonk_relic/SimdMSM/AVX-MSM/demo/381/target/lib/librelic_s.a")

if(NOT MSM_LIBRARY OR NOT RELIC_LIBRARY)
    message(FATAL_ERROR "❌ 找不到静态库！MSM状态: ${MSM_LIBRARY}, RELIC状态: ${RELIC_LIBRARY}")
else()
    message(STATUS "✅ 找到 MSM 库: ${MSM_LIBRARY}")
    message(STATUS "✅ 找到 RELIC 库: ${RELIC_LIBRARY}")
endif()

# ==========================================
# 2. 定义可执行文件
# ==========================================
add_executable(HyperPlonk_relic
    main.cpp
    arithmetic/util.cpp
    # 其他模块源文件...
)

# ==========================================
# 3. 设置头文件包含目录
# ==========================================
target_include_directories(HyperPlonk_relic PRIVATE 
    .                       # 当前目录
    ${GMP_INCLUDE_DIRS}     # GMP 头文件
    ${AVXMSM_INCLUDE_DIR}   # AVX-MSM 头文件
    #arithmetic                  # 算术头文件
    #subroutines/poly_iop        # 子程序头文件
    #subroutines/pcs
    #hyperplonk 
    #transcript                  # transcript头文件
    #field_and_polynomial        # 域和多项式头文件
    # 包含 MSM 的头文件
    "SimdMSM/AVX-MSM/demo/381/include"
    # 包含 RELIC 的头文件
    "SimdMSM/AVX-MSM/include"
    "SimdMSM/AVX-MSM/demo/381/target/include"
    #"/home/lgq/sys/SimdMSM/AVX-MSM/include"
    #"/home/lgq/sys/SimdMSM/AVX-MSM/demo/381/target/include"
)

# ==========================================
# 4. 设置编译器标志
# ==========================================
target_compile_options(HyperPlonk_relic PRIVATE 
    -Wall
    -Wextra
    -Werror
    -fopenmp
    # --- 针对 MSM 优化的硬件指令集 ---
    -mavx512f -mavx512dq -mbmi2 -madx 
)

# ==========================================
# 5. 链接所有静态库/动态库
# ==========================================
target_link_libraries(HyperPlonk_relic PRIVATE 
    ff 
    ${MSM_LIBRARY}       # 先放上层库 (MSM)
    ${RELIC_LIBRARY}     # 再放底层库 (RELIC)
    gmp
    pthread              # MSM 通常依赖多线程
    OpenMP::OpenMP_CXX   # <--- 新增：让 CMake 自动处理 OpenMP 的编译和链接标志
)