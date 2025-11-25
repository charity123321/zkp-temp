cmake_minimum_required(VERSION 3.10)
project(zk-tools)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED gmp)

cmake_minimum_required(VERSION 3.16)
project(zk-tools)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加可执行文件 - 只包含源文件(.cpp)
add_executable(zk-tools
    main.cpp
    
    # 算术模块
    arithmetic/virtual_polynomial.cpp
    arithmetic/util.cpp
    
    # 子程序模块

    
    # 其他模块
    
)

# 包含头文件目录 - 现在需要包含所有头文件所在目录
target_include_directories(zk-tools PRIVATE 
    .                           # 当前目录
    include                     # 如果有统一的include目录
    arithmetic                  # 算术头文件
    subroutines/poly_iop        # 子程序头文件
    transcript                  # transcript头文件
    field_and_polynomial        # 域和多项式头文件
)

# 可选：设置编译器标志
target_compile_options(zk-tools PRIVATE 
    -Wall
    -Wextra
    -Werror
)

target_include_directories(zk-tools PRIVATE ${GMP_INCLUDE_DIRS} .)
target_link_libraries(zk-tools PRIVATE ff ${GMP_LIBRARIES})