cmake_minimum_required(VERSION 3.16)

project(BranaWeylTOV LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(brana_tov
    src/main.cpp
    src/eos.cpp
    src/mcmc.cpp
    src/pchip.cpp
    src/tov.cpp
)

target_include_directories(brana_tov PRIVATE includes)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(brana_tov PRIVATE OpenMP::OpenMP_CXX)
endif()

if (MSVC)
    target_compile_options(brana_tov PRIVATE /W4)
else()
    target_compile_options(brana_tov PRIVATE -Wall -Wextra -pedantic -O2)
endif()
