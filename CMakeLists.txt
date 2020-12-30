cmake_minimum_required(VERSION 3.13)
project(gfold VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
set(ENABLE_ECOS TRUE)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DownloadEpigraph.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DownloadMatplotlib.cmake)

add_executable(gfold main.cpp)
target_link_libraries(gfold
        PRIVATE Eigen3::Eigen
        PRIVATE epigraph
        PRIVATE mpl)
