cmake_minimum_required(VERSION 3.13)

include(FetchContent)
FetchContent_GetProperties(mpl)

if (NOT mpl_POPULATED)
    find_package(Python3 COMPONENTS Development NumPy)
    FetchContent_Declare(
            mpl
            GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    )
    FetchContent_Populate(mpl)

    file(COPY "${mpl_SOURCE_DIR}/matplotlibcpp.h"
            DESTINATION "${mpl_BINARY_DIR}/include/mpl")

    add_library(mpl INTERFACE IMPORTED)
    target_include_directories(mpl
            INTERFACE "${mpl_BINARY_DIR}/include/"
            INTERFACE "${PYTHON_INCLUDE_DIRS}"
            INTERFACE "${Python3_NumPy_INCLUDE_DIRS}")
    target_link_libraries(mpl
            INTERFACE Python3::Python
            INTERFACE Python3::NumPy)
endif ()
