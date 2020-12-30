cmake_minimum_required(VERSION 3.13)

include(FetchContent)
FetchContent_GetProperties(epigraph)

if (NOT epigraph)
    FetchContent_Declare(
            epigraph
            GIT_REPOSITORY https://github.com/EmbersArc/Epigraph.git
    )
    FetchContent_Populate(epigraph)

    add_subdirectory(${epigraph_SOURCE_DIR})
endif ()
