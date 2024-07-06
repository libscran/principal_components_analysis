cmake_minimum_required(VERSION 3.14)

project(scran_principal_component_analysis
    VERSION 1.0.0
    DESCRIPTION "Principal components analysis for single-cell data"
    LANGUAGES CXX)

add_library(scran_principal_component_analysis INTERFACE)
target_compile_features(scran_principal_component_analysis INTERFACE cxx_std_17)
target_include_directories(scran_principal_component_analysis INTERFACE include/scran)

add_subdirectory(extern)
target_link_libraries(scran_principal_component_analysis INTERFACE tatami::tatami tatami::tatami_stats ltla::irlba Eigen3::Eigen)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()
