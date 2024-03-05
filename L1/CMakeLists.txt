
cmake_minimum_required(VERSION 3.12)

project(linear_search)

set(CMAKE_CXX_STANDARD_REQUIRED 11)

add_executable(
    ${PROJECT_NAME}_2d_basic
    src/linear_search_2d.cpp
)

add_executable(
    ${PROJECT_NAME}_2d_Robiins
    src/linear_search_2d_with_Robiins_Monro_Rule.cpp
)
add_executable(
    ${PROJECT_NAME}_2d_Amijio
    src/linear_search_2d_with_amijio.cpp
)


