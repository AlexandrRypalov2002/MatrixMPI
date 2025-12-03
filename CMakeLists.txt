cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

PROJECT(matrix-mpi
    LANGUAGES CXX
)

set(PROJECT_NAME_SERIAL ${PROJECT_NAME}-serial)

find_package(MPI REQUIRED)

add_executable(${PROJECT_NAME})
add_executable(${PROJECT_NAME_SERIAL})

add_subdirectory(math)

target_sources(${PROJECT_NAME}
    PRIVATE
        main.cpp
)

target_sources(${PROJECT_NAME_SERIAL}
    PRIVATE
        main_serial.cpp
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        MPI::MPI_CXX

        MathLib
)

target_link_libraries(${PROJECT_NAME_SERIAL}
    PRIVATE
        MathLib
)