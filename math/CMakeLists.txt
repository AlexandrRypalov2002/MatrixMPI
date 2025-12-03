set(MATH_LIB_NAME MathLib)

add_library(${MATH_LIB_NAME} STATIC)

target_include_directories(${MATH_LIB_NAME}
    PUBLIC
        include
)

target_sources(${MATH_LIB_NAME}
    PUBLIC
        ${CMAKE_SOURCE_DIR}/math/include/matrix/Matrix.h
        ${CMAKE_SOURCE_DIR}/math/include/matrix/factory/MatrixFactory.h

    PRIVATE
        ${CMAKE_SOURCE_DIR}/math/src/matrix/Matrix.cpp
        ${CMAKE_SOURCE_DIR}/math/src/matrix/factory/MatrixFactory.cpp
)
