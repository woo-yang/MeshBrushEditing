cmake_minimum_required(VERSION 3.0)

project(MeshBrushEditing)
set(CMAKE_CXX_STANDARD 17)

if(MSVC)
    add_definitions(-D_USE_MATH_DEFINES)
endif()

FILE(GLOB SOURCES src/*.cpp)

INCLUDE_DIRECTORIES(
    YRender/inc
    $ENV{INCLUDE_DIR}
    $ENV{INCLUDE_DIR}/Eigen
    $ENV{INCLUDE_DIR}/gmp
)


link_directories(
    $ENV{LIB_DIR}
    $ENV{LIB_DIR}/gmp
    $ENV{LIB_DIR}/boost

)

add_subdirectory(common)
add_subdirectory(alg)
add_subdirectory(solver)
add_subdirectory(interface)

add_executable(${PROJECT_NAME} ${SOURCES})


target_link_libraries(${PROJECT_NAME} PUBLIC
    common
    algorithm
    interface
    )
