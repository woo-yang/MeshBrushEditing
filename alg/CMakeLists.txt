file(GLOB header *.h)
file(GLOB source *.cpp)

if(MSVC)
    add_definitions(-D_USE_MATH_DEFINES)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    set(BUILD_SHARED_LIBS TRUE)
endif()



add_library(algorithm SHARED ${header} ${source})
SET_TARGET_PROPERTIES(algorithm PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../)

target_link_libraries(algorithm
    common
    libgmp-10
    libmpfr-4
    libboost_*

)
