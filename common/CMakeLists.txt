file(GLOB header *.h *.hpp)
file(GLOB source *.cpp)

if(MSVC)
    add_definitions(-D_USE_MATH_DEFINES)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    set(BUILD_SHARED_LIBS TRUE)
endif()

add_library(common SHARED ${source} ${header})
SET_TARGET_PROPERTIES(common PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../)
SET_TARGET_PROPERTIES(common PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(common 
    libgmp-10
    libmpfr-4
    libboost_*
)
