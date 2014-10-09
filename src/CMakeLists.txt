cmake_minimum_required (VERSION 2.8) 
project (CeresBA) 
set (CeresBA_VERSION_MAJOR 0)
set (CeresBA_VERSION_MINOR 1)

set(CMAKE_BUILD_TYPE "DEBUG")

MESSAGE( STATUS "Build ${CMAKE_PROJECT_NAME} in mode: " ${CMAKE_BUILD_TYPE} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x")


FIND_PACKAGE(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})


# NORMAL Executable
add_executable (BAdjustBin ceresBAbin.cpp)
TARGET_LINK_LIBRARIES(BAdjustBin ${CERES_LIBRARIES})


# MEX Executable
# NOTE: Adapt path to matlab binary in mex_link.sh  as well!
set(MATLAB_INCLUDE_DIR "/usr/local/MATLAB/MATLAB_Production_Server/R2013a/extern/include/")
#set(MATLAB_INCLUDE_DIR "/usr/pack/matlab-8.3r2014a-fg/extern/include/")

# NOTE: The following path has to be set manually if ceres is not installed in a standard location. Adapt to your local installation!
# If Ceres is installed in a location visible to Matlab, set to "".
#set(CERESLIBSPATH "-L/home/kroegert/local/usrlocal/lib")
set(CERESLIBSPATH "")


include_directories(${MATLAB_INCLUDE_DIR})
set(MEXCERESLIBS "-lglog -lceres")



add_library (BAdjustMex ceresBAmex.cpp)
SET_TARGET_PROPERTIES(BAdjustMex PROPERTIES COMPILE_FLAGS -fPIC)
ADD_CUSTOM_COMMAND(TARGET BAdjustMex
    POST_BUILD
    COMMAND ${CMAKE_SOURCE_DIR}/mex_link.sh ${CMAKE_SOURCE_DIR} BAdjustMex ${CERESLIBSPATH} "${MEXCERESLIBS}"
    )


enable_testing()
add_test( BAtest1 BAdjustBin  ../src/BAproblem_test.txt  /dev/null)



