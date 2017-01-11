#-------------------------------------------------------------------------------
# Checking CMAKE version
cmake_minimum_required( VERSION 2.8 FATAL_ERROR )


set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}-O3 -msse4.2 -m64 -fpermissive -Wno-sign-compare -Wextra -Wall -pedantic -Wno-unused-result" CACHE STRING "" )
set( CMAKE_INSTALL_PREFIX /usr/local/bin )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_definitions(-DHAVE_CONFIG_H)

#-------------------------------------------------------------------------------
# Defining the project
project( rtx2mnc )
set( CMAKE_BUILD_TYPE RELEASE )
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules
	${CMAKE_MODULE_PATH} )

#-------------------------------------------------------------------------------
# Find the package ...
find_package( MINC REQUIRED )

#-------------------------------------------------------------------------------
# Setup include directory for this project
include_directories( ${MINC_INCLUDE_DIR} )

find_package (Threads) 
find_package(DCMTK REQUIRED)
include_directories(${DCMTK_INCLUDE_DIRS})

#-------------------------------------------------------------------------------
# Add the executable, and link to libraries
add_executable( rtx2mnc main.cpp )
target_link_libraries( rtx2mnc ${MINC_minc2_LIBRARY} ${MINC_volume_io2_LIBRARY} ${DCMTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )

#-------------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install( TARGETS rtx2mnc DESTINATION bin )
