#----------------------------------------------------------------------------
# Setup the project
#

## Set the default build type
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Set build type (Debug is default)")

set(DEV_WARNING "cmake build of libsbsdig is still in development. Please test and report issues")

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(sbsdig)

# Setup GNUInstallDirs for installation:
#----------------------------------------------------------------------------
include(GNUInstallDirs)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# The following lines are required so that the executable can dynamically locate the library dependencies independent of 
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#SET(CMAKE_BUILD_TYPE RelWithDebInfo)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#### So far this is like g4sbs. Good.

#----------------------------------------------------------------------------

## Disable in-source builds!
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "cmake cannot build libsbsdig in the root directory: "
    "${PROJECT_SOURCE_DIR}\nsince it conflicts with our old Makefile."
    " Create a build directory anywhere, change to it and re-run cmake."
    )
endif()

#### Yes, we want to keep that

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
#
file(GLOB headers ${PROJECT_SOURCE_DIR}/src/*.h)

##----------------------------------------------------------------------------
## Set the sources which have a corresponding .h file here
#set(sources ${PROJECT_SOURCE_DIR}/src/*.cxx ${PROJECT_SOURCE_DIR}/libsbsgem/*.cxx)
#set(sources src/gmn_tree.cxx src/SBSDigGEMPlane.cxx src/SBSDigAuxi.cxx 
#            src/SPEModel.cxx )

## Make corresponding header files
#string(REGEX REPLACE .cxx .h headers "${sources}")
## Append a few extra header files in the repository.  (But do not append the
## sbs_LinkDef.h header here, that is further down this cmake config.
#list(APPEND headers "src/g4sbs_types.h" "libsbsgem/types.h")
# list(APPEND headers "src/g4sbs_types.h")


string(REPLACE "src/" "" dictheaders1 "${headers}")
string(REPLACE "libsbsgem/" "" dictheaders "${dictheaders1}")

## Setup default variables (copied from Makefile)
#FIXME: should be configure'd:
set(EXTRAWARN ON CACHE BOOL "g++ Wextra warn options if supported")

##----------------------------------------------------------------------------
set(VERBOSE ON CACHE BOOL "Compile extra code for printing verbose messages")
set(TESTCODE ON CACHE BOOL "Compile extra diagnostic code (extra computations and global variables")
set(MCDATA ON CACHE BOOL "Compile support code for MC input data")
set(CXXMAXERRORS 0 CACHE INT "Maximum number of allowed errors before cxx stops")
list(APPEND SBSDIGEXTRADEF_LIST VERBOSE TESTCODE MCDATA)

#----------------------------------------------------------------------------
# Find ROOT 
#ROOT path
# commented out the following line so that this procedure will evolve with ROOT:
# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

add_executable(sbsdig sbsdig.cxx ${sources} ${headers})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
##----------------------------------------------------------------------------
## Setup the PODD Analyzer
#find_package(Podd 1.6 REQUIRED)

##----------------------------------------------------------------------------
## Setup the SBS library
#find_package(SBS REQUIRED)

##----------------------------------------------------------------------------
## Setup the G4SBS library
# find_package(G4SBS REQUIRED)


##----------------------------------------------------------------------------
## Check for SStream support
include(TestForSSTREAM)
if(CMAKE_HAS_ANSI_STRING_STREAM)
  set(HAS_SSTREAM ON)
  list(APPEND SBSDIGEXTRADEF_LIST HAS_SSTREAM)
endif()

##----------------------------------------------------------------------------
## OS Type (not sure what exactly uses this, but since it's in the Makefile
## I'll just copy it here too...
# set(SBSDIG_LINKFLAGS "")
if(UNIX) ## Can be anything from Linux to Mac OS X
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -O -g")
  if(APPLE)
    set(MACVERS ON)
    list(APPEND SBSDIGEXTRADEF_LIST MACVERS)
    set(SBSDIG_LINKFLAGS "${SBS_LINKFLAGS} -Wl,-undefined,dynamic_lookup")
  else()
    set(LINUXVERS ON)
    list(APPEND SBSDIGEXTRADEF_LIST LINUXVERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
    if(EXTRAWARN AND CMAKE_COMPILER_IS_GNUCXX AND NOT
        CMAKE_CXX_COMPILER_VERSION LESS 4)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-missing-field-initializers")
    endif()
    if(CXXMAXERRORS GREATER 0)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=${CXXMAXERRORS}")
    endif()
  endif()
endif() ## I guess no Windows support?



##----------------------------------------------------------------------------
## Copy over some flags we have set in the Makefile
if(VERBOSE)
  list(APPEND VERBOSE)
endif()

#include_directories(${ROOT_INCLUDE_DIR} ${Podd_INCLUDE_DIR} ${SBS_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/libsbsgem)
# include_directories(${ROOT_INCLUDE_DIR} ${G4SBS_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src)
include_directories(${ROOT_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

## Generate the ROOT dictionary
#set(PROJECT_DICTIONARY ${PROJECT_NAME}Dict)
#ROOT_GENERATE_DICTIONARY( ${PROJECT_DICTIONARY} ${dictheaders} OPTIONS -I${PROJECT_SOURCE_DIR}/src -I${PROJECT_SOURCE_DIR}/libsbsgem MODULE ${PROJECT_NAME} LINKDEF src/${PROJECT_NAME}_LinkDef.h )

## Create the library
#add_library(${PROJECT_NAME} SHARED ${sources} ${headers} ${PROJECT_DICTIONARY}.cxx )
set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${SBSDIG_LINKFLAGS}")
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES})
#message(FATAL_ERROR "${PROJECT_NAME} PRIVATE ${SBSEXTRADEFINES}")
## Get a list of preprocessor flags to pass to the targets
foreach(_def ${SBSEXTRADEF_LIST})
  if(${${_def}})
    target_compile_definitions(${PROJECT_NAME} PRIVATE ${_def})
  endif()
endforeach()

configure_file( ${PROJECT_SOURCE_DIR}/sbsdigenv.csh ${CMAKE_CURRENT_BINARY_DIR}/sbsdigenv.csh @ONLY)
configure_file( ${PROJECT_SOURCE_DIR}/sbsdigenv.sh ${CMAKE_CURRENT_BINARY_DIR}/sbsdigenv.sh)

message(AUTHOR_WARNING ${DEV_WARNING})

#install the library under ${CMAKE_INSTALL_PREFIX}/lib
#INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS sbsdig RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
#install the relevant directories: DB, replay, scripts in the top-level installation folder:
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/db DESTINATION ${CMAKE_INSTALL_PREFIX})
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/example DESTINATION ${CMAKE_INSTALL_PREFIX})
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sbsdigenv.sh DESTINATION ${CMAKE_INSTALL_BINDIR})
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sbsdigenv.csh DESTINATION ${CMAKE_INSTALL_BINDIR})
#install all headers under cmake_install_prefix/include
INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# install the ROOT dictionary files: eventually we should learn how to do this in a smarter way. Put the .cxx file under "include"
# and put the _rdict.pcm and .rootmap files under "lib". This seems to work. Don't ask me why
#INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_DICTIONARY}.cxx DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
