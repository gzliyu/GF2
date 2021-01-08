cmake_minimum_required(VERSION 3.2.2 FATAL_ERROR)
###############################################################################################################
# This script will build and run the examples from a folder
# Execute from a command line:
#     ctest -S HDF4_Examples.cmake,OPTION=VALUE -C Release -VV -O test.log
###############################################################################################################

set(CTEST_CMAKE_GENERATOR "Visual Studio 14 2015 Win64")
if("")
  set(CMAKE_GENERATOR_TOOLSET "")
endif()
set(CTEST_DASHBOARD_ROOT ${CTEST_SCRIPT_DIRECTORY})

# handle input parameters to script.
#INSTALLDIR - HDF4 root folder
#CTEST_CONFIGURATION_TYPE - Release, Debug, RelWithDebInfo
#CTEST_SOURCE_NAME - name of source folder; HDF4Examples
if(DEFINED CTEST_SCRIPT_ARG)
  # transform ctest script arguments of the form
  # script.ctest,var1=value1,var2=value2
  # to variables with the respective names set to the respective values
  string(REPLACE "," ";" script_args "${CTEST_SCRIPT_ARG}")
  foreach(current_var ${script_args})
    if("${current_var}" MATCHES "^([^=]+)=(.+)$")
      set("${CMAKE_MATCH_1}" "${CMAKE_MATCH_2}")
    endif()
  endforeach()
endif()

###################################################################
### Following Line is one of [Release, RelWithDebInfo, Debug] #####
set(CTEST_CONFIGURATION_TYPE "$ENV{CMAKE_CONFIG_TYPE}")
if(NOT DEFINED CTEST_CONFIGURATION_TYPE)
  set(CTEST_CONFIGURATION_TYPE "Release")
endif()
set(ADD_BUILD_OPTIONS "${ADD_BUILD_OPTIONS} -DCTEST_CONFIGURATION_TYPE:STRING=${CTEST_CONFIGURATION_TYPE}")
##################################################################

if(NOT DEFINED INSTALLDIR)
  set(INSTALLDIR "C:/Program")
endif()

if(NOT DEFINED CTEST_SOURCE_NAME)
  set(CTEST_SOURCE_NAME "HDF4Examples")
endif()

if(NOT DEFINED HDF_LOCAL)
  set(CDASH_LOCAL "NO")
else()
  set(CDASH_LOCAL "YES")
endif()
if(NOT DEFINED CTEST_SITE)
  set(CTEST_SITE "local")
endif()
if(NOT DEFINED CTEST_BUILD_NAME)
  set(CTEST_BUILD_NAME "examples")
endif()
set(ADD_BUILD_OPTIONS "${ADD_BUILD_OPTIONS} -DSITE:STRING=${CTEST_SITE} -DBUILDNAME:STRING=${CTEST_BUILD_NAME}")

#TAR_SOURCE - name of tarfile
#if(NOT DEFINED TAR_SOURCE)
#  set(CTEST_USE_TAR_SOURCE "HDF4Examples-0.2.4-Source")
#endif()

###############################################################################################################
if(WIN32)
  set(SITE_OS_NAME "Windows")
  set(ENV{HDF4_DIR} "${INSTALLDIR}/cmake")
  set(CTEST_BINARY_NAME ${CTEST_SOURCE_NAME}\\build)
  set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}\\${CTEST_SOURCE_NAME}")
  set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}\\${CTEST_BINARY_NAME}")
else()
  set(ENV{HDF4_DIR} "${INSTALLDIR}/share/cmake")
  set(ENV{LD_LIBRARY_PATH} "${INSTALLDIR}/lib")
  set(CTEST_BINARY_NAME ${CTEST_SOURCE_NAME}/build)
  set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
  set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")
endif()
if(${CDASH_LOCAL})
  set(ADD_BUILD_OPTIONS "${ADD_BUILD_OPTIONS} -DCDASH_LOCAL:BOOL=ON")
endif()
set(ADD_BUILD_OPTIONS "${ADD_BUILD_OPTIONS} -DHDF4_PACKAGE_NAME:STRING=hdf4")

###############################################################################################################
# For any comments please contact cdashhelp@hdfgroup.org
#
###############################################################################################################

if(WIN32)
  include(${CTEST_SCRIPT_DIRECTORY}\\HDF4_Examples_options.cmake)
  include(${CTEST_SCRIPT_DIRECTORY}\\CTestScript.cmake)
else()
  include(${CTEST_SCRIPT_DIRECTORY}/HDF4_Examples_options.cmake)
  include(${CTEST_SCRIPT_DIRECTORY}/CTestScript.cmake)
endif()
