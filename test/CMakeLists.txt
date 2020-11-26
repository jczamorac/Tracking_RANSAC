#-----------------------------------------------------------------------
# - Top Level CMakeLists.txt for mytest Build
#-----------------------------------------------------------------------
#
#-----------------------------------------------------------------------
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

#-----------------------------------------------------------------------
# - Project definition and basic configuration
#   Version handled manually as project(... VERSION ...) is not used
#   in tests/examples which are subprojects. All calls must use
#   the same form.
project(mytest)
set(${PROJECT_NAME}_VERSION_MAJOR  0)
set(${PROJECT_NAME}_VERSION_MINOR  1)
set(${PROJECT_NAME}_VERSION_PATCH  0)
set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

# - Prepend our own CMake Modules to the search path
#
set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/Modules
    ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH
    "$ENV{ROOTSYS}/etc/cmake"
    ${CMAKE_MODULE_PATH})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_definitions( -DTIXML_USE_STL )

include_directories(${CMAKE_SOURCE_DIR}/include)
link_directories(${CMAKE_SOURCE_DIR}/lib)


# - Find ROOT Modules
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -std=c++11 -pthread -Wno-deprecated-declarations")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lncurses" )
include_directories(${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

## Override root's ROOT_GENERATE_DICTIONARY macro to be less error-prone. This
## is done by only allowing for explicitly speficied header locations/include
## directories, instead of using the full search path for this module.
function(ROOT_GENERATE_DICTIONARY dictionary)
  CMAKE_PARSE_ARGUMENTS(ARG "" "" "LINKDEF;OPTIONS" "" ${ARGN})
  ## Get all include directories
  get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
  set(extra_includes "")
  foreach(dir ${dirs})
    set(extra_includes ${extra_includes} "-I${dir}")
  endforeach()
  ## find and call ROOTCLING
  find_program(ROOTCLING rootcling)
  add_custom_command(
    OUTPUT ${dictionary}.cxx ${dictionary}_rdict.pcm
    COMMAND ${ROOTCLING} -f ${dictionary}.cxx -s ${dictionary}.pcm
    ${ARG_OPTIONS} ${ARG_UNPARSED_ARGUMENTS} ${extra_includes} ${ARG_LINKDEF}
    DEPENDS ${includes} ${linkdefs})
endfunction()


# - Configure output Modules
#
#set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

file(GLOB HDRFILES ${PROJECT_SOURCE_DIR}/include/*.h)


add_executable(find_tracks main.cc ${HDRFILES} )
target_link_libraries(find_tracks ${ROOT_LIBRARIES} ${HDRFILES} Trackinglib)
