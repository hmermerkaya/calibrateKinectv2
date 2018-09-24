cmake_minimum_required(VERSION 2.8.9)
project(calibrateKinect CXX)

macro(REMOVE_VTK_DEFINITIONS)
    get_directory_property(_dir_defs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMPILE_DEFINITIONS)
    set(_vtk_definitions)
    foreach(_item ${_dir_defs})
    	message("item ${_item}")
        if(_item MATCHES "vtk*")
            list(APPEND _vtk_definitions -D${_item})
        endif()
    endforeach()
    remove_definitions(${_vtk_definitions})
endmacro(REMOVE_VTK_DEFINITIONS)

set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -std=c++11") 
set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_BUILD_TYPE  Debug)
include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#find_package(ROOT COMPONENTS RooStats MathCore RIO Hist Tree Net)
find_package(ROOT)
add_definitions(${ROOT_CXX_FLAGS})
#message ("root includes: ${ROOT_INCLUDE_DIRS}")
message("root libs:  ${ROOT_LIBRARIES}}")
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})
list(APPEND ROOT_LIBRARIES "/home/hamit/Downloads/root/lib/libGenVector.so")

option(WITH_CUDA "adds cuda support for freenect2" ON)
if(${WITH_CUDA})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_CUDA")
endif()

option(WITH_PCL "adds pcl cloud support" ON)
if(${WITH_PCL})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_PCL")
	find_package(PCL REQUIRED)
	link_directories(${PCL_LIBRARY_DIRS})
	add_definitions(${PCL_DEFINITIONS})
	include_directories(${PCL_INCLUDE_DIRS})
endif()



find_package(OpenCV     REQUIRED)

find_package(freenect2  REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS} ${freenect2_INCLUDE_DIRS})
REMOVE_VTK_DEFINITIONS()

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#ROOT_GENERATE_DICTIONARY(G__Event Event.h LINKDEF EventLinkDef.h)

#---Create a shared library with geneated dictionary
#add_library(Event SHARED Event.cxx G__Event.cxx)
#target_link_libraries(Event ${ROOT_LIBRARIES})

#---Create  a main program using the library


        
add_executable(calibrateKinect src/calibrateKinect.cpp)
target_link_libraries(calibrateKinect ${ROOT_LIBRARIES} ${OpenCV_LIBS} ${freenect2_LIBRARY}  ${PCL_LIBRARIES} )

set_target_properties( calibrateKinect PROPERTIES
        DEBUG_POSTFIX "d"
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

install(TARGETS calibrateKinect
        RUNTIME DESTINATION bin COMPONENT main
        LIBRARY DESTINATION lib COMPONENT main
        ARCHIVE DESTINATION lib COMPONENT main)
