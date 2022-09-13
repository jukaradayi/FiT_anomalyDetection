cmake_minimum_required(VERSION 3.10)
project(fit VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#link_directories(${CMAKE_BINARY_DIR})
#
#target_link_libraries(GLBall mylib)


include_directories(
    include
    ) # add unit test directory ?

#add_library( networkit SHARED IMPORTED )
#set_target_properties( networkit PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/libnetworkit.so )
#target_link_libraries(networkit networkit)



#find_package(networkit REQUIRED)
add_library(fit STATIC
   fit/G_graph.cpp
   fit/H_graph.cpp
   fit/PLM.cpp
   fit/history_graph.cpp
   fit/metrics.cpp
)
link_directories(${CMAKE_BINARY_DIR})

target_link_libraries(fit libnetworkit)

#IF(TBB_FOUND AND SQLite3_FOUND)
#	LINK_LIBRARIES(
#		TBB::tbb
#		SQLite::SQLite3
#	)
#ENDIF()


# In case you put datasets elsewhere, e.g., using WSL
#IF(NOT DEFINED ENV{DATASET_DIR})
#	SET(ENV{DATASET_DIR} "${CMAKE_SOURCE_DIR}/data")
#ENDIF()
#FILE(TO_CMAKE_PATH $ENV{DATASET_DIR} DATASET_DIR)
#ADD_COMPILE_DEFINITIONS(DATASET_DIR="${DATASET_DIR}/")

# Commonly used optimization flags
add_compile_options(-O3)

#add_executable(fit.bin fit/main.cpp)
