cmake_minimum_required(VERSION 3.10)
project(fit VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#link_directories(${CMAKE_BINARY_DIR})
#
#target_link_libraries(GLBall mylib)
set(SOURCE_FILES fit/main.cpp)

add_executable(fit.bin ${SOURCE_FILES})
#add_executable(fit.bin fit/main.cpp)

#add_subdirectory(fit)
#add_library( networkit SHARED IMPORTED )
#set_target_properties( networkit PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/libnetworkit.so )
#target_link_libraries(networkit networkit)



#find_package(networkit REQUIRED)
add_library(fit SHARED
    #fit/SortedCounter.cpp
   fit/G_graph.cpp
   fit/H_graph.cpp
   fit/PLM.cpp
   fit/history_graph.cpp
   fit/metrics.cpp
)
target_include_directories(fit PUBLIC
    include
    ) # add unit test directory ?

link_directories(${CMAKE_BINARY_DIR})

target_link_libraries(fit libnetworkit.so)
target_link_libraries (fit.bin PUBLIC
    fit
)
#IF(TBB_FOUND AND SQLite3_FOUND)
#	LINK_LIBRARIES(
#		TBB::tbb
#		SQLite::SQLite3
#	)
#ENDIF()

# link catch2 for unit tests
find_package(Catch2 3 REQUIRED)
# These tests can use the Catch2-provided main

add_executable(tests tests/test.cpp)
target_link_libraries (tests PUBLIC
    fit
)

target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

# These tests need their own main
#add_executable(custom-main-tests test.cpp )
#target_link_libraries(tests PRIVATE Catch2::Catch2)


# In case you put datasets elsewhere, e.g., using WSL
#IF(NOT DEFINED ENV{DATASET_DIR})
#	SET(ENV{DATASET_DIR} "${CMAKE_SOURCE_DIR}/data")
#ENDIF()
#FILE(TO_CMAKE_PATH $ENV{DATASET_DIR} DATASET_DIR)
#ADD_COMPILE_DEFINITIONS(DATASET_DIR="${DATASET_DIR}/")

# Commonly used optimization flags
add_compile_options(-O3)
target_link_options(${PROJECT_NAME} PUBLIC
    "LINKER:-no-as-needed"
)

