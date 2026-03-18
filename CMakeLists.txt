cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME demo)
project(${PROJECT_NAME})

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

find_package(MSHR REQUIRED)

add_definitions(${MSHR_CXX_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MSHR_CXX_FLAGS}")
include_directories(${MSHR_INCLUDE_DIRS})
include_directories(SYSTEM ${MSHR_EXTERNAL_INCLUDE_DIRS})
link_directories(${MSHR_LIBRARIES_DIRS})

add_executable(poisson Poisson.cpp)
target_link_libraries(poisson dolfin mshr)