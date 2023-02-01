cmake_minimum_required(VERSION 3.18.0)

project(smg CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wreturn-type -Wempty-body -Wno-deprecated -Wl,-no-as-needed")

#ADD_DEFINITIONS(-DNDEBUG)

include_directories(include)

add_executable(smg main.cpp)
