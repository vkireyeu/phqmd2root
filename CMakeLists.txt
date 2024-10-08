cmake_minimum_required(VERSION 2.8.12)

project(phqmd_event)

find_package(ROOT 6.28 REQUIRED )
find_package(ZLIB)


include_directories(${ROOT_INCLUDE_DIRS})
include(${ROOT_USE_FILE})

set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

add_subdirectory(src)
