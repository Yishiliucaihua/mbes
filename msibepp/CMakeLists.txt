cmake_minimum_required (VERSION 2.8)

project (msibepp)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

include_directories (src)

aux_source_directory (src SRC_LIST)

add_executable(msibepp ${SRC_LIST})

target_link_libraries(msibepp m)
