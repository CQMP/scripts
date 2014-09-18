cmake_minimum_required(VERSION 2.8)

project(backcont)

include_directories("/opt/local/include/")
link_directories("/opt/local/lib/")
add_executable(backcont backcont.cpp)

set(CMAKE_CXX_FLAGS  "-W -Wno-sign-compare -O3 -g")
target_link_libraries(backcont boost_program_options-mt)
