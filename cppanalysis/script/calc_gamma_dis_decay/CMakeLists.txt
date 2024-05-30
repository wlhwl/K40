cmake_minimum_required(VERSION 3.21)
project(analysis_k40)

set(CMAKE_CXX_STANDARD 17)

# Specify the ROOT headers directory
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree)
include(${ROOT_USE_FILE})

add_executable(analysis analysis_trigger_root_tts.cpp analysis_trigger_root_tts.h)
target_link_libraries(analysis ${ROOT_LIBRARIES})
