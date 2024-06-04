cmake_minimum_required(VERSION 3.21)
project(fast_whole_det_sim)

set(CMAKE_CXX_STANDARD 17)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp)


#find ROOT
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree)
#find YAML
find_package(yaml-cpp REQUIRED)

if(ROOT_FOUND)
    include(${ROOT_USE_FILE})
else()
    # Specify the ROOT headers directory
    include_directories(/home/ineffablord/documents/ROOT/root-6.28/include)
endif()

add_executable(gen_k40 ${PROJECT_BINARY_DIR}/gen_k40.cc ${sources} ${headers})

if(ROOT_FOUND)
    target_link_libraries(gen_k40 ${ROOT_LIBRARIES})
else()
    # Specify the ROOT libraries directory
    link_directories(/home/ineffablord/documents/ROOT/root-6.28/lib)
endif()
target_link_libraries(gen_k40 ${YAML_CPP_LIBRARIES})
