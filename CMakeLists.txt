cmake_minimum_required(VERSION 3.2)

set(PROJECT_NAME_STR builder_sample)
set(CMAKE_CXX_STANDARD 11)

project (${PROJECT_NAME_STR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

add_executable(${PROJECT_NAME_STR} ${SRC_FILES})
