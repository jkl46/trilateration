cmake_minimum_required(VERSION 3.5)
project(main)

include_directories(include/)

set(SOURCE_FILES
    src/main.cpp
    src/trilaterate.cpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})