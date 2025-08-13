cmake_minimum_required(VERSION 3.15)
project(test_boost LANGUAGES CXX)

# Find Boost with program_options component
find_package(Boost REQUIRED COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(main main.cpp)
    target_link_libraries(main PRIVATE ${Boost_LIBRARIES})
else()
    message(FATAL_ERROR "Boost not found!")
endif()
