cmake_minimum_required(VERSION 3.16)

project(Qt6_SimpleApp LANGUAGES CXX)

# Find the Qt6 package for Widgets module
find_package(Qt6 REQUIRED COMPONENTS Widgets)

add_executable(Qt6_SimpleApp main.cpp)

# Link Qt6 Widgets libraries
target_link_libraries(Qt6_SimpleApp PRIVATE Qt6::Widgets)

# Enable automatic MOC, UIC, and RCC
set_target_properties(Qt6_SimpleApp PROPERTIES
    AUTOMOC ON
    AUTOUIC ON
    AUTORCC ON
)
