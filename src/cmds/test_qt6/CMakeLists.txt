cmake_minimum_required(VERSION 3.16)
project(test_qt6 LANGUAGES CXX)

# Find the Qt6 package for Widgets module
find_package(Qt6 REQUIRED COMPONENTS Widgets)

add_executable(test_qt6 main.cpp)

# Link Qt6 Widgets libraries
target_link_libraries(test_qt6 PRIVATE Qt6::Widgets)

# Enable automatic MOC, UIC, and RCC
set_target_properties(test_qt6 PROPERTIES
    AUTOMOC ON
    AUTOUIC ON
    AUTORCC ON
)
