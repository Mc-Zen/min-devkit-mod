cmake_minimum_required(VERSION 3.19)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Butterfly_Audio_Library VERSION 1.0)

if(PROJECT_IS_TOP_LEVEL)
	option(UNIT_TESTING "Add and build unit tests" OFF)
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()


if(UNIT_TESTING)
	find_package(Catch2 3 REQUIRED)
	enable_testing()
	include(Catch)
endif()


add_subdirectory(src/math)
add_subdirectory(src/wave)
add_subdirectory(src/synth)
add_subdirectory(src/utilities)
add_subdirectory(src/modulation)

set(target Butterfly_Audio_Library)
add_library(${target} INTERFACE)
target_link_libraries(${target} INTERFACE math utilities wave synth synth-modulation)