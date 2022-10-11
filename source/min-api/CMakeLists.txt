# Copyright 2018 The Min-API Authors. All rights reserved.
# Use of this source code is governed by the MIT License found in the License.md file.

cmake_minimum_required(VERSION 3.19)

option(C74_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)


# Add the max sdk, if it exists
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/max-sdk-base/CMakeLists.txt")
    add_subdirectory(max-sdk-base)
else()
    message(FATAL_ERROR "The max-sdk-base repository is missing. Maybe you did not clone recursively?")
endif ()

include(script/add_min_target.cmake)
include(test/add_min_unit_test.cmake)

project(C74MinAPI) # can we rename this?



if (${CMAKE_GENERATOR} MATCHES "Xcode")
 		if (${XCODE_VERSION} VERSION_LESS 10)
  			message(STATUS "Xcode 10 or later is required. Please install from the Mac App Store.")
			return ()
		endif ()
endif ()

if (WIN32)
    add_definitions(
		-DMAXAPI_USE_MSCRT
		-DWIN_VERSION
		-D_USE_MATH_DEFINES
	)
endif ()

file(GLOB_RECURSE MIN_API_HEADERS
	 RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
	 "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)

set(version 0.6.0)
set(target min-api) # uniform naming

add_library(${target} INTERFACE ${MIN_API_HEADERS}
	"script/add_min_target.cmake"
	"script/min-posttarget.cmake"
	"script/min-pretarget.cmake"
	"script/c74_set_target_xcode_warning_flags.cmake"
	"test/add_min_unit_test.cmake"
	"test/min-object-unittest.cmake"
)
add_library(API ALIAS min-api) # alias for backwards compatibility
add_library(min-api-test INTERFACE)
target_include_directories(min-api-test INTERFACE "test")

set_target_properties(${target} PROPERTIES FOLDER "C74_Libraries")
set_target_properties(${target} PROPERTIES VERSION ${version})

target_include_directories(${target} INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/include>"
)


# Create package export
if (PROJECT_IS_TOP_LEVEL)
    option(C74_MINAPI_EXPORT_PACKAGE "Enable package export that can be installed" OFF)
    # How to export the package:
    # --------------------------
    #   mkdir build && cd build
    #   cmake ../ -DC74_MINAPI_EXPORT_PACKAGE=ON
    #   cmake --build .  
    #   cmake --install .  
    #
    # Optionally use --prefix "my/path" in the last command if you don't want to put 
    # it into the default location or can't due to access rights.

    # How to include it in another repository:
    # ----------------------------------------
    # In your cmake file include
    #   find_package(C74MinApi REQUIRED)   
    #   target_link_library(myTarget C74::min-api C74::max-sdk) 
    # or whichever you need. If you didn't install into the default location, add
    #   set(C74MinApi_DIR mypath/lib/cmake/C74MinApi)
    # before the call to find_package. 
    
    # Read more here: https://cmake.org/cmake/help/latest/guide/importing-exporting/index.html


    if (C74_MINAPI_EXPORT_PACKAGE)
        include(GNUInstallDirs)
        include(CMakePackageConfigHelpers)

        set(C74MINAPI_CMAKE_CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/C74MinApi")

        # Create install targets
        install(TARGETS min-api max-sdk-base-headers max-sdk-base
                EXPORT C74Targets
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )

        # Install headers
        install(
            DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/max-sdk-base/c74support"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
            FILES_MATCHING PATTERN "*.h"
        )

        # Install targets
        install(
            EXPORT C74Targets
            FILE C74Targets.cmake
            NAMESPACE C74::
            DESTINATION ${C74MINAPI_CMAKE_CONFIG_DESTINATION}
        )

        # Install CMake scripts
        install(
            FILES
                "max-sdk-base/script/add_max_target.cmake"
                "max-sdk-base/script/macros.cmake"
                "max-sdk-base/script/git-rev.cmake"
                "script/add_min_target.cmake"
                "script/set_target_xcode_warning_flags.cmake"
                "script/min-api.cmake"
                "test/add_min_unit_test.cmake"
                "test/add_min_unit_test.cmake"
          DESTINATION
                ${C74MINAPI_CMAKE_CONFIG_DESTINATION}
        )

        # Package config file
        configure_package_config_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/script/Config.cmake.in"
            "${CMAKE_CURRENT_BINARY_DIR}/C74MinApiConfig.cmake"
            INSTALL_DESTINATION ${C74MINAPI_CMAKE_CONFIG_DESTINATION}
        )

        # Package version file
        write_basic_package_version_file(
          "${CMAKE_CURRENT_BINARY_DIR}/C74MinApiConfigVersion.cmake"
          VERSION "${version}"
          COMPATIBILITY AnyNewerVersion
        )

        # Install package config and version file
        install(FILES
                "${CMAKE_CURRENT_BINARY_DIR}/C74MinApiConfig.cmake"
                "${CMAKE_CURRENT_BINARY_DIR}/C74MinApiConfigVersion.cmake"
            DESTINATION "${C74MINAPI_CMAKE_CONFIG_DESTINATION}"
        )

    endif(C74_MINAPI_EXPORT_PACKAGE)

endif(PROJECT_IS_TOP_LEVEL)


enable_testing()

# min-devkit and probably many user projects add the min-api as a subdirectory,
# however, we want the min-api to be able to act as a standalone project as well.
# so projects including the min-api should use the various .cmake scripts instead
# of including min-api as a subdirectory.
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
else()
    message(DEPRECATION "Adding min-api via add_subdirectory() is deprecated and might break your setup in the future. Please remove the add_subdirectory() line from your CMakeLists.txt file.")
endif()

find_package (Doxygen QUIET)
option (BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if (BUILD_DOCUMENTATION)
    if (NOT DOXYGEN_FOUND)
         message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif ()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    message("Doxygen build started.")

    add_custom_target(doc
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)

    #    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION     share/doc)
endif ()
