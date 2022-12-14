# Copyright 2018 The Min-DevKit Authors. All rights reserved.
# Use of this source code is governed by the MIT License found in the License.md file.

cmake_minimum_required(VERSION 3.0)

set(C74_MIN_API_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../min-api)
include(${C74_MIN_API_DIR}/script/min-pretarget.cmake)


#############################################################
# MAX EXTERNAL
#############################################################


include_directories( 
	"${C74_INCLUDES}"
)


set( SOURCE_FILES
	${PROJECT_NAME}.cpp
)


if (APPLE)
	set(SOURCE_FILES ${SOURCE_FILES} min.environment.mac.cpp)
else ()
	set(SOURCE_FILES ${SOURCE_FILES} min.environment.windows.cpp)
endif ()


add_library( 
	${PROJECT_NAME} 
	MODULE
	${SOURCE_FILES}
)


if (APPLE)
    find_library(CARBON_LIBRARY "Carbon")
    target_link_libraries(${PROJECT_NAME} PUBLIC ${CARBON_LIBRARY})

	find_library(IOKIT_LIBRARY "IOKit")
	target_link_libraries(${PROJECT_NAME} PUBLIC ${IOKIT_LIBRARY})
else()
	set(IPHLPAPI_LIBRARY IPHLPAPI.lib)
	target_link_libraries(${PROJECT_NAME} PUBLIC ${IPHLPAPI_LIBRARY})
endif()

include(${C74_MIN_API_DIR}/script/min-posttarget.cmake)


#############################################################
# UNIT TEST
#############################################################

include(${C74_MIN_API_DIR}/test/min-object-unittest.cmake)
