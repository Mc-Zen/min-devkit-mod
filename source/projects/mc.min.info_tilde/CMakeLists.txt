# Copyright 2018 The Min-DevKit Authors. All rights reserved.
# Use of this source code is governed by the MIT License found in the License.md file.


set(target mc.min.info_tilde)
add_project(${target} SOURCES mc.min.info_tilde.cpp asd.cpp)	


#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED True)
#add_executable(${target} mc.min.info_tilde.cpp)
#target_link_libraries(${target} PUBLIC min-api max-sdk-base)


include("${CMAKE_CURRENT_SOURCE_DIR}/../../min-api/test/min-object-unittest.cmake")
