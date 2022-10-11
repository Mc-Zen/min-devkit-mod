# Copyright 2018 The Min-DevKit Authors. All rights reserved.
# Use of this source code is governed by the MIT License found in the License.md file.

set(target min.beat.random)
add_project(${target} SOURCES min.beat.random.cpp)
add_auto_unit_test(${target})