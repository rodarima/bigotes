# Copyright (c) 2024 Barcelona Supercomputing Center (BSC)
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.20)

project(BIGOTES LANGUAGES C)

add_compile_options(
  -Wall
  -Wextra
  -Wformat
  -Wmissing-prototypes
  -Wstrict-prototypes
  -Wold-style-definition
  -pedantic
  -Werror
)

set(CMAKE_C_VISIBILITY_PRESET hidden)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED TRUE)
set(CMAKE_C_EXTENSIONS FALSE)

add_definitions(-D_POSIX_C_SOURCE=200809L)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Asan UBsan." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

include(GNUInstallDirs)

# Extra build type for AddressSanitizer (Asan)
set(CMAKE_C_FLAGS_ASAN "${CMAKE_C_FLAGS_DEBUG} \
  -fsanitize=address \
  -fno-optimize-sibling-calls \
  -fsanitize-address-use-after-scope \
  -fno-omit-frame-pointer"
  CACHE STRING "Flags used by the C compiler during AddressSanitizer builds." FORCE)

set(CMAKE_C_FLAGS_UBSAN "${CMAKE_C_FLAGS_DEBUG} -fsanitize=undefined"
  CACHE STRING "Flags used by the C compiler during UndefinedBehaviorSanitizer builds." FORCE)

# Required for clock_gettime() in glibc <= 2.17
include(CheckLibraryExists)
check_library_exists(c clock_gettime "" HAVE_CLOCK_GETTIME)

if (NOT HAVE_CLOCK_GETTIME)
  check_library_exists(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME_RT)
  if (HAVE_CLOCK_GETTIME_RT)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lrt")
    set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES} -lrt")
  else()
    message(FATAL_ERROR "cannot find clock_gettime()")
  endif()
endif()

add_executable(bigotes bigotes.c pnorm.c qnorm.c swilk.c dip.c stats.c common.c)
target_link_libraries(bigotes PRIVATE m)
install(TARGETS bigotes RUNTIME DESTINATION bin)
install(FILES bigotes.1 DESTINATION "${CMAKE_INSTALL_MANDIR}/man1")
