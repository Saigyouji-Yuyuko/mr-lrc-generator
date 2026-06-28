cmake_minimum_required(VERSION 3.16)

project(mr_lrc_generator VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MR_LRC_REQUIRE_ISAL "Fail configure when ISA-L is not found" OFF)
find_package(Threads REQUIRED)
find_package(gflags CONFIG QUIET)

if(NOT gflags_FOUND)
    include(FetchContent)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.10 CACHE STRING "" FORCE)
    set(GFLAGS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(GFLAGS_BUILD_gflags_LIB ON CACHE BOOL "" FORCE)
    set(GFLAGS_BUILD_gflags_nothreads_LIB OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        gflags
        GIT_REPOSITORY https://github.com/gflags/gflags.git
        GIT_TAG v2.2.2
    )
    FetchContent_MakeAvailable(gflags)
endif()

set(ISAL_LOCAL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/third_party/isa-l/install/usr")
if(ISAL_INCLUDE_DIR AND NOT EXISTS "${ISAL_INCLUDE_DIR}/erasure_code.h")
    message(WARNING "Ignoring stale ISAL_INCLUDE_DIR: ${ISAL_INCLUDE_DIR}")
    unset(ISAL_INCLUDE_DIR CACHE)
endif()
if(ISAL_LIBRARY AND NOT EXISTS "${ISAL_LIBRARY}")
    message(WARNING "Ignoring stale ISAL_LIBRARY: ${ISAL_LIBRARY}")
    unset(ISAL_LIBRARY CACHE)
endif()
find_path(ISAL_INCLUDE_DIR
    NAMES erasure_code.h
    HINTS
        "${ISAL_LOCAL_PREFIX}/include"
        "${ISAL_LOCAL_PREFIX}/include/isa-l"
    PATH_SUFFIXES isa-l
)
find_library(ISAL_LIBRARY
    NAMES isal
    HINTS
        "${ISAL_LOCAL_PREFIX}/lib"
        "${ISAL_LOCAL_PREFIX}/lib/x86_64-linux-gnu"
)

add_executable(mr-lrc-generator
    src/gf256.cpp
    src/main.cpp
    src/mr_lrc.cpp
)

target_include_directories(mr-lrc-generator PRIVATE src)
target_link_libraries(mr-lrc-generator PRIVATE Threads::Threads)

if(TARGET gflags::gflags)
    target_link_libraries(mr-lrc-generator PRIVATE gflags::gflags)
elseif(TARGET gflags)
    target_link_libraries(mr-lrc-generator PRIVATE gflags)
elseif(TARGET gflags_static)
    target_link_libraries(mr-lrc-generator PRIVATE gflags_static)
else()
    message(FATAL_ERROR "gflags target was not found after configuring gflags")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(mr-lrc-generator PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(ISAL_INCLUDE_DIR AND ISAL_LIBRARY)
    target_include_directories(mr-lrc-generator PRIVATE "${ISAL_INCLUDE_DIR}")
    target_link_libraries(mr-lrc-generator PRIVATE "${ISAL_LIBRARY}")
    target_compile_definitions(mr-lrc-generator PRIVATE MR_LRC_USE_ISAL=1)
    get_filename_component(ISAL_LIBRARY_DIR "${ISAL_LIBRARY}" DIRECTORY)
    set_target_properties(mr-lrc-generator PROPERTIES
        BUILD_RPATH "${ISAL_LIBRARY_DIR}"
        INSTALL_RPATH "${ISAL_LIBRARY_DIR}"
    )
    message(STATUS "Using ISA-L: ${ISAL_LIBRARY}")
else()
    if(MR_LRC_REQUIRE_ISAL)
        message(FATAL_ERROR "ISA-L was not found. Install erasure_code.h/libisal or turn MR_LRC_REQUIRE_ISAL off.")
    endif()
    message(WARNING "ISA-L was not found; using the scalar GF(256) fallback for local builds/tests.")
endif()
