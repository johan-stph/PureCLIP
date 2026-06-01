# Minimum cmake version
cmake_minimum_required (VERSION 3.5...3.31)

# Generate Release builds by default
if (NOT CMAKE_BUILD_TYPE)
    message (STATUS "Setting build type to Release")
    set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type" FORCE)
else()
    message (STATUS "User selected build type: ${CMAKE_BUILD_TYPE}")
endif()

# Name of project and that it is C++ only.
project (PureCLIP CXX)

# Require C++17 explicitly — both compilers on supported platforms support it,
# and it gives consistent behaviour across macOS (AppleClang) and Linux (GCC).
# SeqAn 2.5 is C++14/17 compatible; we opt into 17 for the whole project.
#
# Why not C++20/23?
# SeqAn 2.5 itself is tested against C++23 in its own CI (see seqan/seqan
# .github/workflows/ci_linux.yml) so the library is not the blocker.
# The constraint is on PureCLIP: there is no automated numerical regression
# test suite to verify that the HMM output is bit-for-bit identical after a
# standard bump. C++20 changed several implicit conversion and range rules that
# could silently affect floating-point evaluation order in the density
# functions. C++17 is sufficient for the current codebase and safer until a
# validation dataset + expected-output check is added to CI.
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

# ----------------------------------------------------------------------------
# Apple Clang / arm64 (Apple Silicon) OpenMP support
# Apple's clang does not bundle OpenMP; use libomp installed via Homebrew.
# Install with:  brew install libomp
# ----------------------------------------------------------------------------
if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Try to locate libomp via Homebrew (works for both arm64 and Intel).
    # brew is typically not on cmake's default PATH, so probe known prefixes.
    find_program (BREW_EXEC brew
        PATHS /opt/homebrew/bin /usr/local/bin
        NO_DEFAULT_PATH
    )
    if (BREW_EXEC)
        execute_process (
            COMMAND ${BREW_EXEC} --prefix libomp
            OUTPUT_VARIABLE LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE LIBOMP_BREW_RESULT
        )
    endif()
    # Fall back to well-known Homebrew keg paths if brew lookup failed
    if (NOT LIBOMP_PREFIX OR NOT EXISTS "${LIBOMP_PREFIX}/lib/libomp.dylib")
        foreach (_candidate
                "/opt/homebrew/opt/libomp"
                "/usr/local/opt/libomp")
            if (EXISTS "${_candidate}/lib/libomp.dylib")
                set (LIBOMP_PREFIX "${_candidate}")
                break()
            endif()
        endforeach()
    endif()
    if (EXISTS "${LIBOMP_PREFIX}/lib/libomp.dylib")
        message (STATUS "Apple Clang detected -- using libomp at ${LIBOMP_PREFIX}")
        set (OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set (OpenMP_CXX_LIB_NAMES "omp")
        set (OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LIBOMP_PREFIX}/lib -lomp")
    else()
        message (FATAL_ERROR
            "Apple Clang requires libomp for OpenMP support.\n"
            "Install it with:  brew install libomp\n"
            "Then re-run cmake."
        )
    endif()
endif()

find_package (OpenMP REQUIRED)
include (SeqAn)
include (Boost)
find_package (GSL REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# SeqAn defaults BGZF I/O to 16 threads; PureCLIP manages its own
# threading via -nt and requires serial BAM access for correct results.
# This replicates the single change made in the original skrakau/seqan fork.
set (PURECLIP_COMPILE_DEFS
    SEQAN_BGZF_NUM_THREADS=1
    ${SEQAN_DEFINITIONS}
)

# target_compile_options requires a proper cmake list (one flag per element).
# Both SEQAN_CXX_FLAGS and OpenMP_CXX_FLAGS are space-separated strings so we
# split them before passing to target_compile_options.
separate_arguments (SEQAN_CXX_FLAGS_LIST    UNIX_COMMAND "${SEQAN_CXX_FLAGS}")
separate_arguments (OpenMP_CXX_FLAGS_LIST   UNIX_COMMAND "${OpenMP_CXX_FLAGS}")
set (PURECLIP_COMPILE_OPTIONS ${SEQAN_CXX_FLAGS_LIST} ${OpenMP_CXX_FLAGS_LIST})

set (PURECLIP_INCLUDE_DIRS
    ${SEQAN_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# PureCLIP static build for Linux
set (LINK_TYPE "-dynamic")
if (STATIC_BUILD)
    set (LINK_TYPE "-static")
    if (UNIX AND NOT APPLE)
        message (STATUS "Configuring Linux static build")
        set (CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set (CMAKE_EXE_LINKER_FLAGS "-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    elseif (UNIX AND APPLE)
        set (LINK_TYPE "")
        message (STATUS "Configuring Mac build with static libstdc++")
        string (REPLACE "dylib" "a" GSL_LIBRARIES "${GSL_LIBRARIES}")
        set (CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        set (BUILD_SHARED_LIBRARIES OFF)
    endif()
endif()

# Executables
add_executable (pureclip
    pureclip.cpp
    util.h
    call_sites.h
    parse_alignments.h
    prepro_mle.h
    hmm_1.h
    density_functions.h)

add_executable (winextract winextract.cpp)

foreach (_target pureclip winextract)
    target_include_directories (${_target} SYSTEM PRIVATE ${PURECLIP_INCLUDE_DIRS})
    target_compile_definitions (${_target} PRIVATE ${PURECLIP_COMPILE_DEFS})
    target_compile_options     (${_target} PRIVATE ${PURECLIP_COMPILE_OPTIONS})
    target_link_libraries      (${_target} PRIVATE
        ${Boost_LIBRARIES}
        ${SEQAN_LIBRARIES}
        ${GSL_LIBRARIES})
endforeach()

# Installation
if (PKG_BUILD)
    set (BINDIR ".")
    set (EXTRADIR ".")
else()
    set (BINDIR "bin")
    set (EXTRADIR "share/pureclip")
endif()
install (TARGETS pureclip winextract RUNTIME DESTINATION "${BINDIR}")
install (FILES "../LICENSE.md" DESTINATION "${EXTRADIR}")

# Packaging
if (UNIX AND NOT APPLE)
    set (CPACK_SYSTEM_NAME "linux64${LINK_TYPE}")
elseif (UNIX AND APPLE)
    set (CPACK_SYSTEM_NAME "mac64${LINK_TYPE}")
endif()
set (CPACK_PACKAGE_VERSION_MAJOR "1")
set (CPACK_PACKAGE_VERSION_MINOR "3")
set (CPACK_PACKAGE_VERSION_PATCH "1")
set (CPACK_GENERATOR "TGZ")
set (CPACK_INSTALLED_DIRECTORIES
    "${CMAKE_CURRENT_SOURCE_DIR}/../util;util"
    "${CMAKE_CURRENT_SOURCE_DIR}/../pkg;.")
message (STATUS "Adding file ${CMAKE_CURRENT_SOURCE_DIR}/../pkg/README")
include (CPack)
