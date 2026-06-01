# Minimum cmake version
cmake_minimum_required (VERSION 3.5...3.31)

# Generate Release builds by default
if ( NOT CMAKE_BUILD_TYPE )
    message ( STATUS "Setting build type to Release" )
    set ( CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type" FORCE )
else()
    message ( STATUS "User selected build type: ${CMAKE_BUILD_TYPE}" )
endif()

# Name of project and that it is C++ only.
project (PureCLIP CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

LIST ( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ )

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

find_package ( OpenMP REQUIRED )
include ( SeqAn )
include ( Boost )
find_package (GSL REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# SeqAn defaults BGZF I/O to 16 threads; PureCLIP manages its own
# threading via -nt and requires serial BAM access for correct results.
add_definitions (-DSEQAN_BGZF_NUM_THREADS=1)

# Add include directories.
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS})
include_directories (SYSTEM ${GSL_INCLUDE_DIRS})
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# PureCLIP static build for Linux
set (LINK_TYPE "-dynamic")
if ( STATIC_BUILD )
    set (LINK_TYPE "-static")
    if(UNIX AND NOT APPLE)
        message ( STATUS "Configuring Linux static build" )
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    elseif (UNIX AND APPLE)
        set (LINK_TYPE "")
        message ( STATUS "Configuring Mac build with static libstdc++" )
        string (REPLACE "dylib" "a" GSL_LIBRARIES "${GSL_LIBRARIES}")
        SET(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        SET(BUILD_SHARED_LIBRARIES OFF)
    endif()
endif()

# Add executable and link against SeqAn dependencies.
add_executable (pureclip
                    pureclip.cpp
                    util.h
                    call_sites.h
                    parse_alignments.h
                    prepro_mle.h
                    hmm_1.h
                    density_functions.h)

add_executable (winextract winextract.cpp)

target_link_libraries (pureclip ${Boost_LIBRARIES} ${SEQAN_LIBRARIES} ${GSL_LIBRARIES})
target_link_libraries (winextract ${Boost_LIBRARIES} ${SEQAN_LIBRARIES} ${GSL_LIBRARIES})

# Installation
if ( PKG_BUILD )
    SET ( BINDIR "." )
    SET ( EXTRADIR "." )
else()
    SET ( BINDIR "bin" )
    SET ( EXTRADIR "share/pureclip" )
endif()
install (TARGETS pureclip winextract RUNTIME DESTINATION "${BINDIR}")
install (FILES "../LICENSE.md" DESTINATION "${EXTRADIR}")

# Packaging
if(UNIX AND NOT APPLE)
    SET(CPACK_SYSTEM_NAME "linux64${LINK_TYPE}")
elseif (UNIX AND APPLE)
    SET(CPACK_SYSTEM_NAME "mac64${LINK_TYPE}")
endif()
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "3")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_GENERATOR "TGZ")
SET(CPACK_INSTALLED_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../util;util" "${CMAKE_CURRENT_SOURCE_DIR}/../pkg;.")
message(STATUS "Adding file ${CMAKE_CURRENT_SOURCE_DIR}/../pkg/README")
include ( CPack )
