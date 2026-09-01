if (SEQAN_ROOT)
	if (NOT EXISTS ${SEQAN_ROOT}/share/cmake/seqan/seqan-config.cmake)
		message (FATAL_ERROR "SEQAN_ROOT was specified but '${SEQAN_ROOT}/share/cmake/seqan/seqan-config.cmake' does not exist.")
	endif()
else()
	set ( SEQAN_URL    "https://github.com/seqan/seqan/releases/download/seqan-v2.5.3/seqan-library-2.5.3.zip")
	set ( SEQAN_SHA256 "7da029319f0d0674f5aa1a3939b4e7b89d2a9a8e7a288a182492347ce3b5db3d")
	set ( SEQAN_ZIP_OUT ${CMAKE_CURRENT_BINARY_DIR}/seqan-library-2.5.3.zip )
	set ( SEQAN_ROOT    ${CMAKE_CURRENT_BINARY_DIR}/seqan-library-2.5.3 )

	if (NOT EXISTS ${SEQAN_ROOT}/share/cmake/seqan/seqan-config.cmake)
		# Download zip
		message ("Downloading ${SEQAN_URL}")
		file (DOWNLOAD ${SEQAN_URL} ${SEQAN_ZIP_OUT}
			EXPECTED_HASH SHA256=${SEQAN_SHA256}
			SHOW_PROGRESS STATUS status)
		list (GET status 0 ret)
		if (NOT ret EQUAL 0)
			message (FATAL_ERROR "SeqAn download failed")
		endif()
		# Unpack zip
		message ("Unpacking ${SEQAN_ZIP_OUT}")
		execute_process (COMMAND ${CMAKE_COMMAND} -E tar xf ${SEQAN_ZIP_OUT}
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
		# Remove zip
		if (EXISTS ${SEQAN_ZIP_OUT})
			file (REMOVE ${SEQAN_ZIP_OUT})
		endif()
	endif()

	if (NOT EXISTS ${SEQAN_ROOT}/share/cmake/seqan/seqan-config.cmake)
		message (FATAL_ERROR "Failed to download and unpack '${SEQAN_URL}'")
	endif()
endif ()

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB REQUIRED)

# Load SeqAn 2.5.3 via its CMake config file (uses execute_process, cmake 4.x safe).
set (SeqAn_DIR ${SEQAN_ROOT}/share/cmake/seqan)
find_package (SeqAn CONFIG REQUIRED)
