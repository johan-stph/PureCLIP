find_package (Boost 1.58)

if (NOT Boost_FOUND)
    set (BOOST_URL    "https://boostorg.jfrog.io/artifactory/main/release/1.64.0/source/boost_1_64_0.tar.bz2")
    set (BOOST_SHA256 "7bcc5caace97baa948931d712ea5f37038dbb1c5d89b43ad4def4ed7cb683332")
    set (BOOST_ZIP_OUT ${CMAKE_CURRENT_BINARY_DIR}/boost_1_64_0.tar.bz2)
    set (BOOST_ROOT    ${CMAKE_CURRENT_BINARY_DIR}/boost_1_64_0)

    if (NOT EXISTS ${BOOST_ROOT})
        message ("Downloading ${BOOST_URL}")
        file (DOWNLOAD ${BOOST_URL} ${BOOST_ZIP_OUT}
            EXPECTED_HASH SHA256=${BOOST_SHA256}
            SHOW_PROGRESS STATUS status)
        list (GET status 0 ret)
        if (NOT ret EQUAL 0)
            message (FATAL_ERROR "Boost download failed")
        endif()
        message ("Unpacking ${BOOST_ZIP_OUT}")
        execute_process (COMMAND ${CMAKE_COMMAND} -E tar xf ${BOOST_ZIP_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        if (EXISTS ${BOOST_ZIP_OUT})
            file (REMOVE ${BOOST_ZIP_OUT})
        endif()
    endif()
    find_package (Boost 1.58 REQUIRED)
endif()
