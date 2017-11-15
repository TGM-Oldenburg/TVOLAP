cmake_minimum_required(VERSION 3.1)
project(TVOLAP)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	cmake_policy(SET CMP0054 OLD)
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

#Using c++11 standard
set(CMAKE_CXX_FLAGS "-std=c++0x")

#OS dependent library searches / includes
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    	find_library(STDCPP_LIB libstdc++-6.dll)
    	find_library(WINPTHREAD_LIB libwinpthread-1.dll)
    endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        message(STATUS "Target is Windows, 64 bit toolchain")
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        	find_library(SEH_LIB libgcc_s_seh-1.dll)
    	endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        message(STATUS "Target is Windows, 32 bit toolchain")
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        	find_library(DW2_LIB libgcc_s_dw2-1)
    	endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    message(STATUS "Target is Linux, 64 bit toolchain")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    message(STATUS "Target is Darwin (Max OS X), 64 bit toolchain")
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    message(SEND_ERROR "Unknown system name (does not equal Windows, Linux or Darwin [Mac OS X])")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#Set list of source files
set(TVOLAP_SOURCES
    fft.cpp
    fft.h
    TVOLAP.cpp
    TVOLAP.h
    )

set(EXAMPLE_SOURCES
    fft.cpp
    fft.h
    testTVOLAP.cpp
    )

#Add the library and executable
add_library(TVOLAP SHARED ${TVOLAP_SOURCES})

add_executable(testTVOLAP ${EXAMPLE_SOURCES})

#Link system independent required libraries against the VARy executable
target_link_libraries(testTVOLAP TVOLAP)

#Copy all related dynamic libraries to the binary folder if we are on windows (so we can start the .exe without external includes)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        	file(COPY ${SEH_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        	file(COPY ${DW2_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")

	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    	file(COPY ${STDCPP_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    	file(COPY ${WINPTHREAD_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
	endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#Copy the TVOLAP-header so that it lays next to the dynamic library in the binary folder.
file(COPY "TVOLAP.h" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
