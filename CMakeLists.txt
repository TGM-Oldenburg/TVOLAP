cmake_minimum_required(VERSION 3.0)
project(TVOLAP)

option(BUILD_QT_EXAMPLE "Toggle to build the qt vizualisation example." OFF)

if (BUILD_QT_EXAMPLE)
    #New policy, is a good idea to perform the qt search with nice behaviour
    if(POLICY CMP0020)
        cmake_policy(SET CMP0020 NEW)
    endif(POLICY CMP0020)

    #print qmake path (for user information)
    message(STATUS "qmake excecutabele located in: ${QT_QMAKE_EXECUTABLE}")

    #Include the current directory for the Qt related meta object compilation (MOC)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)

    #Activate MOC for Qt signal and slots (plus other Qt related things)
    set(CMAKE_AUTOMOC ON)
endif(BUILD_QT_EXAMPLE)

#Using c++11 standard
set(CMAKE_CXX_FLAGS "-std=c++0x")

#OS dependent library searches / includes
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    find_library(STDCPP_LIB libstdc++-6.dll)
    find_library(WINPTHREAD_LIB libwinpthread-1.dll)

    if (BUILD_QT_EXAMPLE)
        find_library(QT5CORE_LIB Qt5Core.dll)
        find_library(QT5GUI_LIB Qt5Gui.dll)
        find_library(QT5WIDGETS_LIB Qt5Widgets.dll)
        find_library(QT5PRINT_LIB Qt5PrintSupport.dll)
    endif(BUILD_QT_EXAMPLE)

    if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        message(STATUS "Target is Windows, 64 bit toolchain")
        find_library(SEH_LIB libgcc_s_seh-1.dll)
    else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        message(STATUS "Target is Windows, 32 bit toolchain")
        find_library(DW2_LIB libgcc_s_dw2-1)
    endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    message(STATUS "Target is Linux, 64 bit toolchain")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    message(STATUS "Target is Darwin (Max OS X), 64 bit toolchain")
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    message(SEND_ERROR "Unknown system name (does not equal Windows, Linux or Darwin [Mac OS X])")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#Find qt related things (libraries and include directories)
if (BUILD_QT_EXAMPLE)
    find_package(Qt5Core)
    find_package(Qt5Gui)
    find_package(Qt5Widgets)
    find_package(Qt5PrintSupport)

    #Include all the defined and resolved include paths
    include_directories(${Qt5Core_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5PrintSupport_INCLUDE_PATH})
endif (BUILD_QT_EXAMPLE)

#Set list of source files

set(TVOLAP_SOURCES
    fft.cpp
    fft.h
    TVOLAP.cpp
    TVOLAP.h
    )

set(QT_EXAMPLE_SOURCES
    QtExample.cpp
    QtMainwindow.cpp
    QtMainwindow.h
    qcustomplot.cpp
    qcustomplot.h
    )

set(EXAMPLE_SOURCES
    processExample.cpp
    )

#Add the executable
add_library(TVOLAP SHARED ${TVOLAP_SOURCES})

add_executable(TVOLAP_Example ${EXAMPLE_SOURCES})

#Link system independent required libraries against the VARy executable
target_link_libraries(TVOLAP_Example TVOLAP)

if (BUILD_QT_EXAMPLE)
    add_executable(TVOLAP_QtExample ${QT_EXAMPLE_SOURCES})

    #Link system independent required libraries against the QtExample executable
    target_link_libraries(TVOLAP_QtExample TVOLAP ${Qt5Core_LIBRARIES} ${Qt5Gui_LIBRARIES} ${Qt5Widgets_LIBRARIES} ${Qt5PrintSupport_LIBRARIES})
endif (BUILD_QT_EXAMPLE)

#Copy all related dynamic libraries to the binary folder if we are on windows (so we can start the .exe without external includes)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        file(COPY ${SEH_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        file(COPY ${DW2_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")

    file(COPY ${STDCPP_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${WINPTHREAD_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    #Copy shared qt bibs with the executable, so that it can run standalone
    if (BUILD_QT_EXAMPLE)
        file(COPY ${QT5CORE_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        file(COPY ${QT5GUI_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        file(COPY ${QT5WIDGETS_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        file(COPY ${QT5PRINT_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    endif (BUILD_QT_EXAMPLE)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#Copy the TVOLAP-header so that it lays next to the dynamic library in the binary folder.
file(COPY "TVOLAP.h" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
