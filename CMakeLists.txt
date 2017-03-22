cmake_minimum_required(VERSION 3.0)
project(sanityCheckWOLAP)

#New policy, is a good idea to perform the qt search with nice behaviour
if(POLICY CMP0020)
cmake_policy(SET CMP0020 NEW)
endif()

#Using c++11 standard
set(CMAKE_CXX_FLAGS "-std=c++0x")

#Include the current directory because of MOC-compilation (Qt related code wrapper)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#Activate MOC compilation for Qt signal and slots (plus other Qt related things)
set(CMAKE_AUTOMOC ON)

#OS dependent library searches / includes
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    find_library(STDCPP_LIB libstdc++-6.dll)
    find_library(WINPTHREAD_LIB libwinpthread-1.dll)
    find_library(QT5CORE_LIB Qt5Core.dll)
    find_library(QT5GUI_LIB Qt5Gui.dll)
    find_library(QT5WIDGETS_LIB Qt5Widgets.dll)
    find_library(QT5PRINT_LIB Qt5PrintSupport.dll)
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
find_package(Qt5Core)
find_package(Qt5Gui)
find_package(Qt5Widgets)
find_package(Qt5PrintSupport)

#Include all the defined and resolved include paths
include_directories(${Qt5Core_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5PrintSupport_INCLUDE_PATH})

#Set list of source files (.h) gets included automatically
set(SOURCES
    main.cpp
    mainWindow.cpp
    qcustomplot.cpp
    fft.cpp
    WOLAP.cpp
    )

#Add the executable
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_executable(${PROJECT_NAME} WIN32 ${SOURCES})
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_executable(${PROJECT_NAME} ${SOURCES})
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#Link system independent required libraries against the VARy executable
target_link_libraries(${PROJECT_NAME} ${Qt5Core_LIBRARIES} ${Qt5Gui_LIBRARIES} ${Qt5Widgets_LIBRARIES} ${Qt5PrintSupport_LIBRARIES})

#Copy all related dynamic libraries to the binary folder if we are on windows (so we can start the .exe without external includes)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        file(COPY ${SEH_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        file(COPY ${DW2_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")

    file(COPY ${STDCPP_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${WINPTHREAD_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    #If we build a release, copy shared qt bibs with the executable, so that it can run standalone
    if(CMAKE_BUILD_TYPE MATCHES Release)
        file(COPY ${QT5CORE_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        file(COPY ${QT5GUI_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        file(COPY ${QT5WIDGETS_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        file(COPY ${QT5PRINT_LIB} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    endif(CMAKE_BUILD_TYPE MATCHES Release)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
