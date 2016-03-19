cmake_minimum_required(VERSION 2.8.11)

cmake_policy(SET CMP0020 NEW)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

SET(QT_USE_QTOPENGL TRUE)
set(CMAKE_AUTOMOC ON)

find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)

find_package(OpenGL REQUIRED)

FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CMAKE_SOURCE_DIR})

if(WIN32)
message("Include: ${3rdparty_GL}")
include_directories(SYSTEM ${3rdparty_GL})
endif()


if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -ggdb -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x -ggdb -Wall -Wno-long-long -pedantic -Wextra -Wno-format-extra-args -Wformat=2 -Winit-self -Wmissing-include-dirs -Wsystem-headers -Wfloat-equal -Wconversion -Wlogical-op -Wmissing-noreturn -Wunreachable-code -Winline -Wdisabled-optimization -Woverlength-strings -Wold-style-cast -Woverloaded-virtual -Wsign-promo")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MP")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -ggdb -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x -g -fno-omit-frame-pointer -Wall -Wextra -pedantic -Weverything -Wnosystem-headers")
endif()

set(${PROJECT_NAME}_SRC
    Atoms.cpp
    Cell.cpp
    Coords3D.cpp
    EtchingMenu.cpp
    FloatHelper.cpp
    GeometryPrimitives.cpp
    main.cpp
    MainWindow.cpp
    MaskMenu.cpp
    MaskPreview.cpp
    Render.cpp
    SelectAtomMenu.cpp
    Settings.cpp
    Surface3D.cpp
)

set(${PROJECT_NAME}_HEADERS
    Atoms.h
    Cell.h
    Coords3D.h
    EtchingMenu.h
    FloatHelper.h
    GeometryPrimitives.h
    MainWindow.h
    MaskMenu.h
    MaskPreview.h
    opengl.h
    Render.h
    SelectAtomMenu.h
    Settings.h
	SimulationType.h
    Surface3D.h
)

source_group(${PROJECT_NAME} FILES
    ${${PROJECT_NAME}_SRC}
    ${${PROJECT_NAME}_HEADERS}
)

set(${PROJECT_NAME}_RESOURCES
    resources.qrc
)

QT5_ADD_RESOURCES(${PROJECT_NAME}_RESOURCES_RCC ${${PROJECT_NAME}_RESOURCES})

set(EXECUTABLE_DEPS
    ${${PROJECT_NAME}_SRC}
    ${${PROJECT_NAME}_RESOURCES_RCC}
    ${${PROJECT_NAME}_HEADERS}
)

if(WIN32)
    add_executable(${PROJECT_NAME}
        WIN32
        ${EXECUTABLE_DEPS}
    )
else()
    add_executable(${PROJECT_NAME}
        ${EXECUTABLE_DEPS}
    )
endif()

target_link_libraries(${PROJECT_NAME}
    Qt5::Widgets
	Qt5::OpenGL
    ${OPENGL_LIBRARIES}
)
