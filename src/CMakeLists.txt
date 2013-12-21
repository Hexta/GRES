cmake_minimum_required(VERSION 2.8)

SET(QT_USE_QTOPENGL TRUE)
find_package(Qt4 REQUIRED)
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

find_package(OpenGL REQUIRED)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CMAKE_SOURCE_DIR})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -ggdb")

set(${project}_SRC
    etchingMenu.cpp
    functions.cpp
    main.cpp
    maskPreview.cpp
    selectAtomMenu.cpp
    face.cpp
    geometry.cpp
    maskMenu.cpp
    render.cpp
    settings.cpp
)

set(${project}_HEADERS
    etchingMenu.h
    face.h
    maskMenu.h
    maskPreview.h
    render.h
    selectAtomMenu.h
    settings.h
)



set(${project}_RESOURCES
    resources.qrc
)


QT4_WRAP_CPP(${project}_HEADERS_MOC ${${project}_HEADERS})
QT4_ADD_RESOURCES(${project}_RESOURCES_RCC ${${project}_RESOURCES})

add_executable(${project}
    ${${project}_SRC}
    ${${project}_HEADERS_MOC}
    ${${project}_RESOURCES_RCC}
)


target_link_libraries(${project}
    ${QT_LIBRARIES}
    ${OPENGL_LIBRARIES}
)


