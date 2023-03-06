cmake_minimum_required(VERSION 3.16.3)
set(CMAKE_CXX_STANDARD 20)
project(Thicken)


include_directories("/opt/homebrew/include")

include_directories(/usr/local/include/osqp/)

include_directories(/opt/homebrew/include/eigen3)

include_directories(/opt/osqp-eigen/include/)

include_directories(vcglib/)

LINK_DIRECTORIES(/usr/local/lib/)

find_package(OsqpEigen REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message(STATUS "Release版本")
    set(BuildType "Release")
    set(CMAKE_C_FLAGS "$ENV{CFLAGS} -DNODEBUG -O3 -Wall")
    set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -DNODEBUG -O3 -Wall")
else()

    message(STATUS "Debug版本")
    set(BuildType "Debug")
    set(CMAKE_C_FLAGS "$ENV{CFLAGS} -O0 -g -Wall")
    set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -O0 -g -Wall")
endif()

file (GLOB CppFile
        MeshKernel/*
        )

add_executable(Thicken main.cpp CGALPolygon.cpp DSU.cpp remeshing.cpp MeshBuilder.cpp ${CppFile} vcglib/wrap/ply/plylib.cpp)

if(APPLE)
    target_link_libraries(Thicken
            /opt/homebrew/lib/libgmp.dylib
            -lpthread
            /opt/homebrew/lib/libmpfr.dylib
            ${Boost_LIBRARY_DIRS}
            OsqpEigen::OsqpEigen
            /opt/homebrew/lib/libgflags.dylib
            )
else()
    target_link_libraries(Thicken
            -lgmp
            -lpthread
            -lmpfr
            ${Boost_LIBRARY_DIRS}
            OsqpEigen::OsqpEigen
            -lgflags
            )
endif()
