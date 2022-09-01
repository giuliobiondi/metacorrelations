cmake_minimum_required(VERSION 2.8.12)
project(metacorrelations)
set(CMAKE_CXX_STANDARD 17)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )
find_package(OpenMP)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3" )

add_executable(${PROJECT_NAME} "main.cpp" "preprocesser.cpp" "measures.cpp" "graphtools.cpp" "Classifier.cpp" "Random.cpp" "de.cpp" "Timer.cpp" "Loader.cpp" "Filesystem.cpp" "metrics.cpp" "pathManager.cpp")
include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/snap/snap-core" "${CMAKE_SOURCE_DIR}/snap/glib-core")
add_custom_target(
    snap_lib
    COMMAND make -C snap-core lib
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/snap/
)
add_dependencies(${PROJECT_NAME} snap_lib)
target_link_libraries(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/snap/snap-core/libsnap.a)