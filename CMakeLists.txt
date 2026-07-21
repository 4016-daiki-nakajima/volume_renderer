set(CMAKE_OSX_ARCHITECTURES arm64)
cmake_minimum_required(VERSION 3.12)

project(lajolla)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

if(MSVC)
  add_compile_options(/Wall)
else()
  add_compile_options(-Wall)
endif()

find_package(embree REQUIRED)
include_directories(${EMBREE_INCLUDE_PATH})
# include_directories(embree/embree4)
link_libraries(${EMBREE_LIBRARY})
# link_libraries(embree/build/libembree4.4.dylib)

# set(EMBREE_INCLUDE_PATH "embree/include")
# include_directories(${EMBREE_INCLUDE_PATH})

# set(EMBREE_LIBRARY "embree/build")
# link_directories(${EMBREE_LIBRARY})

include_directories(${CMAKE_SOURCE_DIR}/src)

file(GLOB_RECURSE SRCS CONFIGURE_DEPENDS
     src/*.cpp
     src/*.c
     src/*.h
     src/*.hpp
     src/*.inl)
list(FILTER SRCS EXCLUDE REGEX "src/tests/.*")
list(REMOVE_ITEM SRCS ${CMAKE_SOURCE_DIR}/src/main.cpp)

add_library(lajolla_lib STATIC ${SRCS})
add_executable(lajolla src/main.cpp)
target_link_libraries(lajolla lajolla_lib)
if(MSVC)
  add_custom_command(TARGET lajolla POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
      "${PROJECT_SOURCE_DIR}/embree/bin"
      $<TARGET_FILE_DIR:lajolla>)
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(lajolla Threads::Threads)

enable_testing()

add_executable(test_filter src/tests/filter.cpp)
target_link_libraries(test_filter lajolla_lib)
add_test(filter test_filter)
set_tests_properties(filter PROPERTIES PASS_REGULAR_EXPRESSION "SUCCESS")

add_executable(test_frame src/tests/frame.cpp)
target_link_libraries(test_frame lajolla_lib)
add_test(frame test_frame)
set_tests_properties(frame PROPERTIES PASS_REGULAR_EXPRESSION "SUCCESS")

add_executable(test_image src/tests/image.cpp)
target_link_libraries(test_image lajolla_lib)
add_test(image test_image)
set_tests_properties(image PROPERTIES PASS_REGULAR_EXPRESSION "SUCCESS")

add_executable(test_intersection src/tests/intersection.cpp)
target_link_libraries(test_intersection lajolla_lib)
add_test(intersection test_intersection)
set_tests_properties(intersection PROPERTIES PASS_REGULAR_EXPRESSION "SUCCESS")

add_executable(test_materials src/tests/materials.cpp)
target_link_libraries(test_materials lajolla_lib)
add_test(materials test_materials)
set_tests_properties(materials PROPERTIES PASS_REGULAR_EXPRESSION "SUCCESS")

add_executable(test_matrix src/tests/matrix.cpp)
target_link_libraries(test_matrix lajolla_lib)
add_test(matrix test_matrix)
set_tests_properties(matrix PROPERTIES PASS_REGULAR_EXPRESSION "SUCCESS")

add_executable(test_mipmap src/tests/mipmap.cpp)
target_link_libraries(test_mipmap lajolla_lib)
add_test(mipmap test_mipmap)
set_tests_properties(mipmap PROPERTIES PASS_REGULAR_EXPRESSION "SUCCESS")
