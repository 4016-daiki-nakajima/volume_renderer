find_path(EMBREE_INCLUDE_PATH embree4/rtcore.h
  ${CMAKE_SOURCE_DIR}/embree/include
  /usr/include
  /usr/local/include
  /opt/local/include)

find_library(EMBREE_LIBRARY NAMES embree4 PATHS
  ${CMAKE_SOURCE_DIR}/embree/build
  /usr/lib
  /usr/local/lib
  /opt/local/lib)

if (EMBREE_INCLUDE_PATH AND EMBREE_LIBRARY)
  set(EMBREE_FOUND TRUE)
endif ()