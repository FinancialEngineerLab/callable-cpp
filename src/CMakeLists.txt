file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true "../include/*.hpp" "*.cpp")
set(SOURCES ${SOURCES})

add_library(${CMAKE_PROJECT_NAME}_lib STATIC ${SOURCES})
