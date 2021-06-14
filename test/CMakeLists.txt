enable_testing()
find_package(GTest REQUIRED)

file(GLOB TEST_SRC CONFIGURE_DEPENDS "*.cpp")
add_executable(run_tests ${TEST_SRC})

target_link_libraries(run_tests ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_PROJECT_NAME}_lib)
