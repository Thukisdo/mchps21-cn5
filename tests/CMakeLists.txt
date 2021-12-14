
# Add tests here

find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(dgbmv_test dgbmv_test.cpp)
target_link_libraries(dgbmv_test GTest::gtest GTest::gtest_main)

gtest_discover_tests(dgbmv_test)
