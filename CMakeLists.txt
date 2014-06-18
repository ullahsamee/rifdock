cmake_minimum_required (VERSION 2.6)

add_subdirectory(external/gmock)

include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
set (EXTRA_LIBS ${EXTRA_LIBS} gmock)

include_directories("src")

FILE(GLOB UTIL_CC "src/util/*.cc")
FILE(GLOB NEST_CC "src/nest/*.cc")


project (main_test)
add_executable(main_test src/main/main_test.cc ${UTIL_CC} ${NEST_CC})
target_link_libraries(main_test ${EXTRA_LIBS} )


project (quick_test)
add_executable(quick_test src/main/quick_test.cc)
target_link_libraries(quick_test ${EXTRA_LIBS})
