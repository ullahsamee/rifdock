MESSAGE( "entering apps/rosetta/riflib")

FILE(GLOB L1 "*.cc")
FILE(GLOB L2 "*/*.cc")

list(APPEND EXTRA_LIBS  )

include_directories("..") # scheme source dir
include_directories(".")
link_directories(${CMAKE_LIBRARY_PATH})


add_library(riflib SHARED ${L1} ${L2} ${L3} ${L4} ${L5} )

install ( TARGETS riflib LIBRARY DESTINATION lib )
