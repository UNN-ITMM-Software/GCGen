set(target "${PROJECT_LIBRARY}")

file(GLOB hdrs "${PROJECT_INCLUDE}/*.h*")
file(GLOB srcs "${PROJECT_SOURCE}/*.cpp")

add_library(${target} STATIC ${srcs} ${hdrs})
target_link_libraries(${target} ${LIBRARY_DEPS})

# REPORT
message( STATUS "")
message( STATUS "General configuration for lib_${target}")
message( STATUS "OK")
message( STATUS "")
