cmake_minimum_required(VERSION 3.14)

project(MyART
  VERSION
    1.0
  DESCRIPTION
    "This small project is used to practice High-Performance Computing programming"
  LANGUAGES
    CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(TIFF REQUIRED)
find_package(OpenMP)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(myART include/Image.h include/Image.inl src/Image.cxx src/myART.cxx)

target_include_directories(myART PRIVATE ${TIFF_INCLUDE_DIRS} include)

target_link_libraries(myART ${TIFF_LIBRARIES})

set_target_properties(myART PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED ON
                      CXX_EXTENSIONS OFF
)
