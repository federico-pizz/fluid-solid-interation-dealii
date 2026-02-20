
cmake_minimum_required(VERSION 3.13.4)

project(FSI_NMPDE LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release builds.
if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
else()
  if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    message(WARNING "Forcing CMAKE_BUILD_TYPE=Release (was ${CMAKE_BUILD_TYPE})")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()
message(STATUS)
message(STATUS "Build type (forced): ${CMAKE_BUILD_TYPE}")
message(STATUS)

add_definitions(-DBUILD_TYPE_RELEASE)
# Release flags: optimize and disable asserts via NDEBUG.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

# Locate MPI compiler.
find_package(MPI REQUIRED)

# Locate Boost.
find_package(Boost 1.72.0 REQUIRED
  COMPONENTS filesystem iostreams serialization
  HINTS ${BOOST_DIR} $ENV{BOOST_DIR} $ENV{mkBoostPrefix})
message(STATUS "Using the Boost-${Boost_VERSION} configuration found at ${Boost_DIR}")
message(STATUS)
include_directories(${Boost_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src)

# Locate deal.II and initialize its variables.
find_package(deal.II 9.3.1 REQUIRED
  HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR} $ENV{mkDealiiPrefix})
deal_ii_initialize_cached_variables()

# Add useful compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-conversion -Wmissing-braces -Wnon-virtual-dtor")

# Build object library to compile core sources.
# add_library(FSI_lib OBJECT
#   src/FSI.cpp
# )
# deal_ii_setup_target(FSI_lib)
# target_link_libraries(FSI_lib Boost::filesystem Boost::iostreams Boost::serialization MPI::MPI_CXX)

# Serial executable.
# add_executable(simulate_serial src/simulate_serial.cpp)
# target_link_libraries(simulate_serial FSI_lib)
# deal_ii_setup_target(simulate_serial)


# Build object library to compile parallel sources.
add_library(FSIParallel_lib OBJECT
  src/FSIParallel.cpp
)
deal_ii_setup_target(FSIParallel_lib)
target_link_libraries(FSIParallel_lib Boost::filesystem Boost::iostreams Boost::serialization MPI::MPI_CXX)

# Parallel executable.
add_executable(simulate_parallel src/simulate_parallel.cpp)
target_link_libraries(simulate_parallel FSIParallel_lib)
deal_ii_setup_target(simulate_parallel)

message(STATUS "Configured FSI project with deal.II")
