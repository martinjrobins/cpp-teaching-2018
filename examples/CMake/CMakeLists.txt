cmake_minimum_required (VERSION 3.6)

# Find third-party dependencies
find_package (Boost REQUIRED mpi)
find_package (MPI REQUIRED)

# compile printRank.cpp into an executable target
add_executable (printRank printRank.cpp)

# add include directories to the target
target_include_directories (printRank PUBLIC ${Boost_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})

# link the boost and mpi libraries to the target
target_link_libraries (printRank PUBLIC ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

if (MPI_COMPILE_FLAGS)
  set_target_properties(printRank PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
  set_target_properties(printRank PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()
