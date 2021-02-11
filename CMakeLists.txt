cmake_minimum_required( VERSION 3.10 FATAL_ERROR )

set( ALPACA_ENV "$ENV{ALPACA_ENVIRONMENT}" )
# Set the compilers to be used to build the targets.
# Note, this has to be done before defining the project.
INCLUDE("./cmake/compiler.cmake")

PROJECT( ALPACA )
set( CMAKE_CXX_STANDARD 20 )

INCLUDE( FindMPI REQUIRED )
include_directories( ${MPI_CXX_INCLUDE_DIRS} )

INCLUDE( FindHDF5 REQUIRED )
include_directories( ${HDF5_INCLUDE_DIRS} )

# Define the present git hash as compiler macro.
INCLUDE("./cmake/git_hash.cmake")
# Define warning flags used for compilation.
INCLUDE("./cmake/warning_flags.cmake")
# Define performance flags used for compilation (machine dependent).
INCLUDE("./cmake/performance_flags.cmake")
# Define an option to chosse the dimension of the build.
INCLUDE("./cmake/dimension.cmake")
# Define a target to create the doxygen documentation.
INCLUDE("./cmake/documentation.cmake")

# Include directories to know all necessary ALPACA headers.
INCLUDE_DIRECTORIES(src)
# Include directories of third party libraries
INCLUDE_DIRECTORIES(3rdParty/expression_toolkit)
INCLUDE_DIRECTORIES(3rdParty/tiny_xml)

file(GLOB_RECURSE SOURCE_FILES "src/*.cpp")
list(APPEND SOURCE_FILES "./3rdParty/tiny_xml/tinyxml2.cpp")
file(GLOB_RECURSE SOURCE_FILES_FOR_LIB "src/*.cpp")

# Exclude user_expression.cpp to compile it as a separate library
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*user_expression.cpp$")
list(FILTER SOURCE_FILES_FOR_LIB EXCLUDE REGEX ".*user_expression.cpp$")
list(FILTER SOURCE_FILES_FOR_LIB EXCLUDE REGEX ".*main.cpp$")

include(CheckIPOSupported)
check_ipo_supported( RESULT IPOPOSSIBLE )
option(DBG "Debug" OFF)

# Define the ALPACA executable.
add_executable(ALPACA ${SOURCE_FILES})
# Define the ALPACA library.
add_library(ALPACAlib STATIC ${SOURCE_FILES_FOR_LIB} library/alpaca_runner.h library/alpaca_runner.cpp)
# Define Python module
option(PYMODULE "Python Module" OFF)
if( PYMODULE )
   add_subdirectory(3rdParty/pybind11)
   INCLUDE_DIRECTORIES(3rdParty/pybind11)
   pybind11_add_module(alpacapy ${SOURCE_FILES_FOR_LIB} library/python_module.cpp library/alpaca_runner.h library/alpaca_runner.cpp)
endif( PYMODULE )

# Define a target for unit tests.
file(GLOB_RECURSE TEST_FILES "test/*.cpp")
list(FILTER SOURCE_FILES EXCLUDE REGEX "src/main.cpp$")
list(APPEND TEST_FILES "${SOURCE_FILES}")
add_executable(Paco EXCLUDE_FROM_ALL ${TEST_FILES})

##Include testing libraries into Paco
# First build target for approval tests
add_subdirectory(3rdParty/ApprovalTests.cpp EXCLUDE_FROM_ALL )
target_link_libraries( Paco ApprovalTests::ApprovalTests )

# Then simple includes for the res
target_include_directories( Paco PRIVATE 3rdParty/FakeIt/single_header/catch)
target_include_directories( Paco PRIVATE 3rdParty/Catch2/single_include )
target_include_directories( Paco PRIVATE test )

if( IPOPOSSIBLE AND NOT DBG )
   set_property(TARGET ALPACA PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
   set_property(TARGET ALPACAlib PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
   if( PYMODULE )
      set_property(TARGET alpacapy PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
   endif( PYMODULE )
endif( IPOPOSSIBLE AND NOT DBG )
set_target_properties(ALPACA PROPERTIES COMPILE_FLAGS "${ALPACA_CXX_FLAGS} ${ALPACA_FLOATING_FLAGS}")
set_target_properties(ALPACAlib PROPERTIES COMPILE_FLAGS "${ALPACA_CXX_FLAGS} ${ALPACA_FLOATING_FLAGS}")
if( PYMODULE )
   set_target_properties(alpacapy PROPERTIES COMPILE_FLAGS "${ALPACA_CXX_FLAGS} ${ALPACA_FLOATING_FLAGS}")
endif( PYMODULE )
set( ALLUSEDFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_COMPILE_FLAGS} ${CMAKE_CXX_LINK_FLAGS} ${ALPACA_CXX_FLAGS}" )
if( IPOPOSSIBLE AND NOT DBG )
   set( ALLUSEDFLAGS "${ALLUSEDFLAGS} ${CMAKE_CXX_COMPILE_OPTIONS_IPO}" )
endif( IPOPOSSIBLE AND NOT DBG )
message( "-- Alpaca compiles and links using the following flags (without warranty, check VERBOSE build)\n-- ${ALLUSEDFLAGS}" )
message( "-- Alpaca compiles and links using the following floating flags:" )
message("${ALPACA_FLOATING_FLAGS}")
target_compile_definitions(ALPACA PUBLIC TEST_VIRTUAL=)
target_compile_definitions(ALPACAlib PUBLIC TEST_VIRTUAL=)
if( PYMODULE )
   target_compile_definitions(alpacapy PRIVATE TEST_VIRTUAL=)
endif( PYMODULE )

target_link_libraries( ALPACA ${MPI_CXX_LIBRARIES} )
target_link_libraries( ALPACA ${HDF5_LIBRARIES} )
target_link_libraries( ALPACAlib ${HDF5_LIBRARIES} )
if( PYMODULE )
   target_link_libraries( alpacapy PRIVATE ${HDF5_LIBRARIES} )
endif( PYMODULE )

if( MPI_CXX_COMPILE_FLAGS )
  set_target_properties( ALPACA PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif( MPI_CXX_COMPILE_FLAGS )
if( MPI_CXX_LINK_FLAGS )
  set_target_properties( ALPACA PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif( MPI_CXX_LINK_FLAGS )


# Define the UserExpressions executable, encapsulating the expression toolkit from the rest of the program
find_library(UserExpr NAMES UserExpressions PATHS ${CMAKE_CURRENT_BINARY_DIR})
if(NOT UserExpr)
   MESSAGE( STATUS "UserExpressions not found - will be build" )
   add_library(UserExpressions STATIC src/user_expression.cpp )
   if( IPOPOSSIBLE AND NOT DBG )
      set_property(TARGET UserExpressions PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
   endif( IPOPOSSIBLE AND NOT DBG )
   target_link_libraries(ALPACA UserExpressions)
   target_link_libraries(ALPACAlib UserExpressions)
   if( PYMODULE )
      target_link_libraries(alpacapy PRIVATE UserExpressions)
   endif( PYMODULE )
   target_link_libraries(Paco UserExpressions)
else(NOT UserExpr)
   MESSAGE( STATUS "UserExpressions found at ${UserExpr}" )
   target_link_libraries(ALPACA ${UserExpr})
   target_link_libraries(ALPACAlib ${UserExpr})
   if( PYMODULE )
      target_link_libraries(alpacapy PRIVATE ${UserExpr})
   endif( PYMODULE )
   target_link_libraries(Paco ${UserExpr})
endif(NOT UserExpr)

set_target_properties(ALPACAlib PROPERTIES EXCLUDE_FROM_ALL TRUE)
if( PYMODULE )
   set_target_properties(alpacapy PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif( PYMODULE )

# Use additional flags for test coverage extraction
set(PACO_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage")
if(UNIX AND NOT APPLE)
   set(GCOV_LIBRARY "-lgcov")
else(UNIX AND NOT APPLE)
   set(GCOV_LIBRARY "--coverage")
endif(UNIX AND NOT APPLE)

set_target_properties(Paco PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(Paco PROPERTIES COMPILE_FLAGS "${PACO_CXX_FLAGS}")
target_compile_definitions(Paco PUBLIC TEST_VIRTUAL=virtual)

target_link_libraries( Paco ${MPI_CXX_LIBRARIES} )
target_link_libraries( Paco ${HDF5_LIBRARIES} )
target_link_libraries( Paco ${GCOV_LIBRARY} )

if( MPI_CXX_COMPILE_FLAGS )
  set_target_properties( Paco PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif( MPI_CXX_COMPILE_FLAGS )
if( MPI_CXX_LINK_FLAGS )
  set_target_properties( Paco PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif( MPI_CXX_LINK_FLAGS )

install(TARGETS ALPACAlib DESTINATION lib)
install(FILES library/alpaca_runner.h DESTINATION include)

