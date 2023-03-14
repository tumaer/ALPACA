message( STATUS "Set performance flags" )

set(STANDARD_FLAGS "-m64 -g -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STANDARD_FLAGS}")

# Make the performance option avaible
if(ALPACA_ENV STREQUAL "LRZ" OR ALPACA_ENV STREQUAL "SNG")
    option(PERFORMANCE "Performance" ON)
else(ALPACA_ENV STREQUAL "LRZ" OR ALPACA_ENV STREQUAL "SNG")
    option(PERFORMANCE "Performance" OFF)
endif(ALPACA_ENV STREQUAL "LRZ" OR ALPACA_ENV STREQUAL "SNG")

# Idependently of environment set the performance flag and machine performance flags
if(PERFORMANCE)
    set(MACHINE_FLAGS "${MACHINE_FLAGS} -D PERFORMANCE")
endif(PERFORMANCE)

# Set additional machine performance flags for the supermuc
set(MACHINE_PERFORMANCE_FLAGS "")
if(ALPACA_ENV STREQUAL "SNG")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MACHINE_PERFORMANCE_FLAGS "${MACHINE_PERFORMANCE_FLAGS} -vecabi=cmdtarget -xHost -qopt-zmm-usage=high -qopt-report=5")
    else(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MACHINE_PERFORMANCE_FLAGS "${MACHINE_PERFORMANCE_FLAGS} -ftree-vectorize -mavx512f -mavx512cd -ffast-math -march=skylake-avx512 -fopt-info-vec ")
    endif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
endif(ALPACA_ENV STREQUAL "SNG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MACHINE_FLAGS}")

# Handle deubg mode
option(DBG "Debug" OFF)
if(DBG)
   set(ALPACA_CXX_FLAGS "-O0")
   set(PACO_CXX_FLAGS "-O0")
else(DBG)
   set(ALPACA_CXX_FLAGS "-O3 -DNDEBUG")
   set(ALPACA_CXX_FLAGS "${ALPACA_CXX_FLAGS} ${MACHINE_PERFORMANCE_FLAGS}")
   set(PACO_CXX_FLAGS "")
endif(DBG)

# Option to allow enabling full symmetry operations
option(SYMMETRY "Symmetry preserving flags" ON)

# Set the floating point flags for intel compiler
set(ALPACA_FLOATING_FLAGS "")
if(ALPACA_ENV STREQUAL "LRZ" OR ALPACA_ENV STREQUAL "SNG")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        # If symmetry is set add precise model to keep floating point operation order (allows exception handling)
        if(SYMMETRY)
            set(ALPACA_FLOATING_FLAGS "-fp-model precise")
        # Otherwise allow full compiler optimization
        else(SYMMETRY)
            set(ALPACA_FLOATING_FLAGS "-fp-model fast=2")
        endif(SYMMETRY)
    endif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
endif(ALPACA_ENV STREQUAL "LRZ" OR ALPACA_ENV STREQUAL "SNG")