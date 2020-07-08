MESSAGE( STATUS "Set warning flags" )

set(STANDARD_WARNINGS "-Wall -pedantic -W -Wformat -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wreturn-type -Wno-unused-function")
if( APPLE )
    set(STANDARD_WARNINGS "${STANDARD_WARNINGS} -mmacosx-version-min=10.15")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STANDARD_WARNINGS}")
