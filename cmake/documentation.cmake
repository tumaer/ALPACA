# check if Doxygen is installed
find_package(Doxygen)
if(DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    MESSAGE( STATUS "Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( documentation
	COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation of ALPACA with Doxygen"
        VERBATIM )
    set_target_properties( documentation PROPERTIES EXCLUDE_FROM_ALL TRUE )
else(DOXYGEN_FOUND)
  MESSAGE( AUTHOR_WARNING "Doxygen need to be installed to generate the doxygen documentation")
endif(DOXYGEN_FOUND)
