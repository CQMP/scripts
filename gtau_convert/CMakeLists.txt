cmake_minimum_required(VERSION 2.8)

# global options
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

project(alps_gf)

# dependencies
# 1. lapack
find_package(LAPACK)
link_libraries(${LAPACK_LIBRARIES})

# 2. python (optional)
option(BuildPython "Build python modules" ON)
if (BuildPython) 
    find_package (PythonInterp)
    # tell cmake to find libs for found PythonInterp
    set(Python_ADDITIONAL_VERSIONS "${PYTHON_VERSION_STRING}")
    # There is a bug, involving cmake calling PYTHON framework on APPLE, which always happens to be the system python
    if (APPLE)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print sysconfig.get_config_var(\"LIBDIR\")," OUTPUT_VARIABLE P_LIB OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print sysconfig.get_python_inc()," OUTPUT_VARIABLE P_INC OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(PYTHON_INCLUDE_DIR "${P_INC}")
        set(PYTHON_LIBRARY "${P_LIB}/libpython${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}.dylib") # dirty hack
        unset(P_INC) 
        unset(P_LIB)
    endif(APPLE)

    find_package (PythonLibs REQUIRED)
    include_directories(${PYTHON_INCLUDE_DIRS})
    message(STATUS "Python includes: " ${PYTHON_INCLUDE_DIRS})
    message(STATUS "Python libraries: " ${PYTHON_LIBRARIES})
    link_libraries(${PYTHON_LIBRARIES})
    #set(BOLD_HYB_BUILD_PYTHON True)
    list(APPEND boost_components python)
endif(BuildPython)

# 3. Boost (with/without boost::python)
list(APPEND boost_components filesystem system program_options) 
find_package (Boost REQUIRED COMPONENTS ${boost_components})
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# build library 
include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include)
add_library(alps_gf src/fouriertransform.C)
link_libraries(alps_gf)

# build executables and python
set(execs gtau_convert)
foreach(exec ${execs})
    add_executable(${exec} prog/${exec}.cpp)
    install ( TARGETS ${exec} DESTINATION bin )

    if (BuildPython) 
    # build python module
    add_library(py${exec} MODULE prog/${exec}.cpp)
    set_target_properties(py${exec} PROPERTIES PREFIX "")
    target_compile_definitions(py${exec} PUBLIC -DBUILD_PYTHON_MODULE)
    install ( TARGETS py${exec} DESTINATION lib/python )
    endif()
endforeach()
