# Supported CMake versions

cmake_minimum_required(VERSION 3.23 FATAL_ERROR)

# Project configuration

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(MatrixPermanent VERSION 0.0.1 LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed; make a build directory and run CMake from it\n")
endif()

if(ENV{PERMANENT_TUNE})
  set(PERMANENT_TUNE ON)
endif()

if(SKBUILD OR ENV{PERMANENT_PYTHON})
  set(PERMANENT_PYTHON ON)
endif()

# Dependencies

if(PROJECT_IS_TOP_LEVEL AND PERMANENT_PYTHON)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
elseif(PERMANENT_TUNE)
  find_package(Python REQUIRED COMPONENTS Interpreter)
endif()

# Tuning header target

set(PERMANENT_TUNING_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/permanent/tuning.h")
set(PERMANENT_TUNING_DATA "${CMAKE_BINARY_DIR}/tuning.csv")

if(PERMANENT_TUNE)
  add_executable(PermanentTune src/tuning.cc)
  target_include_directories(PermanentTune PRIVATE include)

  add_custom_target(TuningData ALL
    COMMAND
      PermanentTune "${PERMANENT_TUNING_DATA}"
    BYPRODUCTS
      "${PERMANENT_TUNING_DATA}")

  add_custom_target(TuningHeader ALL
    COMMAND
    ${Python_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_tuning_header.py"
        "${PERMANENT_TUNING_DATA}"
        "${PERMANENT_TUNING_HEADER}"
    BYPRODUCTS
      "${PERMANENT_TUNING_HEADER}")
  add_dependencies(TuningHeader TuningData)
else()
  add_custom_target(TuningHeader ALL
    COMMAND
      cp "-f"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/permanent/tuning.default.h"
        "${PERMANENT_TUNING_HEADER}"
    BYPRODUCTS
      "${PERMANENT_TUNING_HEADER}")
endif()

# Header library target

add_library(permanent_headers INTERFACE)
add_library(MatrixPermanent::headers ALIAS permanent_headers)
add_dependencies(permanent_headers TuningHeader)
target_include_directories(permanent_headers
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(PROJECT_IS_TOP_LEVEL)
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)
  include(CPack)
  write_basic_package_version_file(
    "${PROJECT_NAME}ConfigVersion.cmake"
    VERSION
        ${PROJECT_VERSION}
    COMPATIBILITY
        SameMajorVersion)
  configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION
      ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(TARGETS permanent_headers
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  install(
    EXPORT
      ${PROJECT_NAME}_Targets
    FILE
      ${PROJECT_NAME}Targets.cmake
    NAMESPACE
      ${PROJECT_NAME}::
    DESTINATION
      ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(
    FILES
      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION
      ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(
    DIRECTORY
      ${PROJECT_SOURCE_DIR}/include/permanent
    DESTINATION
      include)
endif()

# Python library target

if(PROJECT_IS_TOP_LEVEL AND PERMANENT_PYTHON)
  execute_process(
    COMMAND
    ${Python_EXECUTABLE} "-c" "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE
      NumPy_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  python_add_library(permanent MODULE src/main.cc WITH_SOABI)
  target_include_directories(permanent PRIVATE include ${NumPy_INCLUDE_DIRS})
  if (UNIX)
    if (APPLE)
      set_target_properties(permanent
        PROPERTIES
          LINK_FLAGS "-Wl,-dylib,-undefined,dynamic_lookup")
    else()
      set_target_properties(permanent
        PROPERTIES
          LINK_FLAGS "-Wl,--allow-shlib-undefined")
    endif()
  endif()
  install(TARGETS permanent DESTINATION .)
endif()
