cmake_minimum_required(VERSION 3.21)          # 3.21 → first to ship module-preview + good WIN32 fixes
cmake_policy(SET CMP0091 NEW)                 # honour CMAKE_<LANG>_FLAGS_<CONFIG>

# ------------------------------------------------------------------------------
# Project + version
# ------------------------------------------------------------------------------
include(cmake/Version.cmake)                  # sets FASTHAZARD_VERSION_*
project(FastHazard
        VERSION ${FASTHAZARD_VERSION}
        LANGUAGES C CXX ASM)


# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

function(print_target_contents tgt)
    if(NOT TARGET ${tgt})
        message(STATUS "Target '${tgt}' does not exist yet")
        return()
    endif()

    get_target_property(libs  ${tgt} INTERFACE_LINK_LIBRARIES)
    get_target_property(opts  ${tgt} INTERFACE_LINK_OPTIONS)

    message(STATUS "── ${tgt} · INTERFACE_LINK_LIBRARIES")
    foreach(L IN LISTS libs)
        message(STATUS "   ${L}")
    endforeach()

    message(STATUS "── ${tgt} · INTERFACE_LINK_OPTIONS")
    foreach(O IN LISTS opts)
        message(STATUS "   ${O}")
    endforeach()
endfunction()


# ------------------------------------------------------------------------------
# User-facing options
# ------------------------------------------------------------------------------
option(USE_BUILTIN_EXTERNAL_BINS
       "Link against the pre-built static libraries in /external/lib.
        ON by default on Windows, OFF elsewhere."
       ${WIN32})
option(ENABLE_INSTALLER "Generate a platform installer with CPack" ON)

set(FF_HAS_VISUAL         ON  CACHE BOOL "Build visualisation back-end")
set(FF_HAS_SCRIPT         ON  CACHE BOOL "Build AngelScript back-end")
set(FF_HAS_OPENCV         ON  CACHE BOOL "Enable OpenCV features")
set(FF_BUILD_STATIC       ON  CACHE BOOL "Build everything with –static")

# ---------------------------------------------------------------------------
#  Feature flags  (default = ON  →  #define FF_HAS_<FLAG> 1)
# ---------------------------------------------------------------------------
foreach(flag VISUAL OPENCV SCRIPT AUTH GDAL CURL MODEL OMP)
    option(FF_HAS_${flag} "Enable ${flag} support" ON)          # creates a cache entry

    # Emit a compile-time define with 0/1
    #   $<BOOL:…>  →  "1" if TRUE, "0" if FALSE
    add_compile_definitions(FF_HAS_${flag}=$<BOOL:${FF_HAS_${flag}}>)
endforeach()

# ------------------------------------------------------------------------------
# Build-type defaults (Debug = –g –O3)   ---------------------------------------
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Debug/Release/RelWithDebInfo/MinSizeRel" FORCE)
endif()

set(CMAKE_CXX_STANDARD            20)
set(CMAKE_CXX_STANDARD_REQUIRED   ON)
set(CMAKE_CXX_EXTENSIONS          OFF)

set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0 ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
if(FF_BUILD_STATIC)
  string(APPEND CMAKE_CXX_FLAGS_DEBUG   " -static -static-libstdc++ -static-libgcc -fno-PIE -no-pie -Wl,--default-image-base-low -Wl,--dynamicbase=no")
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -static -static-libstdc++ -static-libgcc -Wl,--default-image-base-low")
endif()

# ------------------------------------------------------------------------------
# Platform / tool-chain tweaks
# ------------------------------------------------------------------------------
if(NOT FF_BUILD_STATIC)
    find_package(OpenMP)                             # optional
    if(OpenMP_CXX_FOUND)
      add_compile_definitions(FF_HAS_OPENMP=1)
      set(CMAKE_THREAD_LIBS_INIT    ${OpenMP_CXX_LIBRARIES})
      set(CMAKE_CXX_FLAGS           "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()
# ------------------------------------------------------------------------------
# External dependencies
# ------------------------------------------------------------------------------
list(APPEND FH_LINK_LIBS)        # collect libs that every target needs
if(FF_HAS_VISUAL)
find_package(OpenGL REQUIRED)
add_library(glad STATIC "external/glad/src/glad.c")
target_include_directories(glad PUBLIC include)   # glad’s `glad/glad.h`
add_library(FH::glad ALIAS glad)
list(APPEND FH_LINK_LIBS FH::glad)
endif()
find_package(Threads REQUIRED)
list(APPEND FH_LINK_LIBS Threads::Threads)
if(USE_BUILTIN_EXTERNAL_BINS)
  include(cmake/BuiltinDeps.cmake)               # populates FH_LINK_LIBS
else()
  find_package(OpenCV REQUIRED   COMPONENTS core imgproc      QUIET)
  find_package(GDAL   REQUIRED                                   )
  find_package(CURL   REQUIRED                                   )
  find_package(OpenSSL REQUIRED                                 )

  list(APPEND FH_LINK_LIBS
       OpenCV::opencv_core OpenCV::opencv_imgproc
       GDAL::GDAL CURL::libcurl OpenSSL::SSL)
endif()
message("Print external contents")
print_target_contents(FH_externallib)
set_property(GLOBAL PROPERTY FH_3P_TARGET FH::externallib)

if(FF_HAS_VISUAL)
    message("add glfw dir")
    add_subdirectory(external/glfw)

endif()
# ------------------------------------------------------------------------------
# Sub-directories  -------------------------------------------------------------
# ------------------------------------------------------------------------------
set(FH_SUBDIRS
    external common util math geospatial io algorithms web  model scripting visualization)

foreach(dir IN LISTS FH_SUBDIRS)
  add_subdirectory(${dir})
endforeach()

set(FH_ALL_INCLUDE_DIRS)
foreach(dir IN LISTS FH_SUBDIRS)
    list(APPEND FH_ALL_INCLUDE_DIRS  "${CMAKE_CURRENT_SOURCE_DIR}/${dir}")
endforeach()

foreach(dir IN LISTS FH_SUBDIRS)
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/${dir}")
endforeach()

add_library(FH_includes INTERFACE)

target_include_directories(FH_includes INTERFACE ${FH_ALL_INCLUDE_DIRS})
# Later:   target_link_libraries(<any target>  PUBLIC FH_includes)


# ------------------------------------------------------------------------------
# Link libraries into the two executables in /app
# ------------------------------------------------------------------------------
# Every sub-dir creates a STATIC library named FH::<name> (see template below).
# Collect the ones that exist; this way the CMakeLists in /app does not need
# manual maintenance when you add/remove modules.
get_property(all_targets GLOBAL PROPERTY FH_ALL_TARGETS)
foreach(t IN LISTS all_targets)
  list(APPEND FH_ALL_LIBS FH::${t})
endforeach()

message("${FH_ALL_LIBS}")

add_subdirectory(app)

install(TARGETS fastflood RUNTIME DESTINATION bin) #fh
# ------------------------------------------------------------------------------
# Optional installer with CPack  ----------------------------------------------
# ------------------------------------------------------------------------------
if(ENABLE_INSTALLER)
  include(InstallRequiredSystemLibraries)
  install(TARGETS fastflood
          RUNTIME DESTINATION bin
          COMPONENT runtime) #fh

  # bundle any runtime dlls or extra files you still need 
  if(USE_BUILTIN_EXTERNAL_BINS AND WIN32)
    install(DIRECTORY external/bin/windows/ DESTINATION bin COMPONENT runtime)
  endif()

  set(CPACK_PACKAGE_NAME               "FastHazard")
  set(CPACK_PACKAGE_VENDOR             "Fast Hazard B.V.")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "FastHazard modelling suite")
  set(CPACK_PACKAGE_VERSION            ${FASTHAZARD_VERSION})
  set(CPACK_RESOURCE_FILE_LICENSE      "${CMAKE_SOURCE_DIR}/LICENSE")
  set(CPACK_PACKAGE_EXECUTABLES        fh "FastHazard"  fastflood "FastFlood CLI")
  if(WIN32)
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_NSIS_MODIFY_PATH ON)
  else()
    set(CPACK_GENERATOR "TGZ;ZIP")
  endif()
  include(CPack)
endif()

# ------------------------------------------------------------------------------
# Experimental switch for C++20 modules
# ------------------------------------------------------------------------------
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.28)
  set(CMAKE_EXPERIMENTAL_CXX_MODULE_CMAKE_API "cxx20" CACHE STRING "" FORCE)
endif()




# ------------------------------------------------------------------------------
# Unit testing
# ------------------------------------------------------------------------------
option(FH_BUILD_TESTS "Build tests" ON)
if(FH_BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_subdirectory("test")
endif()

# ------------------------------------------------------------------------------
# Final flags to be set for the default options
# ------------------------------------------------------------------------------
SET(GLFW_BUILD_DOCS OFF)
SET(GLFW_BUILD_TESTS OFF)
SET(GLFW_BUILD_EXAMPLES OFF)
SET(GLFW_INSTALL OFF)
