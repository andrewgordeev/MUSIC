# Install script for directory: /home/andrew/Documents/Work/music-hydro/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/andrew/Documents/Work/music-hydro")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/libmusic.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/libmusic.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/libmusic.so"
         RPATH "/home/andrew/Documents/Work/music-hydro")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/andrew/Documents/Work/music-hydro/libmusic.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/andrew/Documents/Work/music-hydro" TYPE SHARED_LIBRARY FILES "/home/andrew/Documents/Work/music-hydro/src/libmusic.so")
  if(EXISTS "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/libmusic.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/libmusic.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/libmusic.so"
         OLD_RPATH ":::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/home/andrew/Documents/Work/music-hydro")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/libmusic.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/MUSIChydro" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/MUSIChydro")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/MUSIChydro"
         RPATH "/home/andrew/Documents/Work/music-hydro")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/andrew/Documents/Work/music-hydro/MUSIChydro")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/andrew/Documents/Work/music-hydro" TYPE EXECUTABLE FILES "/home/andrew/Documents/Work/music-hydro/src/MUSIChydro")
  if(EXISTS "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/MUSIChydro" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/MUSIChydro")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/MUSIChydro"
         OLD_RPATH "/home/andrew/Documents/Work/music-hydro/src:"
         NEW_RPATH "/home/andrew/Documents/Work/music-hydro")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/andrew/Documents/Work/music-hydro/MUSIChydro")
    endif()
  endif()
endif()

