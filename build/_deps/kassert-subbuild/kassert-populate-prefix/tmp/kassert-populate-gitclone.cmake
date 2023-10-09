
if(NOT "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-gitinfo.txt" IS_NEWER_THAN "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"  clone --no-checkout --config "advice.detachedHead=false" "https://github.com/kamping-site/kassert.git" "kassert-src"
    WORKING_DIRECTORY "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/kamping-site/kassert.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  checkout e683aef --
  WORKING_DIRECTORY "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'e683aef'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git"  submodule update --recursive --init 
    WORKING_DIRECTORY "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-src"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-gitinfo.txt"
    "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-gitclone-lastrun.txt'")
endif()

