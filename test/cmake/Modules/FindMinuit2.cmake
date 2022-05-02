find_path(MINUIT2_INCLUDE_DIR Minuit2/MnUserFcn.h
          HINTS $ENV{HOME}/local/include $ENV{HOME}/include)

find_library(MINUIT2_LIBRARY NAMES Minuit2
             HINTS $ENV{HOME}/local/lib $ENV{HOME}/lib)

set (MINUIT2_LIBRARIES ${MINUIT2_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Minuit2  DEFAULT_MSG
                                  MINUIT2_LIBRARY MINUIT2_INCLUDE_DIR)

mark_as_advanced(MINUIT2_INCLUDE_DIR MINUIT2_LIBRARY)
