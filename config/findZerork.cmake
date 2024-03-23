
# check to see if the ZERORK_DIR was specified, if not build
if (NOT (DEFINED ENV{ZERORK_DIR}))
    message(STATUS "ZERORK_DIR not set.  Downloading and building zerork...")

    FetchContent_Declare(zerork
            GIT_REPOSITORY https://github.com/LLNL/zero-rk.git
            GIT_TAG 71befd3f6ff6f49a66e07b82659b88d2bc1beb8c #points to zerork master branch 3/22/24
            )
    FetchContent_MakeAvailable(zerork)

    set_include_directories_as_system(zerork)

    install(TARGETS zerork_cfd_plugin_gpu zerork_cfd_plugin ckconverter zerork_vectormath zerorkutilities zerork zerork_cuda spify
            EXPORT ablateTargets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            )


elseif (DEFINED ENV{ZERORK_DIR})
    message(STATUS "Found ZERORK_DIR, using prebuilt zerork")


    add_library(zerork_cfd_plugin INTERFACE IMPORTED GLOBAL)
    target_include_directories(zerork_cfd_plugin INTERFACE "$ENV{ZERORK_DIR}/include")
    target_link_libraries(zerork_cfd_plugin INTERFACE "$ENV{ZERORK_DIR}/lib/libzerork_cfd_plugin.so")


#    add_library(zerork_cfd_plugin ALIAS zerork)


elseif ()
    find_package(zerork REQUIRED)
endif ()

