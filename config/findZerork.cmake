
# check to see if the ZERORK_DIR was specified, if not build
if (NOT (DEFINED ENV{ZERORK_DIR}))
    message(STATUS "ZERORK_DIR not set.  Downloading and building zerork...")

    # CPU build
    if (NOT (DEFINED ENV{ABLATE_GPU}))
        message(STATUS "Builing zerork for CPUs.")

        FetchContent_Declare(zerork
                GIT_REPOSITORY https://github.com/LLNL/zero-rk.git
                GIT_TAG b032c7d3ce4120aa9192908303c2671a4c0170f1  #git main branch for both cpu and CUDA
        )
        FetchContent_MakeAvailable(zerork)

        set_include_directories_as_system(zerork)

        install(TARGETS zerork_cfd_plugin zerork_vectormath ckconverter zerorkutilities zerork spify
                EXPORT ablateTargets
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )


        # exclude some of the zero rk builds by default that are not needed and causing issues on macOS
        set_target_properties(zerork_flame_api zerork_flame_api_tester.x zerork_flame_api_mpi zerork_flame_api_tester_mpi.x premixed_steady_flame_solver.x premixed_steady_flame_solver_mpi.x PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

        add_library(ZERORK::zerork_cfd_plugin ALIAS zerork_cfd_plugin)

        #Nvidia build
    elseif ($ENV{ABLATE_GPU} STREQUAL "CUDA")
        message(STATUS "Builing zerork for Nvidia GPUs with cuda.")

        FetchContent_Declare(zerork
                GIT_REPOSITORY https://github.com/LLNL/zero-rk.git
                GIT_TAG b032c7d3ce4120aa9192908303c2671a4c0170f1  #git main branch for both cpu and CUDA
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

        add_library(ZERORK::zerork_cfd_plugin ALIAS zerork_cfd_plugin_gpu)

    elseif ($ENV{ABLATE_GPU} STREQUAL "ROCM")
        message(STATUS "Builing zerork for AMD GPUs and hip.")

        FetchContent_Declare(zerork
                GIT_REPOSITORY https://github.com/LLNL/zero-rk.git
                GIT_TAG adb46511232e94925fe2b8ee972e30e3635c425d  #points to zerork hip branch
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

        add_library(ZERORK::zerork_cfd_plugin ALIAS zerork_cfd_plugin_gpu)

    endif ()

elseif (DEFINED ENV{ZERORK_DIR})
    message(STATUS "Found ZERORK_DIR, using prebuilt zerork")

    add_library(zerork_cfd_plugin INTERFACE IMPORTED GLOBAL)
    target_include_directories(zerork_cfd_plugin INTERFACE "$ENV{ZERORK_DIR}/include")
    target_link_libraries(zerork_cfd_plugin INTERFACE "$ENV{ZERORK_DIR}/lib/libzerork_cfd_plugin.so")

    add_library(ZERORK::zerork_cfd_plugin ALIAS zerork_cfd_plugin)

endif ()

