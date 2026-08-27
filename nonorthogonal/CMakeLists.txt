cmake_minimum_required(VERSION 3.10)
project(run_mhd LANGUAGES C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

# This application subtree intentionally has one frozen production target.
# The standard orthogonal targets remain available from the repository root.
set(COORD_TYPE "3" CACHE STRING "Frozen non-orthogonal coordinate backend")
set(PROBLEM "22" CACHE STRING "Frozen Earth M-I-coupled problem")
if(NOT COORD_TYPE STREQUAL "3" OR NOT PROBLEM STREQUAL "22")
    message(FATAL_ERROR
            "nonorthogonal/ supports only the frozen COORD_TYPE=3, PROBLEM=22 production target")
endif()

# Find MPI and require at least MPI 1.2
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

option(GAMERA_STRICT_INNER_MAGNETIC_WALL
       "Copy first-active-shell magnetic flux and Cartesian B into every inner ghost"
       ON)
option(GAMERA_STRICT_INNER_GAS_WALL
       "Store the prescribed wall drift velocity directly in every inner gas ghost"
       ON)
option(GAMERA_YINYANG_MFE_INTERFACE
       "Use the experimental MFE-style ghost-only Yin-Yang MHD interface"
       OFF)
option(GAMERA_YINYANG_HDIV_RECONCILE
       "Reconcile active Yin-Yang magnetic face fluxes with a discrete H(div) projection"
       OFF)
option(GAMERA_YINYANG_HDIV_OPTIMIZED
       "Use patch-local communication, cached metric diagonals, and OpenMP kernels in the Yin-Yang H(div) projection"
       OFF)
option(GAMERA_YINYANG_HDIV_DISTRIBUTED
       "Distribute Yin-Yang H(div) faces, vectors, and PCG work over the patch Cartesian communicator"
       OFF)
option(GAMERA_YINYANG_HDIV_PROFILE
       "Log per-call phase timings for the Yin-Yang H(div) projection"
       OFF)
option(GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
       "Run expensive distributed-vs-global H(div) ownership/operator oracles"
       OFF)
option(GAMERA_YINYANG_SPARSE_OVERSET
       "Exchange only requested Yin-Yang donor cells over MPI distributed graphs"
       OFF)
option(GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
       "Emit structured sparse Yin-Yang overset timing and payload records"
       OFF)
option(GAMERA_BENCH_PROFILE_JSONL
       "Buffer per-step benchmark timings and write one JSONL file per MPI rank after the timed loop"
       OFF)
option(GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS
       "Build the standalone single-species H+ bounce-averaged inner-magnetosphere core and tests"
       OFF)
option(GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
       "Enable the experimental default-OFF component-matched H+ feedback driver gate"
       OFF)
option(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
       "Enable the experimental default-OFF online H+ diagnostic callback/scheduler seam (no MHD feedback)"
       OFF)
option(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK
       "Enable the experimental default-OFF live online-H+ pressure-only MHD feedback seam"
       OFF)
# This subtree is the frozen global-magnetosphere production release. The
# unfinished H+ / inner-magnetosphere research branches are intentionally not
# distributed here, so fail clearly instead of reaching missing source files.
if(GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS OR
   GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK OR
   GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC OR
   GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK)
    message(FATAL_ERROR
            "Experimental H+ / inner-magnetosphere features are not part of the frozen production release")
endif()
if(PROBLEM STREQUAL "19")
    set(_GAMERA_DEFAULT_INNER_RADIUS_RE "4.0")
    set(_GAMERA_DEFAULT_OUTER_RADIUS_RE "100.0")
else()
    set(_GAMERA_DEFAULT_INNER_RADIUS_RE "2.5")
    set(_GAMERA_DEFAULT_OUTER_RADIUS_RE "200.0")
endif()
set(GAMERA_BOW_SHOCK_INNER_RADIUS_RE "${_GAMERA_DEFAULT_INNER_RADIUS_RE}" CACHE STRING
    "Inner radial boundary in RE for stretched Yin-Yang bow-shock/Earth problems")
set(GAMERA_BOW_SHOCK_OUTER_RADIUS_RE "${_GAMERA_DEFAULT_OUTER_RADIUS_RE}" CACHE STRING
    "Outer radial boundary in RE for stretched Yin-Yang bow-shock/Earth problems")
set(GAMERA_BOW_SHOCK_RADIAL_STRETCH "3.0" CACHE STRING
    "Exponential radial-map stretch for stretched Yin-Yang bow-shock/Earth problems")
set(GAMERA_EARTH_RADIAL_MAP_VERSION "4" CACHE STRING
    "Earth magnetosphere radial map: 1=legacy exponential, 2=archived 4--100 RE, 3=archived 3--15--50--200 RE, 4=2.5--14--200 RE C1 candidate production profile")
if(NOT GAMERA_EARTH_RADIAL_MAP_VERSION STREQUAL "1" AND
   NOT GAMERA_EARTH_RADIAL_MAP_VERSION STREQUAL "2" AND
   NOT GAMERA_EARTH_RADIAL_MAP_VERSION STREQUAL "3" AND
   NOT GAMERA_EARTH_RADIAL_MAP_VERSION STREQUAL "4")
    message(FATAL_ERROR "GAMERA_EARTH_RADIAL_MAP_VERSION must be 1, 2, 3, or 4")
endif()
if(GAMERA_YINYANG_MFE_INTERFACE AND GAMERA_YINYANG_HDIV_RECONCILE)
    message(FATAL_ERROR
            "GAMERA_YINYANG_MFE_INTERFACE and GAMERA_YINYANG_HDIV_RECONCILE are mutually exclusive interface algorithms")
endif()
if(GAMERA_YINYANG_HDIV_DISTRIBUTED AND
   (NOT GAMERA_YINYANG_HDIV_RECONCILE OR
    NOT GAMERA_YINYANG_HDIV_OPTIMIZED))
    message(FATAL_ERROR
            "GAMERA_YINYANG_HDIV_DISTRIBUTED requires GAMERA_YINYANG_HDIV_RECONCILE=ON and GAMERA_YINYANG_HDIV_OPTIMIZED=ON")
endif()
if(GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY AND
   NOT GAMERA_YINYANG_HDIV_DISTRIBUTED)
    message(FATAL_ERROR
            "GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY requires GAMERA_YINYANG_HDIV_DISTRIBUTED=ON")
endif()
if(GAMERA_YINYANG_SPARSE_OVERSET_PROFILE AND
   NOT GAMERA_YINYANG_SPARSE_OVERSET)
    message(FATAL_ERROR
            "GAMERA_YINYANG_SPARSE_OVERSET_PROFILE requires GAMERA_YINYANG_SPARSE_OVERSET=ON")
endif()
if(GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK AND
   NOT GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS)
    message(FATAL_ERROR
            "GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK requires GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS=ON")
endif()
if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC AND
   NOT GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS)
    message(FATAL_ERROR
            "GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC requires GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS=ON")
endif()
if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC AND
   (NOT COORD_TYPE STREQUAL "3" OR NOT PROBLEM STREQUAL "22"))
    message(FATAL_ERROR
            "GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC requires COORD_TYPE=3 and M-I-coupled PROBLEM=22")
endif()
if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK AND
   NOT GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC)
    message(FATAL_ERROR
            "GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK requires GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC=ON")
endif()
if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK AND
   NOT GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS)
    message(FATAL_ERROR
            "GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK requires GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS=ON")
endif()
if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK AND
   GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK)
    message(FATAL_ERROR
            "Live online-H+ pressure feedback and prescribed matched-H+ feedback are mutually exclusive")
endif()
if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK AND
   (NOT COORD_TYPE STREQUAL "3" OR NOT PROBLEM STREQUAL "22"))
    message(FATAL_ERROR
            "GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK requires COORD_TYPE=3 and M-I-coupled PROBLEM=22")
endif()

# Find MPI and require at least MPI 1.2
find_package(MPI REQUIRED COMPONENTS C)  # Cxx is for C++, omit if not needed
find_package(HDF5 REQUIRED COMPONENTS C)

# OpenMP configuration
if(APPLE)
    find_library(OMP_LIB omp PATHS /opt/homebrew/opt/libomp/lib REQUIRED)
    find_path(OMP_INCLUDE_DIR omp.h PATHS /opt/homebrew/opt/libomp/include REQUIRED)

    message(STATUS "Found OpenMP library: ${OMP_LIB}")
    message(STATUS "Found OpenMP include dir: ${OMP_INCLUDE_DIR}")

    add_library(OpenMP::OpenMP_C UNKNOWN IMPORTED)
    set_target_properties(OpenMP::OpenMP_C PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${OMP_INCLUDE_DIR}"
            IMPORTED_LOCATION "${OMP_LIB}"
    )
    set(OpenMP_C_FOUND TRUE)
else()
    find_package(OpenMP REQUIRED)
endif()


if(MPI_C_FOUND AND MPI_C_COMPILER)
    message(STATUS "MPI C compiler found: ${MPI_C_COMPILER}")
    message(STATUS "MPI C include path: ${MPI_C_INCLUDE_PATH}")
    message(STATUS "MPI C libraries: ${MPI_C_LIBRARIES}")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")

    if(CMAKE_BUILD_TYPE MATCHES Release)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fomit-frame-pointer")

        find_package(OpenMP)
        if(OpenMP_C_FOUND)
            message(STATUS "OpenMP found! Version: ${OpenMP_C_VERSION}")
            message(STATUS "OpenMP flags: ${OpenMP_C_FLAGS}")
            message(STATUS "OpenMP include dirs: ${OpenMP_C_INCLUDE_DIRS}")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        endif()


    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g") # Debug模式
    endif()
    # Enable testing
    enable_testing()

    include_directories(${CMAKE_SOURCE_DIR}/logc)
    include_directories(${CMAKE_SOURCE_DIR}/src)
    file(GLOB COMMON_SRC_FILES "src/*.c" "logc/*.c" "math_lib/*.c")

    # Keep the reusable non-orthogonal GAMERA kernels in a standalone library.
    # COORD_TYPE=3 links this library to the legacy MPI/OpenMP driver, while the
    # unit tests exercise the numerical layer without HDF5 or a full run.
    add_library(gamera_nonorthogonal_kernels STATIC
            src/nonorthogonal/nonorthogonal_absolute_cadence.c
            src/nonorthogonal/nonorthogonal_advance.c
            src/nonorthogonal/nonorthogonal_background.c
            src/nonorthogonal/nonorthogonal_flux.c
            src/nonorthogonal/nonorthogonal_geometry.c
            src/nonorthogonal/nonorthogonal_grid.c
            src/nonorthogonal/nonorthogonal_initialization.c
            src/nonorthogonal/nonorthogonal_legacy_adapter.c
            src/nonorthogonal/nonorthogonal_mesh.c
            src/nonorthogonal/nonorthogonal_mi.c
            src/nonorthogonal/nonorthogonal_mi_dpb.c
            src/nonorthogonal/nonorthogonal_electron_precipitation.c
            src/nonorthogonal/nonorthogonal_operators.c
            src/nonorthogonal/nonorthogonal_reconstruction.c
            src/nonorthogonal/nonorthogonal_state.c
            src/nonorthogonal/nonorthogonal_step.c
            src/nonorthogonal/nonorthogonal_storage.c
            src/nonorthogonal/nonorthogonal_sweep.c
            src/nonorthogonal/nonorthogonal_yinyang.c)
    target_include_directories(gamera_nonorthogonal_kernels PUBLIC
            ${CMAKE_SOURCE_DIR}/src/nonorthogonal)
    target_include_directories(gamera_nonorthogonal_kernels PRIVATE
            ${CMAKE_SOURCE_DIR}/src)
    # These kernels deliberately reject non-finite geometry/state inputs.
    # Preserve those checks even when the legacy Release build enables
    # -ffast-math globally.
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(gamera_nonorthogonal_kernels PRIVATE
                -fno-fast-math)
    endif()
    if(NOT MSVC)
        target_link_libraries(gamera_nonorthogonal_kernels PUBLIC m)
    endif()
    if(APPLE)
        target_compile_options(gamera_nonorthogonal_kernels PRIVATE
                "-Xpreprocessor" "-fopenmp")
        target_include_directories(gamera_nonorthogonal_kernels PRIVATE
                ${OMP_INCLUDE_DIR})
        target_link_libraries(gamera_nonorthogonal_kernels PUBLIC ${OMP_LIB})
    else()
        target_link_libraries(gamera_nonorthogonal_kernels PUBLIC
                OpenMP::OpenMP_C)
    endif()

    if(GAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS)
        add_library(gamera_inner_magnetosphere STATIC
                src/inner_magnetosphere/gamera_inner_magnetosphere.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_coupling.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_dipole.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_injection.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_online_hplus_coordinates.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_online_hplus_initialization.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_losses.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_live_tracer.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_mapped_shell.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_one_way.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_online_hplus.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_pressure_mapping.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_pressure_feedback.c
                src/inner_magnetosphere/gamera_inner_magnetosphere_target_provider.c)
        target_include_directories(gamera_inner_magnetosphere PUBLIC
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(gamera_inner_magnetosphere PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
            set_property(SOURCE
                    src/inner_magnetosphere/gamera_inner_magnetosphere_online_hplus_initialization.c
                    APPEND PROPERTY COMPILE_OPTIONS -std=c99)
        endif()
        if(NOT MSVC)
            target_link_libraries(gamera_inner_magnetosphere PUBLIC m)
        endif()
        add_library(gamera_nonorthogonal_inner_magnetosphere STATIC
                src/nonorthogonal/nonorthogonal_inner_magnetosphere_feedback.c
                src/nonorthogonal/nonorthogonal_inner_magnetosphere_online_scheduler.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_feedback.c)
        target_include_directories(gamera_nonorthogonal_inner_magnetosphere PUBLIC
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
        target_link_libraries(gamera_nonorthogonal_inner_magnetosphere PUBLIC
                gamera_nonorthogonal_kernels
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(gamera_nonorthogonal_inner_magnetosphere PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_executable(test_inner_magnetosphere_hplus
                tests/test_inner_magnetosphere_hplus.c)
        target_link_libraries(test_inner_magnetosphere_hplus PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_hplus PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_hplus
                COMMAND test_inner_magnetosphere_hplus)
        add_executable(test_inner_magnetosphere_target_provider
                tests/test_inner_magnetosphere_target_provider.c)
        target_link_libraries(test_inner_magnetosphere_target_provider PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_target_provider PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_target_provider
                COMMAND test_inner_magnetosphere_target_provider)
        add_executable(test_inner_magnetosphere_online_hplus
                tests/test_inner_magnetosphere_online_hplus.c)
        target_link_libraries(test_inner_magnetosphere_online_hplus PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_online_hplus PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_online_hplus
                COMMAND test_inner_magnetosphere_online_hplus)
        add_executable(test_inner_magnetosphere_dipole
                tests/test_inner_magnetosphere_dipole.c)
        target_link_libraries(test_inner_magnetosphere_dipole PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_dipole PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_dipole
                COMMAND test_inner_magnetosphere_dipole)
        add_executable(test_inner_magnetosphere_live_tracer
                tests/test_inner_magnetosphere_live_tracer.c)
        target_link_libraries(test_inner_magnetosphere_live_tracer PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_live_tracer PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_live_tracer
                COMMAND test_inner_magnetosphere_live_tracer)
        add_library(gamera_inner_magnetosphere_live_mpi STATIC
                src/inner_magnetosphere/gamera_inner_magnetosphere_live_mpi_producer.c)
        target_include_directories(gamera_inner_magnetosphere_live_mpi PUBLIC
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
        target_link_libraries(gamera_inner_magnetosphere_live_mpi PUBLIC
                gamera_inner_magnetosphere MPI::MPI_C)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(gamera_inner_magnetosphere_live_mpi PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_executable(test_inner_magnetosphere_live_mpi_producer
                tests/test_inner_magnetosphere_live_mpi_producer.c)
        target_link_libraries(
                test_inner_magnetosphere_live_mpi_producer PRIVATE
                gamera_inner_magnetosphere_live_mpi)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_inner_magnetosphere_live_mpi_producer PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_live_mpi_producer_2rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_inner_magnetosphere_live_mpi_producer>
                        ${MPIEXEC_POSTFLAGS})
        add_test(NAME inner_magnetosphere_live_mpi_producer_4rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_inner_magnetosphere_live_mpi_producer>
                        ${MPIEXEC_POSTFLAGS})
        add_test(NAME inner_magnetosphere_live_mpi_producer_8rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_inner_magnetosphere_live_mpi_producer>
                        ${MPIEXEC_POSTFLAGS})
        set_tests_properties(
                inner_magnetosphere_live_mpi_producer_2rank
                inner_magnetosphere_live_mpi_producer_4rank
                inner_magnetosphere_live_mpi_producer_8rank
                PROPERTIES TIMEOUT 60)
        add_executable(test_nonorthogonal_live_driver_default_off
                tests/test_nonorthogonal_live_driver_default_off.c)
        target_include_directories(
                test_nonorthogonal_live_driver_default_off PRIVATE
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_nonorthogonal_live_driver_default_off PRIVATE
                    -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_live_driver_default_off
                COMMAND test_nonorthogonal_live_driver_default_off)
        add_executable(test_inner_magnetosphere_one_way
                tests/test_inner_magnetosphere_one_way.c)
        target_link_libraries(test_inner_magnetosphere_one_way PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_one_way PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_one_way
                COMMAND test_inner_magnetosphere_one_way)
        add_executable(test_inner_magnetosphere_injection
                tests/test_inner_magnetosphere_injection.c)
        target_link_libraries(test_inner_magnetosphere_injection PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_injection PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_injection
                COMMAND test_inner_magnetosphere_injection)
        add_executable(test_inner_magnetosphere_pressure_feedback
                tests/test_inner_magnetosphere_pressure_feedback.c)
        target_link_libraries(test_inner_magnetosphere_pressure_feedback PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_pressure_feedback PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_pressure_feedback
                COMMAND test_inner_magnetosphere_pressure_feedback)
        add_executable(test_inner_magnetosphere_pressure_mapping
                tests/test_inner_magnetosphere_pressure_mapping.c)
        target_link_libraries(test_inner_magnetosphere_pressure_mapping PRIVATE
                gamera_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_pressure_mapping PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_pressure_mapping
                COMMAND test_inner_magnetosphere_pressure_mapping)
        add_executable(test_nonorthogonal_inner_magnetosphere_feedback
                tests/test_nonorthogonal_inner_magnetosphere_feedback.c)
        target_link_libraries(test_nonorthogonal_inner_magnetosphere_feedback PRIVATE
                gamera_nonorthogonal_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_inner_magnetosphere_feedback PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_inner_magnetosphere_feedback
                COMMAND test_nonorthogonal_inner_magnetosphere_feedback)
        add_executable(test_nonorthogonal_online_feedback_controller
                tests/test_nonorthogonal_online_feedback_controller.c)
        target_link_libraries(test_nonorthogonal_online_feedback_controller PRIVATE
                gamera_nonorthogonal_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_online_feedback_controller PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_online_feedback_controller
                COMMAND test_nonorthogonal_online_feedback_controller)
        add_executable(test_nonorthogonal_inner_magnetosphere_online_scheduler
                tests/test_nonorthogonal_inner_magnetosphere_online_scheduler.c)
        target_link_libraries(
                test_nonorthogonal_inner_magnetosphere_online_scheduler PRIVATE
                gamera_nonorthogonal_inner_magnetosphere)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_nonorthogonal_inner_magnetosphere_online_scheduler
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_inner_magnetosphere_online_scheduler
                COMMAND test_nonorthogonal_inner_magnetosphere_online_scheduler)
        add_executable(test_nonorthogonal_absolute_cadence
                tests/test_nonorthogonal_absolute_cadence.c)
        target_link_libraries(test_nonorthogonal_absolute_cadence PRIVATE
                gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_absolute_cadence PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_absolute_cadence
                COMMAND test_nonorthogonal_absolute_cadence)
        find_package(HDF5 REQUIRED COMPONENTS C)
        add_library(gamera_inner_magnetosphere_restart_hdf5 STATIC
                src/inner_magnetosphere/gamera_inner_magnetosphere_online_hplus_restart_hdf5.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_include_directories(gamera_inner_magnetosphere_restart_hdf5
                PUBLIC ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere
                PRIVATE ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                        ${HDF5_INCLUDE_DIRS})
        target_link_libraries(gamera_inner_magnetosphere_restart_hdf5
                PUBLIC gamera_inner_magnetosphere
                PRIVATE ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(gamera_inner_magnetosphere_restart_hdf5
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_executable(test_inner_magnetosphere_online_hplus_restart_hdf5
                tests/test_inner_magnetosphere_online_hplus_restart_hdf5.c)
        target_include_directories(
                test_inner_magnetosphere_online_hplus_restart_hdf5
                PRIVATE ${HDF5_INCLUDE_DIRS})
        target_link_libraries(
                test_inner_magnetosphere_online_hplus_restart_hdf5
                PRIVATE gamera_inner_magnetosphere_restart_hdf5
                        ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_inner_magnetosphere_online_hplus_restart_hdf5
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_online_hplus_restart_hdf5
                COMMAND test_inner_magnetosphere_online_hplus_restart_hdf5)
        if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK)
            # Strict rank-local pressure-feedback restart payload.  Keep the
            # HDF5/MPI dependency and its tests entirely out of feedback-OFF
            # builds; the production executable links this library below only
            # when the compile-time feedback gate is enabled.
            add_library(gamera_nonorthogonal_online_feedback_restart STATIC
                    src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5.c
                    src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_feedback_restart_mpi.c
                    src/nonorthogonal_backend/nonorthogonal_sha256.c)
            target_include_directories(
                    gamera_nonorthogonal_online_feedback_restart PUBLIC
                    ${CMAKE_SOURCE_DIR}/src
                    ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                    ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                    ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere
                    ${HDF5_INCLUDE_DIRS})
            target_link_libraries(
                    gamera_nonorthogonal_online_feedback_restart PUBLIC
                    gamera_nonorthogonal_inner_magnetosphere
                    MPI::MPI_C ${HDF5_LIBRARIES} m)
            if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
                target_compile_options(
                        gamera_nonorthogonal_online_feedback_restart PRIVATE
                        -fno-fast-math -Wall -Wextra -Werror -pedantic)
            endif()

            add_executable(
                    test_nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5
                    tests/test_nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5.c)
            target_link_libraries(
                    test_nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5
                    PRIVATE gamera_nonorthogonal_online_feedback_restart)
            if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
                target_compile_options(
                        test_nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5
                        PRIVATE -fno-fast-math -Wall -Wextra -Werror
                                -pedantic)
            endif()
            add_test(
                    NAME nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5
                    COMMAND test_nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5)

            add_executable(
                    test_nonorthogonal_inner_magnetosphere_online_feedback_restart_mpi
                    tests/test_nonorthogonal_inner_magnetosphere_online_feedback_restart_mpi.c)
            target_link_libraries(
                    test_nonorthogonal_inner_magnetosphere_online_feedback_restart_mpi
                    PRIVATE gamera_nonorthogonal_online_feedback_restart)
            if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
                target_compile_options(
                        test_nonorthogonal_inner_magnetosphere_online_feedback_restart_mpi
                        PRIVATE -fno-fast-math -Wall -Wextra -Werror
                                -pedantic)
            endif()
            add_test(
                    NAME nonorthogonal_inner_magnetosphere_online_feedback_restart_2rank
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                            ${MPIEXEC_PREFLAGS}
                            $<TARGET_FILE:test_nonorthogonal_inner_magnetosphere_online_feedback_restart_mpi>
                            ${MPIEXEC_POSTFLAGS})
            set_tests_properties(
                    nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5
                    nonorthogonal_inner_magnetosphere_online_feedback_restart_2rank
                    PROPERTIES TIMEOUT 60)
        endif()
        add_library(gamera_inner_magnetosphere_snapshot_hdf5 STATIC
                src/inner_magnetosphere/gamera_inner_magnetosphere_online_hplus_snapshot_hdf5.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_include_directories(gamera_inner_magnetosphere_snapshot_hdf5
                PUBLIC ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere
                PRIVATE ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                        ${HDF5_INCLUDE_DIRS})
        target_link_libraries(gamera_inner_magnetosphere_snapshot_hdf5
                PUBLIC gamera_inner_magnetosphere
                PRIVATE ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(gamera_inner_magnetosphere_snapshot_hdf5
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_executable(test_inner_magnetosphere_online_hplus_snapshot_hdf5
                tests/test_inner_magnetosphere_online_hplus_snapshot_hdf5.c)
        target_include_directories(
                test_inner_magnetosphere_online_hplus_snapshot_hdf5 PRIVATE
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(
                test_inner_magnetosphere_online_hplus_snapshot_hdf5 PRIVATE
                gamera_inner_magnetosphere_snapshot_hdf5 ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_inner_magnetosphere_online_hplus_snapshot_hdf5
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_online_hplus_snapshot_hdf5
                COMMAND test_inner_magnetosphere_online_hplus_snapshot_hdf5)
        add_executable(test_inner_magnetosphere_online_hplus_replay
                tests/test_inner_magnetosphere_online_hplus_replay.c
                tools/inner_magnetosphere_mapped_shell_hdf5.c)
        target_include_directories(test_inner_magnetosphere_online_hplus_replay
                PRIVATE ${CMAKE_SOURCE_DIR}/tools ${HDF5_INCLUDE_DIRS})
        target_link_libraries(test_inner_magnetosphere_online_hplus_replay
                PRIVATE gamera_inner_magnetosphere ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_inner_magnetosphere_online_hplus_replay
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_online_hplus_replay
                COMMAND test_inner_magnetosphere_online_hplus_replay
                ${CMAKE_SOURCE_DIR}/replay/physical_boundary_64_6h_closed/maps
                ${CMAKE_SOURCE_DIR}/replay/physical_boundary_64_6h_closed/component_matched_v1/cadence_states)
        add_executable(test_inner_magnetosphere_online_hplus_initialization
                tests/test_inner_magnetosphere_online_hplus_initialization.c
                tools/inner_magnetosphere_mapped_shell_hdf5.c)
        target_include_directories(
                test_inner_magnetosphere_online_hplus_initialization PRIVATE
                ${CMAKE_SOURCE_DIR}/tools ${HDF5_INCLUDE_DIRS})
        target_link_libraries(
                test_inner_magnetosphere_online_hplus_initialization PRIVATE
                gamera_inner_magnetosphere ${HDF5_LIBRARIES})
        set_target_properties(
                test_inner_magnetosphere_online_hplus_initialization
                PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED YES
                           C_EXTENSIONS NO)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_inner_magnetosphere_online_hplus_initialization PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME inner_magnetosphere_online_hplus_initialization
                COMMAND test_inner_magnetosphere_online_hplus_initialization
                ${CMAKE_SOURCE_DIR}/replay/physical_boundary_64_6h_closed/maps/map_000000.h5
                ${CMAKE_SOURCE_DIR}/replay/physical_boundary_64_6h_closed/component_matched_v1/cadence_states/state_000000.h5)
        add_executable(replay_inner_magnetosphere_hplus
                tools/replay_inner_magnetosphere_hplus.c
                tools/inner_magnetosphere_mapped_shell_hdf5.c)
        target_include_directories(replay_inner_magnetosphere_hplus PRIVATE
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(replay_inner_magnetosphere_hplus PRIVATE
                gamera_inner_magnetosphere ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(replay_inner_magnetosphere_hplus PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_executable(replay_inner_magnetosphere_one_way
                tools/replay_inner_magnetosphere_one_way.c
                tools/inner_magnetosphere_mapped_shell_hdf5.c)
        target_include_directories(replay_inner_magnetosphere_one_way PRIVATE
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(replay_inner_magnetosphere_one_way PRIVATE
                gamera_inner_magnetosphere ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(replay_inner_magnetosphere_one_way PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_executable(replay_inner_magnetosphere_injection
                tools/replay_inner_magnetosphere_injection.c
                tools/inner_magnetosphere_mapped_shell_hdf5.c)
        target_include_directories(replay_inner_magnetosphere_injection PRIVATE
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(replay_inner_magnetosphere_injection PRIVATE
                gamera_inner_magnetosphere ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(replay_inner_magnetosphere_injection PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_executable(replay_inner_magnetosphere_pressure_feedback
                tools/replay_inner_magnetosphere_pressure_feedback.c)
        target_include_directories(replay_inner_magnetosphere_pressure_feedback PRIVATE
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(replay_inner_magnetosphere_pressure_feedback PRIVATE
                gamera_inner_magnetosphere ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(replay_inner_magnetosphere_pressure_feedback PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_executable(replay_inner_magnetosphere_matched_feedback
                tools/replay_inner_magnetosphere_matched_feedback.c
                tools/inner_magnetosphere_mapped_shell_hdf5.c)
        target_include_directories(replay_inner_magnetosphere_matched_feedback PRIVATE
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(replay_inner_magnetosphere_matched_feedback PRIVATE
                gamera_inner_magnetosphere ${HDF5_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(replay_inner_magnetosphere_matched_feedback PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
    endif()

    option(GAMERA_BUILD_NONORTHOGONAL_TESTS
            "Build unit tests for the non-orthogonal GAMERA kernels" ON)
    if(GAMERA_BUILD_NONORTHOGONAL_TESTS)
        add_executable(test_nonorthogonal_kernels
                tests/test_nonorthogonal_kernels.c)
        target_link_libraries(test_nonorthogonal_kernels PRIVATE
                gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_kernels PRIVATE
                    -fno-fast-math)
        endif()
        add_test(NAME nonorthogonal_kernels COMMAND test_nonorthogonal_kernels)

        add_executable(test_nonorthogonal_legacy_adapter
                tests/test_nonorthogonal_legacy_adapter.c
                src/array.c
                logc/log.c)
        target_include_directories(test_nonorthogonal_legacy_adapter PRIVATE
                ${CMAKE_SOURCE_DIR}/src
                ${MPI_C_INCLUDE_PATH})
        target_link_libraries(test_nonorthogonal_legacy_adapter PRIVATE
                gamera_nonorthogonal_kernels
                ${MPI_C_LIBRARIES})
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_legacy_adapter PRIVATE
                    -fno-fast-math)
        endif()
        add_test(NAME nonorthogonal_legacy_adapter
                COMMAND test_nonorthogonal_legacy_adapter)

        add_executable(test_nonorthogonal_problems
                tests/test_nonorthogonal_problems.c)
        target_link_libraries(test_nonorthogonal_problems PRIVATE
                gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_problems PRIVATE
                    -fno-fast-math)
        endif()
        add_test(NAME nonorthogonal_problems COMMAND test_nonorthogonal_problems)

        add_executable(test_nonorthogonal_radial_map
                tests/test_nonorthogonal_radial_map.c
                src/nonorthogonal_radial_map.c)
        target_include_directories(test_nonorthogonal_radial_map PRIVATE
                ${CMAKE_SOURCE_DIR}/src)
        if(NOT MSVC)
            target_link_libraries(test_nonorthogonal_radial_map PRIVATE m)
        endif()
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_radial_map PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_radial_map
                COMMAND test_nonorthogonal_radial_map)

        add_executable(test_solar_wind
                tests/test_solar_wind.c
                src/solar_wind.c)
        target_include_directories(test_solar_wind PRIVATE
                ${CMAKE_SOURCE_DIR}/src)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_solar_wind PRIVATE -fno-fast-math)
        endif()
        if(NOT MSVC)
            target_link_libraries(test_solar_wind PRIVATE m)
        endif()
        add_test(NAME solar_wind COMMAND test_solar_wind)

        add_executable(test_nonorthogonal_background
                tests/test_nonorthogonal_background.c)
        target_link_libraries(test_nonorthogonal_background PRIVATE
                gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_background PRIVATE
                    -fno-fast-math)
        endif()
        add_test(NAME nonorthogonal_background
                COMMAND test_nonorthogonal_background)

        add_executable(test_nonorthogonal_mi
                tests/test_nonorthogonal_mi.c)
        target_link_libraries(test_nonorthogonal_mi PRIVATE
                gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_mi PRIVATE
                    -fno-fast-math)
        endif()
        add_test(NAME nonorthogonal_mi COMMAND test_nonorthogonal_mi)

        add_executable(test_nonorthogonal_electron_precipitation
                tests/test_nonorthogonal_electron_precipitation.c)
        target_link_libraries(test_nonorthogonal_electron_precipitation PRIVATE
                gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_electron_precipitation
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_electron_precipitation
                COMMAND test_nonorthogonal_electron_precipitation)

        add_executable(test_nonorthogonal_mi_dpb
                tests/test_nonorthogonal_mi_dpb.c)
        target_link_libraries(test_nonorthogonal_mi_dpb PRIVATE
                gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_mi_dpb PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_mi_dpb
                COMMAND test_nonorthogonal_mi_dpb)

        add_executable(test_nonorthogonal_mi_snapshot
                tests/test_nonorthogonal_mi_snapshot.c
                src/nonorthogonal_backend/nonorthogonal_mi_coupling.c)
        target_include_directories(test_nonorthogonal_mi_snapshot PRIVATE
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend)
        target_compile_definitions(test_nonorthogonal_mi_snapshot PRIVATE
                GAMERA_MI_SNAPSHOT_TESTING)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_mi_snapshot PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        if(NOT MSVC)
            target_link_libraries(test_nonorthogonal_mi_snapshot PRIVATE m)
        endif()
        add_test(NAME nonorthogonal_mi_snapshot
                COMMAND test_nonorthogonal_mi_snapshot)

        add_executable(test_nonorthogonal_mi_collective_mpi
                tests/test_nonorthogonal_mi_collective_mpi.c
                src/nonorthogonal_backend/nonorthogonal_mi_collective.c)
        target_include_directories(test_nonorthogonal_mi_collective_mpi PRIVATE
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend)
        target_link_libraries(test_nonorthogonal_mi_collective_mpi PRIVATE
                MPI::MPI_C m)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_mi_collective_mpi PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_mi_collective_2rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_nonorthogonal_mi_collective_mpi>
                        ${MPIEXEC_POSTFLAGS})
        set_tests_properties(nonorthogonal_mi_collective_2rank
                PROPERTIES TIMEOUT 15)

        add_executable(test_nonorthogonal_absolute_cadence
                tests/test_nonorthogonal_absolute_cadence.c)
        target_link_libraries(test_nonorthogonal_absolute_cadence PRIVATE
                gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_absolute_cadence PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_absolute_cadence
                COMMAND test_nonorthogonal_absolute_cadence)

        add_executable(test_nonorthogonal_sha256
                tests/test_nonorthogonal_sha256.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_include_directories(test_nonorthogonal_sha256 PRIVATE
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_sha256 PRIVATE
                    -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_sha256 COMMAND test_nonorthogonal_sha256)

        add_executable(test_nonorthogonal_mi_restart_mpi
                tests/test_nonorthogonal_mi_restart_mpi.c
                src/nonorthogonal_backend/nonorthogonal_mi_restart.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_include_directories(test_nonorthogonal_mi_restart_mpi PRIVATE
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(test_nonorthogonal_mi_restart_mpi PRIVATE
                MPI::MPI_C ${HDF5_LIBRARIES} m)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_mi_restart_mpi PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_mi_restart_2rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_nonorthogonal_mi_restart_mpi>
                        ${CMAKE_CURRENT_BINARY_DIR}
                        ${MPIEXEC_POSTFLAGS})
        set_tests_properties(nonorthogonal_mi_restart_2rank
                PROPERTIES TIMEOUT 30)

        if(GAMERA_YINYANG_SPARSE_OVERSET)
            add_executable(test_yinyang_sparse_plan_mpi
                    tests/test_yinyang_sparse_plan_mpi.c
                    src/nonorthogonal_backend/nonorthogonal_yinyang_sparse_plan.c)
            target_include_directories(test_yinyang_sparse_plan_mpi PRIVATE
                    ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend)
            target_compile_definitions(test_yinyang_sparse_plan_mpi PRIVATE
                    GAMERA_YINYANG_SPARSE_OVERSET)
            target_link_libraries(test_yinyang_sparse_plan_mpi PRIVATE
                    MPI::MPI_C)
            add_test(NAME yinyang_sparse_plan_mpi
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                            ${MPIEXEC_PREFLAGS}
                            $<TARGET_FILE:test_yinyang_sparse_plan_mpi>
                            ${MPIEXEC_POSTFLAGS})
            add_test(NAME yinyang_sparse_plan_invalid_exchange_aborts
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                            ${MPIEXEC_PREFLAGS}
                            $<TARGET_FILE:test_yinyang_sparse_plan_mpi>
                            --abort-invalid-exchange ${MPIEXEC_POSTFLAGS})
            add_test(NAME yinyang_sparse_plan_null_build_aborts
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                            ${MPIEXEC_PREFLAGS}
                            $<TARGET_FILE:test_yinyang_sparse_plan_mpi>
                            --abort-null-build ${MPIEXEC_POSTFLAGS})
            add_test(NAME yinyang_sparse_plan_reused_build_aborts
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                            ${MPIEXEC_PREFLAGS}
                            $<TARGET_FILE:test_yinyang_sparse_plan_mpi>
                            --abort-reused-build ${MPIEXEC_POSTFLAGS})
            set_tests_properties(
                    yinyang_sparse_plan_invalid_exchange_aborts
                    yinyang_sparse_plan_null_build_aborts
                    yinyang_sparse_plan_reused_build_aborts
                    PROPERTIES WILL_FAIL TRUE TIMEOUT 10)
        endif()
    endif()

    # Select implementation files based on coordinate system type
    if(COORD_TYPE STREQUAL "0")
        # Cartesian coordinate system
        message(STATUS "Using Cartesian coordinate system")
        file(GLOB CARTESIAN_SRC "src/cartesian/*.c")
        set(SRC_FILES ${COMMON_SRC_FILES} ${CARTESIAN_SRC})

    elseif(COORD_TYPE STREQUAL "2")
        # Spherical coordinate system
        message(STATUS "Using Spherical coordinate system")
        file(GLOB SPHERICAL_SRC "src/spherical/*.c")
        set(SRC_FILES ${COMMON_SRC_FILES} ${SPHERICAL_SRC})
    elseif(COORD_TYPE STREQUAL "3")
        message(STATUS "Using non-orthogonal Cartesian-state backend")
        file(GLOB NONORTHOGONAL_BACKEND_SRC "src/nonorthogonal_backend/*.c")
        # This controller is linked through
        # gamera_nonorthogonal_inner_magnetosphere when its experimental gate
        # is enabled; do not also compile it into the executable from the glob.
        list(REMOVE_ITEM NONORTHOGONAL_BACKEND_SRC
                "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_feedback.c"
                "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_feedback_restart_hdf5.c"
                "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_feedback_restart_mpi.c")
        # The prescribed H+ target loader is an experimental verification
        # path, not part of the production/default-OFF driver.  Keep its HDF5
        # and inner-magnetosphere dependencies completely out of an OFF build.
        if(NOT GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK)
            list(REMOVE_ITEM NONORTHOGONAL_BACKEND_SRC
                    "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_runtime.c")
        endif()
        if(NOT GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK AND
           NOT GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC AND
           NOT PROBLEM STREQUAL "22")
            list(REMOVE_ITEM NONORTHOGONAL_BACKEND_SRC
                    "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_sha256.c")
        endif()
        if(NOT GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC)
            list(REMOVE_ITEM NONORTHOGONAL_BACKEND_SRC
                    "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_checkpoint_manifest_v2.c"
                    "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_runtime.c"
                    "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_live_adapter.c"
                    "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_live_driver.c"
                    "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_checkpoint.c"
                    "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_restart_mpi.c")
            if(NOT PROBLEM STREQUAL "22")
                list(REMOVE_ITEM NONORTHOGONAL_BACKEND_SRC
                        "${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend/nonorthogonal_mi_restart.c")
            endif()
        endif()
        set(SRC_FILES ${COMMON_SRC_FILES} ${NONORTHOGONAL_BACKEND_SRC}
                src/cartesian/cartesian_RingAverage.c)
    else()
        message(FATAL_ERROR "Invalid or no COORD_TYPE specified. Use '0' Cartesian, '2' spherical, or '3' non-orthogonal.")
    endif()

    list(FILTER SRC_FILES EXCLUDE REGEX "/\\..*")
    # Add your main executable target
    add_executable(${PROJECT_NAME} ${SRC_FILES})  # Adjust paths as necessary
    if(GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK)
        if(NOT COORD_TYPE STREQUAL "3" OR NOT PROBLEM STREQUAL "22")
            message(FATAL_ERROR
                    "GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK is restricted to COORD_TYPE=3, PROBLEM=22")
        endif()
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK)
        target_link_libraries(${PROJECT_NAME}
                gamera_nonorthogonal_inner_magnetosphere)
        add_executable(test_nonorthogonal_sha256
                tests/test_nonorthogonal_sha256.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_include_directories(test_nonorthogonal_sha256 PRIVATE
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_sha256 PRIVATE
                    -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_sha256 COMMAND test_nonorthogonal_sha256)
        add_executable(test_nonorthogonal_prescribed_stream_loader
                tests/test_nonorthogonal_prescribed_stream_loader.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_runtime.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_compile_definitions(test_nonorthogonal_prescribed_stream_loader
                PRIVATE GAMERA_NO_PRESCRIBED_STREAM_LOADER_TEST_API)
        target_include_directories(test_nonorthogonal_prescribed_stream_loader
                PRIVATE
                ${CMAKE_SOURCE_DIR}/src
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(test_nonorthogonal_prescribed_stream_loader
                PRIVATE ${HDF5_LIBRARIES} MPI::MPI_C m)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_prescribed_stream_loader
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        find_package(Python3 REQUIRED COMPONENTS Interpreter)
        add_test(NAME nonorthogonal_prescribed_stream_loader
                COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tests/test_nonorthogonal_prescribed_stream_loader.py
                $<TARGET_FILE:test_nonorthogonal_prescribed_stream_loader>
                ${CMAKE_SOURCE_DIR}/replay/physical_boundary_64_6h_closed/component_matched_v1/cadence_states
                ${CMAKE_SOURCE_DIR}/replay/physical_boundary_64_6h_closed/maps)
    endif()
    if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC)
        if(NOT COORD_TYPE STREQUAL "3" OR NOT PROBLEM STREQUAL "22")
            message(FATAL_ERROR
                    "GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC is restricted to COORD_TYPE=3, PROBLEM=22")
        endif()
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC)
        if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK)
            target_compile_definitions(${PROJECT_NAME} PRIVATE
                    GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK)
            target_link_libraries(${PROJECT_NAME}
                    gamera_nonorthogonal_online_feedback_restart)
        endif()
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            set_property(SOURCE
                    src/nonorthogonal_backend/nonorthogonal_checkpoint_manifest_v2.c
                    src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_live_driver.c
                    src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_checkpoint.c
                    src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_restart_mpi.c
                    src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_runtime.c
                    src/nonorthogonal_backend/nonorthogonal_mi_restart.c
                    APPEND PROPERTY COMPILE_OPTIONS -fno-fast-math)
        endif()
        target_link_libraries(${PROJECT_NAME}
                gamera_nonorthogonal_inner_magnetosphere
                gamera_inner_magnetosphere_live_mpi
                gamera_inner_magnetosphere_snapshot_hdf5
                gamera_inner_magnetosphere_restart_hdf5)
        add_executable(test_nonorthogonal_inner_magnetosphere_online_runtime
                tests/test_nonorthogonal_inner_magnetosphere_online_runtime.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_runtime.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_include_directories(
                test_nonorthogonal_inner_magnetosphere_online_runtime PRIVATE
                ${CMAKE_SOURCE_DIR}/src
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
        target_link_libraries(
                test_nonorthogonal_inner_magnetosphere_online_runtime PRIVATE
                gamera_nonorthogonal_inner_magnetosphere MPI::MPI_C)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_nonorthogonal_inner_magnetosphere_online_runtime
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_inner_magnetosphere_online_runtime
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_nonorthogonal_inner_magnetosphere_online_runtime>
                        ${MPIEXEC_POSTFLAGS})
        add_executable(test_nonorthogonal_inner_magnetosphere_live_adapter_mpi
                tests/test_nonorthogonal_inner_magnetosphere_live_adapter_mpi.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_live_adapter.c)
        target_include_directories(
                test_nonorthogonal_inner_magnetosphere_live_adapter_mpi PRIVATE
                ${CMAKE_SOURCE_DIR}/src
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
        target_link_libraries(
                test_nonorthogonal_inner_magnetosphere_live_adapter_mpi PRIVATE
                gamera_nonorthogonal_kernels
                gamera_inner_magnetosphere_live_mpi MPI::MPI_C)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_nonorthogonal_inner_magnetosphere_live_adapter_mpi
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_inner_magnetosphere_live_adapter_2rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_nonorthogonal_inner_magnetosphere_live_adapter_mpi>
                        ${MPIEXEC_POSTFLAGS})
        # Compile the actual production live-driver translation unit against
        # narrow link-time mocks.  This exercises its MPI orchestration without
        # launching MHD or weakening the production implementation with a test
        # hook.
        add_executable(test_nonorthogonal_inner_magnetosphere_live_driver_mpi
                tests/test_nonorthogonal_inner_magnetosphere_live_driver_mpi.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_live_driver.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_compile_definitions(
                test_nonorthogonal_inner_magnetosphere_live_driver_mpi PRIVATE
                GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC)
        target_include_directories(
                test_nonorthogonal_inner_magnetosphere_live_driver_mpi PRIVATE
                ${CMAKE_SOURCE_DIR}/src
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
        target_link_libraries(
                test_nonorthogonal_inner_magnetosphere_live_driver_mpi PRIVATE
                MPI::MPI_C m)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_nonorthogonal_inner_magnetosphere_live_driver_mpi
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic
                            -Wno-strict-prototypes)
        endif()
        add_test(NAME nonorthogonal_inner_magnetosphere_live_driver_2rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_nonorthogonal_inner_magnetosphere_live_driver_mpi>
                        ${MPIEXEC_POSTFLAGS})
        if(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK)
            # Compile the same production live-driver translation unit with
            # the online feedback seam enabled.  Its narrow mocks exercise
            # fresh t=0 activation, collective prepare/commit, runtime-OFF
            # identity, and the deliberately unsupported restart path.
            add_executable(
                    test_nonorthogonal_inner_magnetosphere_live_driver_feedback_mpi
                    tests/test_nonorthogonal_inner_magnetosphere_live_driver_mpi.c
                    src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_live_driver.c
                    src/nonorthogonal_backend/nonorthogonal_sha256.c
                    src/inner_magnetosphere/gamera_inner_magnetosphere_online_hplus_coordinates.c)
            target_compile_definitions(
                    test_nonorthogonal_inner_magnetosphere_live_driver_feedback_mpi
                    PRIVATE
                    GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
                    GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK)
            target_include_directories(
                    test_nonorthogonal_inner_magnetosphere_live_driver_feedback_mpi
                    PRIVATE
                    ${CMAKE_SOURCE_DIR}/src
                    ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                    ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                    ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
            target_link_libraries(
                    test_nonorthogonal_inner_magnetosphere_live_driver_feedback_mpi
                    PRIVATE MPI::MPI_C m)
            if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
                target_compile_options(
                        test_nonorthogonal_inner_magnetosphere_live_driver_feedback_mpi
                        PRIVATE -fno-fast-math -Wall -Wextra -Werror
                                -pedantic -Wno-strict-prototypes)
            endif()
            add_test(
                    NAME nonorthogonal_inner_magnetosphere_live_driver_feedback_2rank
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                            ${MPIEXEC_PREFLAGS}
                            $<TARGET_FILE:test_nonorthogonal_inner_magnetosphere_live_driver_feedback_mpi>
                            ${MPIEXEC_POSTFLAGS})
        endif()
        add_executable(test_nonorthogonal_inner_magnetosphere_online_restart_mpi
                tests/test_nonorthogonal_inner_magnetosphere_online_restart_mpi.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_restart_mpi.c)
        target_include_directories(
                test_nonorthogonal_inner_magnetosphere_online_restart_mpi PRIVATE
                ${CMAKE_SOURCE_DIR}/src
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
        target_link_libraries(
                test_nonorthogonal_inner_magnetosphere_online_restart_mpi PRIVATE
                gamera_nonorthogonal_inner_magnetosphere MPI::MPI_C m)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_nonorthogonal_inner_magnetosphere_online_restart_mpi
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_inner_magnetosphere_online_restart_2rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_nonorthogonal_inner_magnetosphere_online_restart_mpi>
                        ${MPIEXEC_POSTFLAGS})
        add_executable(test_nonorthogonal_inner_magnetosphere_online_checkpoint_mpi
                tests/test_nonorthogonal_inner_magnetosphere_online_checkpoint_mpi.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_checkpoint.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_restart_mpi.c
                src/nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_runtime.c
                src/nonorthogonal_backend/nonorthogonal_checkpoint_manifest_v2.c)
        target_include_directories(
                test_nonorthogonal_inner_magnetosphere_online_checkpoint_mpi
                PRIVATE
                ${CMAKE_SOURCE_DIR}/src
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(
                test_nonorthogonal_inner_magnetosphere_online_checkpoint_mpi
                PRIVATE
                gamera_nonorthogonal_inner_magnetosphere
                gamera_inner_magnetosphere_restart_hdf5
                MPI::MPI_C ${HDF5_LIBRARIES} m)
        if(GAMERA_YINYANG_MFE_INTERFACE)
            target_compile_definitions(
                    test_nonorthogonal_inner_magnetosphere_online_checkpoint_mpi
                    PRIVATE GAMERA_YINYANG_MFE_INTERFACE)
        endif()
        if(GAMERA_YINYANG_HDIV_RECONCILE)
            target_compile_definitions(
                    test_nonorthogonal_inner_magnetosphere_online_checkpoint_mpi
                    PRIVATE GAMERA_YINYANG_HDIV_RECONCILE)
        endif()
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_nonorthogonal_inner_magnetosphere_online_checkpoint_mpi
                    PRIVATE -fno-fast-math -Wall -Wextra -Werror -pedantic
                            -Wno-strict-prototypes)
        endif()
        add_test(NAME nonorthogonal_inner_magnetosphere_online_checkpoint_2rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_nonorthogonal_inner_magnetosphere_online_checkpoint_mpi>
                        ${CMAKE_CURRENT_BINARY_DIR}
                        ${MPIEXEC_POSTFLAGS})
        add_executable(test_nonorthogonal_mi_restart_mpi
                tests/test_nonorthogonal_mi_restart_mpi.c
                src/nonorthogonal_backend/nonorthogonal_mi_restart.c
                src/nonorthogonal_backend/nonorthogonal_sha256.c)
        target_include_directories(test_nonorthogonal_mi_restart_mpi PRIVATE
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${HDF5_INCLUDE_DIRS})
        target_link_libraries(test_nonorthogonal_mi_restart_mpi PRIVATE
                MPI::MPI_C ${HDF5_LIBRARIES} m)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(test_nonorthogonal_mi_restart_mpi PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_mi_restart_2rank
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:test_nonorthogonal_mi_restart_mpi>
                        ${CMAKE_CURRENT_BINARY_DIR}
                        ${MPIEXEC_POSTFLAGS})
        add_executable(test_nonorthogonal_checkpoint_manifest_v2
                tests/test_nonorthogonal_checkpoint_manifest_v2.c
                src/nonorthogonal_backend/nonorthogonal_checkpoint_manifest_v2.c)
        target_include_directories(
                test_nonorthogonal_checkpoint_manifest_v2 PRIVATE
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal
                ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
                ${CMAKE_SOURCE_DIR}/src/inner_magnetosphere)
        if(NOT MSVC)
            target_link_libraries(
                    test_nonorthogonal_checkpoint_manifest_v2 PRIVATE m)
        endif()
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(
                    test_nonorthogonal_checkpoint_manifest_v2 PRIVATE
                    -fno-fast-math -Wall -Wextra -Werror -pedantic)
        endif()
        add_test(NAME nonorthogonal_checkpoint_manifest_v2
                COMMAND test_nonorthogonal_checkpoint_manifest_v2)
        find_package(Python3 REQUIRED COMPONENTS Interpreter)
        add_test(NAME nonorthogonal_online_restart_frontend_failstop
                COMMAND ${Python3_EXECUTABLE}
                        ${CMAKE_SOURCE_DIR}/tests/test_nonorthogonal_online_restart_frontend.py
                        ${MPIEXEC_EXECUTABLE}
                        ${MPIEXEC_NUMPROC_FLAG}
                        $<TARGET_FILE:${PROJECT_NAME}>)
        set_tests_properties(
                nonorthogonal_inner_magnetosphere_online_runtime
                nonorthogonal_inner_magnetosphere_live_adapter_2rank
                nonorthogonal_inner_magnetosphere_live_driver_2rank
                nonorthogonal_inner_magnetosphere_online_restart_2rank
                nonorthogonal_inner_magnetosphere_online_checkpoint_2rank
                nonorthogonal_checkpoint_manifest_v2
                nonorthogonal_online_restart_frontend_failstop
                PROPERTIES TIMEOUT 60)
    endif()
    if(GAMERA_STRICT_INNER_MAGNETIC_WALL)
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_STRICT_INNER_MAGNETIC_WALL)
    endif()
    if(GAMERA_STRICT_INNER_GAS_WALL)
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_STRICT_INNER_GAS_WALL)
    endif()
    if(GAMERA_YINYANG_MFE_INTERFACE)
        # Experimental overset-interface mode: exchange primitive fluid states
        # only into angular ghosts and let each patch compute/evolve its own
        # active edge EMF with constrained transport, following Yin-Yang MFE.
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_YINYANG_MFE_INTERFACE)
    endif()
    if(GAMERA_YINYANG_HDIV_RECONCILE)
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_YINYANG_HDIV_RECONCILE)
        if(GAMERA_YINYANG_HDIV_OPTIMIZED)
            target_compile_definitions(${PROJECT_NAME} PRIVATE
                    GAMERA_YINYANG_HDIV_OPTIMIZED)
        endif()
        if(GAMERA_YINYANG_HDIV_DISTRIBUTED)
            target_compile_definitions(${PROJECT_NAME} PRIVATE
                    GAMERA_YINYANG_HDIV_DISTRIBUTED)
        endif()
        if(GAMERA_YINYANG_HDIV_PROFILE)
            target_compile_definitions(${PROJECT_NAME} PRIVATE
                    GAMERA_YINYANG_HDIV_PROFILE)
        endif()
        if(GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY)
            target_compile_definitions(${PROJECT_NAME} PRIVATE
                    GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY)
        endif()
    endif()
    if(GAMERA_YINYANG_SPARSE_OVERSET)
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_YINYANG_SPARSE_OVERSET)
    endif()
    if(GAMERA_YINYANG_SPARSE_OVERSET_PROFILE)
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_YINYANG_SPARSE_OVERSET_PROFILE)
    endif()
    if(GAMERA_BENCH_PROFILE_JSONL)
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_BENCH_PROFILE_JSONL)
    endif()
    target_compile_definitions(${PROJECT_NAME} PRIVATE
            GAMERA_COORDINATE_BACKEND_ID=${COORD_TYPE}
            GAMERA_PROBLEM_ID=${PROBLEM}
            GAMERA_BOW_SHOCK_INNER_RADIUS_RE=${GAMERA_BOW_SHOCK_INNER_RADIUS_RE}
            GAMERA_BOW_SHOCK_OUTER_RADIUS_RE=${GAMERA_BOW_SHOCK_OUTER_RADIUS_RE}
            GAMERA_BOW_SHOCK_RADIAL_STRETCH=${GAMERA_BOW_SHOCK_RADIAL_STRETCH}
            GAMERA_EARTH_RADIAL_MAP_VERSION=${GAMERA_EARTH_RADIAL_MAP_VERSION})

    # Resolve HDF5 before adding MPI headers.  Anaconda's MPI prefix may ship a
    # second hdf5.h, so the selected HDF5 header must win include lookup.
    find_package(HDF5 REQUIRED COMPONENTS C)
    target_include_directories(${PROJECT_NAME} PRIVATE
            ${HDF5_INCLUDE_DIRS}
            ${CMAKE_SOURCE_DIR}/src
            ${CMAKE_SOURCE_DIR}/logc
            ${CMAKE_SOURCE_DIR}/math_lib
            ${CMAKE_SOURCE_DIR}/src/nonorthogonal
            ${CMAKE_SOURCE_DIR}/src/nonorthogonal_backend
    )
    target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE
            ${MPI_C_INCLUDE_PATH})

    # Define the appropriate macro for the selected problem
    if(COORD_TYPE STREQUAL "3" AND NOT (PROBLEM STREQUAL "6" OR PROBLEM STREQUAL "7" OR PROBLEM STREQUAL "10" OR PROBLEM STREQUAL "11" OR PROBLEM STREQUAL "12" OR PROBLEM STREQUAL "13" OR PROBLEM STREQUAL "14" OR PROBLEM STREQUAL "15" OR PROBLEM STREQUAL "16" OR PROBLEM STREQUAL "17" OR PROBLEM STREQUAL "18" OR PROBLEM STREQUAL "19" OR PROBLEM STREQUAL "20" OR PROBLEM STREQUAL "21" OR PROBLEM STREQUAL "22"))
        message(FATAL_ERROR "COORD_TYPE=3 accepts PROBLEM=6 (warped OT), 7 (spherical 3D blast), and 10-22 (Yin-Yang validation/application ladder).")
    elseif((PROBLEM STREQUAL "6" OR PROBLEM STREQUAL "7" OR PROBLEM STREQUAL "10" OR PROBLEM STREQUAL "11" OR PROBLEM STREQUAL "12" OR PROBLEM STREQUAL "13" OR PROBLEM STREQUAL "14" OR PROBLEM STREQUAL "15" OR PROBLEM STREQUAL "16" OR PROBLEM STREQUAL "17" OR PROBLEM STREQUAL "18" OR PROBLEM STREQUAL "19" OR PROBLEM STREQUAL "20" OR PROBLEM STREQUAL "21" OR PROBLEM STREQUAL "22") AND NOT COORD_TYPE STREQUAL "3")
        message(FATAL_ERROR "PROBLEM=6/7/10/11/12/13/14/15/16/17/18/19/20/21 require COORD_TYPE=3.")
    elseif(PROBLEM STREQUAL "9" AND NOT COORD_TYPE STREQUAL "2")
        message(FATAL_ERROR "PROBLEM=9 (spherical HD pole blast) requires COORD_TYPE=2.")
    elseif(PROBLEM STREQUAL "1")
        target_compile_definitions(${PROJECT_NAME} PRIVATE PROBLEM_OT_cartesian)
        message(STATUS "Compiling for PROBLEM 1: PROBLEM_OT_cartesian")
    elseif(PROBLEM STREQUAL "2")
        target_compile_definitions(${PROJECT_NAME} PRIVATE PROBLEM_FieldLoop_cartesian)
        message(STATUS "Compiling for PROBLEM 2: PROBLEM_FieldLoop_cartesian")
    elseif(PROBLEM STREQUAL "3")
        target_compile_definitions(${PROJECT_NAME} PRIVATE PROBLEM_FieldLoop_spherical)
        message(STATUS "Compiling for PROBLEM 3: PROBLEM_FieldLoop_spherical")
    elseif(PROBLEM STREQUAL "4")
        target_compile_definitions(${PROJECT_NAME} PRIVATE PROBLEM_Blast_cartesian)
        message(STATUS "Compiling for PROBLEM 4: PROBLEM_Blast_cartesian")
    elseif(PROBLEM STREQUAL "5")
        target_compile_definitions(${PROJECT_NAME} PRIVATE PROBLEM_Blast_spherical)
        message(STATUS "Compiling for PROBLEM 5: PROBLEM_Blast_spherical")
    elseif(PROBLEM STREQUAL "6")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND PROBLEM_NONORTHOGONAL_OT)
        message(STATUS "Compiling for PROBLEM 6: warped non-orthogonal OT")
    elseif(PROBLEM STREQUAL "7")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                PROBLEM_NONORTHOGONAL_BLAST_SPHERICAL)
        message(STATUS "Compiling for PROBLEM 7: spherical non-orthogonal 3D blast")
    elseif(PROBLEM STREQUAL "9")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                PROBLEM_HD_BLAST_SPHERICAL_RING)
        message(STATUS "Compiling for PROBLEM 9: spherical HD pole blast with Ring Average")
    elseif(PROBLEM STREQUAL "10")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                PROBLEM_NONORTHOGONAL_YINYANG_HD_BLAST)
        message(STATUS "Compiling for PROBLEM 10: non-orthogonal Yin-Yang HD blast")
    elseif(PROBLEM STREQUAL "11")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                PROBLEM_NONORTHOGONAL_YINYANG_HD_BLAST_LONG)
        message(STATUS "Compiling for PROBLEM 11: travelling non-orthogonal Yin-Yang HD blast")
    elseif(PROBLEM STREQUAL "12")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_NONORTHOGONAL_HAS_CELL_SOURCE
                PROBLEM_NONORTHOGONAL_YINYANG_ENTROPY)
        message(STATUS "Compiling for PROBLEM 12: one-orbit Yin-Yang spherical entropy wave")
    elseif(PROBLEM STREQUAL "13")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                PROBLEM_NONORTHOGONAL_YINYANG_HD_OBSTACLE)
        message(STATUS "Compiling for PROBLEM 13: strict-HD Yin-Yang sphere obstacle")
    elseif(PROBLEM STREQUAL "14")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                PROBLEM_NONORTHOGONAL_YINYANG_MHD_UNIFORM)
        message(STATUS "Compiling for PROBLEM 14: uniform Yin-Yang MHD/CT gate")
    elseif(PROBLEM STREQUAL "15")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
                GAMERA_TIME_DEPENDENT_WIND
                PROBLEM_NONORTHOGONAL_YINYANG_MHD_OBSTACLE)
        message(STATUS "Compiling for PROBLEM 15: Yin-Yang IMF sphere obstacle")
    elseif(PROBLEM STREQUAL "16")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_CELL_SOURCE
                GAMERA_TIME_DEPENDENT_WIND
                GAMERA_IONOSPHERE_SOURCE_LOSS
                PROBLEM_NONORTHOGONAL_YINYANG_IONOSPHERE)
        message(STATUS "Compiling for PROBLEM 16: Yin-Yang ionosphere source/loss")
    elseif(PROBLEM STREQUAL "17")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_CELL_SOURCE
                GAMERA_TIME_DEPENDENT_WIND
                GAMERA_IONOSPHERE_SOURCE_LOSS
                GAMERA_CRUSTAL_FIELD
                PROBLEM_NONORTHOGONAL_YINYANG_CRUSTAL)
        message(STATUS "Compiling for PROBLEM 17: Yin-Yang localized crustal field")
    elseif(PROBLEM STREQUAL "18")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_ADJUST
                GAMERA_NONORTHOGONAL_USE_BORIS
                GAMERA_TIME_DEPENDENT_WIND
                GAMERA_EARTH_DIPOLE_BACKGROUND
                PROBLEM_NONORTHOGONAL_EARTH_MAGNETOSPHERE)
        message(STATUS "Compiling for PROBLEM 18: Yin-Yang Earth dipole background-field magnetosphere")
    elseif(PROBLEM STREQUAL "19")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_RADIAL_MAP
                GAMERA_NONORTHOGONAL_USE_BORIS
                GAMERA_TIME_DEPENDENT_WIND
                GAMERA_BOW_SHOCK
                PROBLEM_NONORTHOGONAL_YINYANG_BOW_SHOCK)
        message(STATUS "Compiling for PROBLEM 19: Yin-Yang Mach-10 unmagnetized bow shock")
    elseif(PROBLEM STREQUAL "20")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_RADIAL_MAP
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_ADJUST
                GAMERA_NONORTHOGONAL_USE_BORIS
                GAMERA_TIME_DEPENDENT_WIND
                GAMERA_EARTH_DIPOLE_BACKGROUND
                GAMERA_BOW_SHOCK
                GAMERA_EARTH_MAGNETOSPHERE_STRETCHED
                PROBLEM_NONORTHOGONAL_EARTH_MAGNETOSPHERE_STRETCHED)
        message(STATUS "Compiling for PROBLEM 20: stretched-grid Yin-Yang Earth magnetosphere")
    elseif(PROBLEM STREQUAL "21")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_RADIAL_MAP
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_ADJUST
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_RESIDUAL_INIT
                GAMERA_NONORTHOGONAL_USE_BORIS
                GAMERA_TIME_DEPENDENT_WIND
                GAMERA_EARTH_DIPOLE_BACKGROUND
                GAMERA_BOW_SHOCK
                GAMERA_EARTH_MAGNETOSPHERE_STRETCHED
                GAMERA_EARTH_UPSTREAM_STARTUP
                PROBLEM_NONORTHOGONAL_EARTH_MAGNETOSPHERE_STRETCHED)
        message(STATUS "Compiling for PROBLEM 21: upstream-launched Yin-Yang Earth magnetosphere")
    elseif(PROBLEM STREQUAL "22")
        target_compile_definitions(${PROJECT_NAME} PRIVATE
                GAMERA_NONORTHOGONAL_BACKEND
                GAMERA_YINYANG_BACKEND
                GAMERA_YINYANG_MHD_BACKEND
                GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
                GAMERA_NONORTHOGONAL_HAS_RADIAL_MAP
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_ADJUST
                GAMERA_NONORTHOGONAL_HAS_BACKGROUND_RESIDUAL_INIT
                GAMERA_NONORTHOGONAL_USE_BORIS
                GAMERA_TIME_DEPENDENT_WIND
                GAMERA_EARTH_DIPOLE_BACKGROUND
                GAMERA_BOW_SHOCK
                GAMERA_EARTH_MAGNETOSPHERE_STRETCHED
                GAMERA_EARTH_UPSTREAM_STARTUP
                GAMERA_MI_COUPLING
                PROBLEM_NONORTHOGONAL_EARTH_MAGNETOSPHERE_STRETCHED)
        message(STATUS "Compiling for PROBLEM 22: upstream-launched electrostatic M-I-coupled Yin-Yang Earth magnetosphere")
    else()
        message(FATAL_ERROR "Invalid PROBLEM. Use 1-5 for legacy, 6-7/10-22 for non-orthogonal tests, or 9 for the spherical HD Ring Average test.")
    endif()

    if(COORD_TYPE STREQUAL "3")
        target_link_libraries(${PROJECT_NAME} gamera_nonorthogonal_kernels)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(${PROJECT_NAME} PRIVATE -fno-fast-math)
        endif()
    endif()

    # Link against MPI libraries
    target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES})  # For C

    # Link against OpenMP libraries
    if(APPLE)
        target_compile_options(${PROJECT_NAME} PRIVATE "-Xpreprocessor" "-fopenmp")
        target_include_directories(${PROJECT_NAME} PRIVATE ${OMP_INCLUDE_DIR})
        target_link_libraries(${PROJECT_NAME} ${OMP_LIB})
    else()
        target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_C)
    endif()

    # Link against HDF5 libraries
    target_link_libraries(${PROJECT_NAME} ${HDF5_C_LIBRARIES})
    # Anaconda's MPI include directory can also contain an unrelated hdf5.h.
    # Put the include directory belonging to the selected HDF5 library first,
    # otherwise the executable can compile against one ABI and load another.

    # Copy compile_commands.json to source directory after build
    add_custom_command(
            TARGET ${PROJECT_NAME}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/compile_commands.json ${CMAKE_SOURCE_DIR}/
    )
else()
    message(FATAL_ERROR "MPI was not found, but it is required.")
endif()
