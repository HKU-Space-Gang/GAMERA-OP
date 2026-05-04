cmake_minimum_required(VERSION 3.10)
project(run_mhd LANGUAGES C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(APPLE)
    # add include
    include_directories(/opt/homebrew/include)

    # link directories
    link_directories(/opt/homebrew/lib)
endif()

# Find MPI and require at least MPI 1.2
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find MPI and require at least MPI 1.2
find_package(MPI REQUIRED COMPONENTS C)  # Cxx is for C++, omit if not needed

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
    else()
        message(FATAL_ERROR "Invalid or no COORD_TYPE specified. Please set COORD_TYPE to '0' for Cartesian or '2' for Spherical.")
    endif()

    list(FILTER SRC_FILES EXCLUDE REGEX "/\\..*")
    # Add your main executable target
    add_executable(${PROJECT_NAME} ${SRC_FILES})  # Adjust paths as necessary

    target_include_directories(${PROJECT_NAME} PRIVATE
            ${CMAKE_SOURCE_DIR}/src
            ${CMAKE_SOURCE_DIR}/logc
            ${CMAKE_SOURCE_DIR}/math_lib
            ${MPI_C_INCLUDE_PATH}  # 添加 MPI 头文件路径
    )

    # Define the appropriate macro for the selected problem
    if(PROBLEM STREQUAL "1")
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
    else()
        message(FATAL_ERROR "Invalid or no PROBLEM specified. Please set PROBLEM to '1', '2', '3', '4', or '5'.")
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
    find_package(HDF5 REQUIRED COMPONENTS C)
    target_link_libraries(${PROJECT_NAME} ${HDF5_C_LIBRARIES})
    target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})  # 添加 HDF5 头文件路径

    # Copy compile_commands.json to source directory after build
    add_custom_command(
            TARGET ${PROJECT_NAME}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/compile_commands.json ${CMAKE_SOURCE_DIR}/
    )
else()
    message(FATAL_ERROR "MPI was not found, but it is required.")
endif()
