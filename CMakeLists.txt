cmake_minimum_required(VERSION 3.10)
project(MEFS C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find dependencies using PkgConfig
find_package(PkgConfig REQUIRED)
pkg_check_modules(CFITSIO REQUIRED cfitsio)
pkg_check_modules(FFTW3 REQUIRED fftw3)

add_executable(mefs_psf
    src/mefs_psf.c
    src/mefs_utils.c
    src/mefs_fits.c
    src/mefs_optics.c
    src/mefs_simulation.c
)

# Include paths
target_include_directories(mefs_psf PRIVATE
    ${CFITSIO_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
)

# Link directories and libraries
target_link_directories(mefs_psf PRIVATE
    ${CFITSIO_LIBRARY_DIRS}
    ${FFTW3_LIBRARY_DIRS}
)

find_package(OpenMP REQUIRED)

target_link_libraries(mefs_psf PRIVATE
    ${CFITSIO_LIBRARIES}
    ${FFTW3_LIBRARIES}
    OpenMP::OpenMP_C
    m
)

# Scene generator tool
add_executable(mefs_create_scene src/mefs_create_scene.c)
target_link_libraries(mefs_create_scene PRIVATE m)

# Score grading tool
add_executable(mefs_score src/mefs_score.c)
target_include_directories(mefs_score PRIVATE
    ${CFITSIO_INCLUDE_DIRS}
)
target_link_directories(mefs_score PRIVATE
    ${CFITSIO_LIBRARY_DIRS}
)
target_link_libraries(mefs_score PRIVATE
    ${CFITSIO_LIBRARIES}
    m
)
