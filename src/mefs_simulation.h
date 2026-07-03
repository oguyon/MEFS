#ifndef MEFS_SIMULATION_H
#define MEFS_SIMULATION_H

#include <fftw3.h>
#include <stddef.h>

/**
 * @brief Accumulate the IFS intensity pattern for a single point source.
 */
void accumulate_ifs_source(
    const fftw_complex *psf,
    int                 grid_size,
    int                 L,
    int                 M,
    int                 P,
    fftw_complex       *lenslet_padded,
    fftw_complex       *lenslet_padded_shifted,
    fftw_complex       *farfield_shifted,
    fftw_complex       *farfield,
    fftw_plan           ifs_plan,
    double              pinhole_diam_ldl,
    double             *ifsraw_im_accum);

/**
 * @brief Run the end-to-end IFS simulation for a single point source.
 */
int run_ifs_simulation(
    const fftw_complex *psf,
    int                 grid_size,
    double              diameter,
    double              lenslet_size_ld,
    int                 lenslet_count,
    int                 ifs_fft_size,
    const char         *ifs_complex_file,
    const char         *ifs_intensity_file,
    const char         *ifs_flux_file,
    double              pinhole_diam_ldl,
    const char         *ifs_mask_file,
    double              init_pupil_energy);

/**
 * @brief Compute the complex MEFS reference mode.
 */
void compute_mefs_reference_mode(
    const fftw_complex *pupil,
    int                 rows,
    int                 cols,
    double              diameter,
    double              coronagraph_param,
    double              mefs_x,
    double              mefs_y,
    fftw_complex       *mefs_mode);

/**
 * @brief Construct the parameter verification string for MEFS files.
 */
void get_mefs_param_string(
    char   *buf,
    size_t  buf_len,
    double  coro_param,
    double  diameter,
    int     grid_size,
    double  mefs_x,
    double  mefs_y);

/**
 * @brief Check if cached MEFS reference parameters match current execution.
 */
int check_mefs_file_parameters(
    const char *filename,
    double      coro_param,
    double      diameter,
    int         grid_size,
    double      mefs_x,
    double      mefs_y);

/**
 * @brief Construct the parameter verification string for scene files.
 */
void get_scene_param_string(
    char   *buf,
    size_t  buf_len,
    double  coro_param,
    double  diameter,
    int     grid_size,
    int     enable_ifs,
    double  lenslet_size,
    int     lenslet_count,
    int     ifs_fft_size,
    double  pinhole_diam);

/**
 * @brief Check if cached scene log parameters match current execution.
 */
int check_scene_file_parameters(
    const char *filename,
    double      coro_param,
    double      diameter,
    int         grid_size,
    int         enable_ifs,
    double      lenslet_size,
    int         lenslet_count,
    int         ifs_fft_size,
    double      pinhole_diam);

/**
 * @brief Load computed flux sums from existing MEFS text output files.
 */
int load_mefs_scene_data(
    const char *filename,
    double     *out_uncoupled,
    double     *out_coupled);

/**
 * @brief Run simulation in MEFS matched-filter mode for a scene file.
 */
int run_mefs_scene_simulation(
    const char   *scene_path,
    const char   *name_suffix,
    fftw_complex *pupil,
    int           rows,
    int           cols,
    double        diameter,
    int           grid_size,
    double        coronagraph_param,
    double        mefs_x,
    double        mefs_y,
    double        init_pupil_energy,
    double       *out_total_uncoupled_flux,
    double       *out_total_coupled_flux);

/**
 * @brief Run simulation (IFS or MEFS mode) for an individual scene file or all scenes.
 */
int run_scene_simulation_file(
    const char   *scene_path,
    const char   *name_suffix,
    fftw_complex *pupil,
    int           rows,
    int           cols,
    double        diameter,
    int           grid_size,
    int           enable_ifs,
    double        lenslet_size_ld,
    int           lenslet_count,
    int           ifs_fft_size,
    double        pinhole_diam_ldl,
    const char   *ifs_mask_file,
    double        coronagraph_param,
    double        init_pupil_energy,
    int           run_mefs_mode,
    double        mefs_x,
    double        mefs_y);

#endif /* MEFS_SIMULATION_H */
