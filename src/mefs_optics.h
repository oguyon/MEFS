#ifndef MEFS_OPTICS_H
#define MEFS_OPTICS_H

#include <fftw3.h>

/**
 * @brief Generate a circular pupil with sub-pixel anti-aliasing.
 *
 * @param rows Number of rows in the grid.
 * @param cols Number of columns in the grid.
 * @param diameter Pupil diameter in pixels.
 * @return fftw_complex* Allocated complex array containing the pupil amplitude, or NULL.
 */
fftw_complex *generate_circular_pupil(
    int    rows,
    int    cols,
    double diameter);

/**
 * @brief Estimate the pupil diameter in pixels from the amplitude of the input field.
 *
 * @param pupil Complex grid wavefront representation.
 * @param rows Number of rows.
 * @param cols Number of columns.
 * @return double Estimated pupil diameter in pixels, or -1.0 if not found.
 */
double estimate_pupil_diameter(
    const fftw_complex *pupil,
    int                 rows,
    int                 cols);

/**
 * @brief Calculate starlight leakage subtraction weights for the coronagraph.
 *
 * @param coronagraph_param Continuous parameter (0.0 to 3.0).
 * @param weights Output double array of size 6 containing computed mode subtraction weights.
 */
void get_coronagraph_weights(
    double  coronagraph_param,
    double *weights);

#endif /* MEFS_OPTICS_H */
