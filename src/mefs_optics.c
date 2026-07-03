#include "mefs_optics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

/**
 * @brief Calculate starlight leakage subtraction weights for the coronagraph.
 */
void get_coronagraph_weights(
    double  coronagraph_param,
    double *weights)
{
    for (int m = 0; m < 6; m++)
    {
        weights[m] = 0.0;
    }

    if (coronagraph_param <= 0.0)
    {
        return;
    }

    if (coronagraph_param <= 1.0)
    {
        weights[0] = coronagraph_param;
    }
    else if (coronagraph_param <= 2.0)
    {
        double x = coronagraph_param - 1.0;
        weights[0] = 1.0;
        weights[1] = x;
        weights[2] = x;
    }
    else
    {
        double x = coronagraph_param - 2.0;
        if (x > 1.0)
        {
            x = 1.0;
        }
        weights[0] = 1.0;
        weights[1] = 1.0;
        weights[2] = 1.0;
        weights[3] = x;
        weights[4] = x;
        weights[5] = x;
    }
} // void get_coronagraph_weights

/**
 * @brief Generate a circular pupil with sub-pixel anti-aliasing.
 */
fftw_complex *generate_circular_pupil(
    int    rows,
    int    cols,
    double diameter)
{
    fftw_complex *pupil = fftw_alloc_complex(rows * cols);
    if (pupil == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Failed to allocate memory for circular pupil.\n",
                __func__, __LINE__);
        return NULL;
    }

    double cx = (cols - 1) / 2.0;
    double cy = (rows - 1) / 2.0;
    double R = diameter / 2.0;

    for (int ii = 0; ii < rows; ii++)
    {
        double dy = ii - cy;
        for (int jj = 0; jj < cols; jj++)
        {
            double dx = jj - cx;
            double dist = sqrt(dx * dx + dy * dy);
            int idx = ii * cols + jj;

            if (dist < R - 0.7072)
            {
                pupil[idx][0] = 1.0;
            }
            else if (dist > R + 0.7072)
            {
                pupil[idx][0] = 0.0;
            }
            else
            {
                /* Edge pixel: Integrate sub-pixels (20x20 grid) */
                int count = 0;
                int K = 20;

                for (int ky = 0; ky < K; ky++)
                {
                    double sy = ii - 0.5 + (ky + 0.5) / K - cy;
                    for (int kx = 0; kx < K; kx++)
                    {
                        double sx = jj - 0.5 + (kx + 0.5) / K - cx;
                        if (sx * sx + sy * sy <= R * R)
                        {
                            count++;
                        }
                    }
                }
                pupil[idx][0] = (double)count / (K * K);
            }
            pupil[idx][1] = 0.0;
        }
    } // for (int ii = 0; ii < rows; ii++)
    return pupil;
} // fftw_complex *generate_circular_pupil

/**
 * @brief Estimate the pupil diameter in pixels from the amplitude of the input field.
 */
double estimate_pupil_diameter(
    const fftw_complex *pupil,
    int                 rows,
    int                 cols)
{
    int min_r = rows;
    int max_r = -1;
    int min_c = cols;
    int max_c = -1;

    for (int ii = 0; ii < rows; ii++)
    {
        for (int jj = 0; jj < cols; jj++)
        {
            int idx = ii * cols + jj;
            double r = pupil[idx][0];
            double im = pupil[idx][1];
            double amp = sqrt(r * r + im * im);
            if (amp > 0.05) /* 5% threshold */
            {
                if (ii < min_r) min_r = ii;
                if (ii > max_r) max_r = ii;
                if (jj < min_c) min_c = jj;
                if (jj > max_c) max_c = jj;
            }
        }
    } // for (int ii = 0; ii < rows; ii++)

    if (max_r == -1 || max_c == -1)
    {
        return -1.0;
    }

    double diam_r = max_r - min_r + 1;
    double diam_c = max_c - min_c + 1;
    return (diam_r + diam_c) / 2.0;
} // double estimate_pupil_diameter
