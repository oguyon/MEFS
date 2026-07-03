#include "mefs_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <fftw3.h>

/* Define global color code pointers */
const char *c_reset   = "";
const char *c_bold    = "";
const char *c_cyan    = "";
const char *c_b_green = "";
const char *c_green   = "";
const char *c_magenta = "";
const char *c_yellow  = "";
const char *c_b_red   = "";
const char *c_grey    = "";

/**
 * @brief Initialize ANSI color code variables based on the environment.
 */
void init_colors(void)
{
    if (getenv("NO_COLOR") == NULL)
    {
        c_reset   = "\x1b[0m";
        c_bold    = "\x1b[1m";
        c_cyan    = "\x1b[1;36m";
        c_b_green = "\x1b[1;32m";
        c_green   = "\x1b[32m";
        c_magenta = "\x1b[35m";
        c_yellow  = "\x1b[33m";
        c_b_red   = "\x1b[1;31m";
        c_grey    = "\x1b[90m";
    }
}

/**
 * @brief Print a progress bar to stdout.
 */
void print_progress_bar(
    int count,
    int total)
{
    int width = 40;
    double progress = (double)(count + 1) / total;
    int filled = (int)(progress * width);

    printf("\r%sProgress:%s %s[%s", c_cyan, c_reset, c_grey, c_reset);
    for (int ii = 0; ii < width; ii++)
    {
        if (ii < filled)
        {
            printf("%s=%s", c_green, c_reset);
        }
        else if (ii == filled)
        {
            printf("%s>%s", c_green, c_reset);
        }
        else
        {
            printf(" ");
        }
    }
    printf("%s]%s %d/%d (%.1f%%)", c_grey, c_reset, count + 1, total, progress * 100.0);
    fflush(stdout);
} // void print_progress_bar

/**
 * @brief Check if file A has a newer modification time than file B.
 */
int is_file_newer(
    const char *file_a,
    const char *file_b)
{
    struct stat stat_a;
    struct stat stat_b;
    if (stat(file_a, &stat_a) != 0)
    {
        return 0;
    }
    if (stat(file_b, &stat_b) != 0)
    {
        return 1;
    }
    return stat_a.st_mtime > stat_b.st_mtime;
} // int is_file_newer

/**
 * @brief Shift the zero-frequency component to the center of the spectrum.
 */
void fftshift(
    const fftw_complex *in,
    fftw_complex       *out,
    int                 rows,
    int                 cols)
{
    int r_shift = rows / 2;
    int c_shift = cols / 2;

    for (int ii = 0; ii < rows; ii++)
    {
        int target_r = (ii + r_shift) % rows;
        for (int jj = 0; jj < cols; jj++)
        {
            int target_c = (jj + c_shift) % cols;
            int src_idx = ii * cols + jj;
            int dst_idx = target_r * cols + target_c;

            out[dst_idx][0] = in[src_idx][0];
            out[dst_idx][1] = in[src_idx][1];
        }
    }
} // void fftshift

/**
 * @brief Solve the linear system G * c = b using Gaussian elimination.
 */
int solve_linear_system(
    double *G,
    double *b,
    int     N)
{
    for (int ii = 0; ii < N; ii++)
    {
        /* Find pivot row */
        int pivot = ii;
        double max_val = fabs(G[ii * N + ii]);
        for (int row = ii + 1; row < N; row++)
        {
            double val = fabs(G[row * N + ii]);
            if (val > max_val)
            {
                max_val = val;
                pivot = row;
            }
        }

        /* Swap rows in G and b */
        if (pivot != ii)
        {
            for (int col = 0; col < N; col++)
            {
                double temp = G[ii * N + col];
                G[ii * N + col] = G[pivot * N + col];
                G[pivot * N + col] = temp;
            }
            double temp = b[ii];
            b[ii] = b[pivot];
            b[pivot] = temp;
        }

        /* Check singularity */
        if (fabs(G[ii * N + ii]) < 1e-12)
        {
            return -1; /* Singular matrix */
        }

        /* Eliminate columns below diagonal */
        for (int row = ii + 1; row < N; row++)
        {
            double factor = G[row * N + ii] / G[ii * N + ii];
            for (int col = ii; col < N; col++)
            {
                G[row * N + col] -= factor * G[ii * N + col];
            }
            b[row] -= factor * b[ii];
        }
    } // for (int ii = 0; ii < N; ii++)

    /* Back-substitution */
    for (int ii = N - 1; ii >= 0; ii--)
    {
        for (int col = ii + 1; col < N; col++)
        {
            b[ii] -= G[ii * N + col] * b[col];
        }
        b[ii] /= G[ii * N + ii];
    } // for (int ii = N - 1; ii >= 0; ii--)

    return 0;
} // int solve_linear_system
