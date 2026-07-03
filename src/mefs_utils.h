#ifndef MEFS_UTILS_H
#define MEFS_UTILS_H

#include <fftw3.h>

/* ANSI color escape sequences */
extern const char *c_reset;
extern const char *c_bold;
extern const char *c_cyan;
extern const char *c_b_green;
extern const char *c_green;
extern const char *c_magenta;
extern const char *c_yellow;
extern const char *c_b_red;
extern const char *c_grey;

/**
 * @brief Initialize ANSI color code variables based on the environment.
 */
void init_colors(void);

/**
 * @brief Print a progress bar to stdout.
 *
 * @param count Current completed point index.
 * @param total Total number of points to process.
 */
void print_progress_bar(
    int count,
    int total);

/**
 * @brief Check if file A has a newer modification time than file B.
 *
 * @param file_a Path to file A.
 * @param file_b Path to file B.
 * @return int 1 if file A is newer, 0 otherwise.
 */
int is_file_newer(
    const char *file_a,
    const char *file_b);

/**
 * @brief Shift the zero-frequency component to the center of the spectrum.
 *
 * @param in Input complex grid data.
 * @param out Output shifted complex grid data.
 * @param rows Number of rows.
 * @param cols Number of columns.
 */
void fftshift(
    const fftw_complex *in,
    fftw_complex       *out,
    int                 rows,
    int                 cols);

/**
 * @brief Solve the linear system G * c = b using Gaussian elimination.
 *
 * @param G System matrix of size N x N (flattened, modified in place).
 * @param b Right-hand side vector of size N (holds solution c on return).
 * @param N Dimension of the system.
 * @return int 0 on success, -1 if the matrix is singular.
 */
int solve_linear_system(
    double *G,
    double *b,
    int     N);

#endif /* MEFS_UTILS_H */
