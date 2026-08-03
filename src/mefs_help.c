/**
 * @file mefs_help.c
 * @brief Top-level help page for the MEFS toolkit.
 *
 * Prints a colorful overview of all executables, a typical
 * simulation workflow, build instructions, and a condensed
 * option reference for each tool.
 *
 * Respects the NO_COLOR environment variable
 * (see https://no-color.org).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── ANSI color codes ─────────────────────────── */

static const char *c_reset   = "";
static const char *c_bold    = "";
static const char *c_dim     = "";
static const char *c_cyan    = "";
static const char *c_b_cyan  = "";
static const char *c_b_green = "";
static const char *c_green   = "";
static const char *c_magenta = "";
static const char *c_yellow  = "";
static const char *c_b_red   = "";
static const char *c_grey    = "";
static const char *c_blue    = "";
static const char *c_b_blue  = "";
static const char *c_white   = "";

static void init_colors(void)
{
    if (getenv("NO_COLOR") == NULL)
    {
        c_reset   = "\x1b[0m";
        c_bold    = "\x1b[1m";
        c_dim     = "\x1b[2m";
        c_cyan    = "\x1b[1;36m";
        c_b_cyan  = "\x1b[96m";
        c_b_green = "\x1b[1;32m";
        c_green   = "\x1b[32m";
        c_magenta = "\x1b[35m";
        c_yellow  = "\x1b[33m";
        c_b_red   = "\x1b[1;31m";
        c_grey    = "\x1b[90m";
        c_blue    = "\x1b[34m";
        c_b_blue  = "\x1b[1;34m";
        c_white   = "\x1b[1;37m";
    }
}

/* ── Decorative helpers ───────────────────────── */

/**
 * @brief Print a horizontal rule.
 */
static void hr(void)
{
    printf("  %s", c_dim);
    for (int ii = 0; ii < 68; ii++)
    {
        putchar('-');
    }
    printf("%s\n", c_reset);
}

/**
 * @brief Print a section heading.
 *
 * @param icon  Unicode icon string (e.g. a star emoji).
 * @param title Section title text.
 */
static void section(
    const char *icon,
    const char *title)
{
    printf("\n  %s%s  %s%s\n", c_cyan, icon, title, c_reset);
    hr();
}

/* ── Section printers ─────────────────────────── */

static void print_banner(void)
{
    printf("\n");
    printf("  %s", c_b_blue);
    printf("  __  __ _____ _____ ____  \n");
    printf("  |  \\/  | ____|  ___/ ___| \n");
    printf("  | |\\/| |  _| | |_  \\___ \\ \n");
    printf("  | |  | | |___|  _|  ___) |\n");
    printf("  |_|  |_|_____|_|   |____/ \n");
    printf("  %s\n", c_reset);
    printf("  %sMatched Electric Field Spectrograph%s\n",
           c_white, c_reset);
    printf("  %sSimulation & Analysis Toolkit%s\n",
           c_grey, c_reset);
    printf("\n");
}

static void print_overview(void)
{
    section("\xe2\x9c\xa8", "OVERVIEW");

    printf("  MEFS is an end-to-end numerical "
           "toolkit for modeling\n");
    printf("  high-contrast, high-resolution "
           "exoplanet spectroscopy\n");
    printf("  instruments (HWO / ELT class).\n\n");

    printf("  It uses a %stwo-stage starlight "
           "suppression%s strategy:\n",
           c_bold, c_reset);
    printf("    %s1.%s %sCoronagraph%s"
           "   ~10^7 contrast suppression\n",
           c_yellow, c_reset, c_b_green, c_reset);
    printf("    %s2.%s %sPhotonics%s"
           "     Lantern and/or chip for\n"
           "                     "
           "nulling & modal decomposition\n\n",
           c_yellow, c_reset, c_b_green, c_reset);
}

static void print_executables(void)
{
    section("\xf0\x9f\x94\xa7", "EXECUTABLES");

    /* mefs-create-scene */
    printf("  %s\xe2\x96\xb6%s  %smefs-create-scene%s\n",
           c_yellow, c_reset, c_b_green, c_reset);
    printf("     Astronomical scene generator.\n");
    printf("     Creates point-source lists for"
           " planets, zodiacal\n");
    printf("     dust, and exozodiacal disks.\n\n");

    /* mefs-psf */
    printf("  %s\xe2\x96\xb6%s  %smefs-psf%s\n",
           c_yellow, c_reset, c_b_green, c_reset);
    printf("     Core optical simulation engine.\n");
    printf("     Fourier propagation, coronagraph,"
           " IFS lenslet\n");
    printf("     array, and MEFS matched-filter"
           " projection.\n\n");

    /* mefs-score */
    printf("  %s\xe2\x96\xb6%s  %smefs-score%s\n",
           c_yellow, c_reset, c_b_green, c_reset);
    printf("     Performance grading tool.\n");
    printf("     Evaluates SNR, contrast curves,"
           " and detection\n");
    printf("     metrics from IFS simulation"
           " outputs.\n\n");
}

static void print_workflow(void)
{
    section("\xf0\x9f\x9a\x80", "TYPICAL WORKFLOW");

    printf("  %sStep 1%s  Generate scenes\n",
           c_b_cyan, c_reset);
    printf("    %s$%s %smefs-create-scene"
           " --benchmark%s\n\n",
           c_grey, c_reset, c_green, c_reset);

    printf("  %sStep 2%s  Run PSF simulation\n",
           c_b_cyan, c_reset);
    printf("    %s$%s %smefs-psf --scene"
           " scene.planet.txt --mefs_x 1.0%s\n\n",
           c_grey, c_reset, c_green, c_reset);

    printf("  %sStep 3%s  Grade performance\n",
           c_b_cyan, c_reset);
    printf("    %s$%s %smefs-score --scene"
           " scene.planet.txt%s\n\n",
           c_grey, c_reset, c_green, c_reset);
}

static void print_options_reference(void)
{
    section("\xf0\x9f\x93\x8b",
            "QUICK OPTION REFERENCE");

    /* ── mefs-create-scene ── */
    printf("  %smefs-create-scene%s %s[OPTIONS]%s\n",
           c_b_green, c_reset, c_magenta, c_reset);

    printf("    %s--type%s %-8s "
           "%splanet | disk | exozodi%s\n",
           c_green, c_reset, "<type>",
           c_grey, c_reset);
    printf("    %s-x%s / %s-y%s %-6s "
           "%sPosition (lam/D)%s\n",
           c_green, c_reset, c_green, c_reset,
           "<val>", c_grey, c_reset);
    printf("    %s-r%s %-12s %sDisk radius (lam/D)%s\n",
           c_green, c_reset, "<val>",
           c_grey, c_reset);
    printf("    %s-f%s %-12s %sTotal flux%s\n",
           c_green, c_reset, "<val>",
           c_grey, c_reset);
    printf("    %s-o%s %-12s %sOutput scene file%s\n",
           c_green, c_reset, "<file>",
           c_grey, c_reset);
    printf("    %s-a%s %-12s %sAppend mode%s\n",
           c_green, c_reset, "",
           c_grey, c_reset);
    printf("    %s-h%s %-12s %sFull help%s\n\n",
           c_green, c_reset, "",
           c_grey, c_reset);

    /* ── mefs-psf ── */
    printf("  %smefs-psf%s %s[OPTIONS]%s\n",
           c_b_green, c_reset, c_magenta, c_reset);

    printf("    %s-i%s %-12s %sInput pupil FITS%s\n",
           c_green, c_reset, "<file>",
           c_grey, c_reset);
    printf("    %s-o%s %-12s %sOutput complex PSF%s\n",
           c_green, c_reset, "<file>",
           c_grey, c_reset);
    printf("    %s-n%s %-12s %sPupil grid size%s\n",
           c_green, c_reset, "<size>",
           c_grey, c_reset);
    printf("    %s-s%s %-12s %sEnable IFS mode%s\n",
           c_green, c_reset, "",
           c_grey, c_reset);
    printf("    %s-c%s %-12s "
           "%sCoronagraph order (1-3)%s\n",
           c_green, c_reset, "<ord>",
           c_grey, c_reset);
    printf("    %s-x%s / %s-y%s %-6s "
           "%sSource offset (lam/D)%s\n",
           c_green, c_reset, c_green, c_reset,
           "<val>", c_grey, c_reset);
    printf("    %s-S%s %-12s %sScene definition file%s\n",
           c_green, c_reset, "<file>",
           c_grey, c_reset);
    printf("    %s-M%s %-12s %sEnable MEFS projection%s\n",
           c_green, c_reset, "",
           c_grey, c_reset);
    printf("    %s-h%s %-12s %sFull help%s\n\n",
           c_green, c_reset, "",
           c_grey, c_reset);

    /* ── mefs-score ── */
    printf("  %smefs-score%s %s[OPTIONS]%s\n",
           c_b_green, c_reset, c_magenta, c_reset);

    printf("    %s-p%s %-12s %sPlanet flux FITS%s\n",
           c_green, c_reset, "<file>",
           c_grey, c_reset);
    printf("    %s-s%s %-12s "
           "%sScenesum FITS%s\n",
           c_green, c_reset, "<file>",
           c_grey, c_reset);
    printf("    %s-o%s %-12s %sOutput ASCII file%s\n",
           c_green, c_reset, "<file>",
           c_grey, c_reset);
    printf("    %s-h%s %-12s %sFull help%s\n\n",
           c_green, c_reset, "",
           c_grey, c_reset);
}

static void print_build(void)
{
    section("\xf0\x9f\x8f\x97\xef\xb8\x8f ",
            "BUILDING");

    printf("  %smkdir -p build && cd build%s\n",
           c_green, c_reset);
    printf("  %scmake ..%s\n",
           c_green, c_reset);
    printf("  %smake -j$(nproc)%s\n\n",
           c_green, c_reset);

    printf("  %sRequires:%s CMake >= 3.10,"
           " CFITSIO, FFTW3, OpenMP\n\n",
           c_bold, c_reset);
}

static void print_links(void)
{
    section("\xf0\x9f\x94\x97", "LINKS");

    printf("  %sDocs%s   "
           "https://oguyon.github.io/MEFS/\n",
           c_b_cyan, c_reset);
    printf("  %sRepo%s   "
           "https://github.com/oguyon/MEFS\n\n",
           c_b_cyan, c_reset);
}

static void print_color_status(void)
{
    int colors_on = (getenv("NO_COLOR") == NULL);

    printf("  %sColor:%s %s%s%s",
           c_bold, c_reset,
           c_yellow,
           colors_on ? "ENABLED" : "DISABLED",
           c_reset);

    if (colors_on)
    {
        printf("  %s(set NO_COLOR to disable)%s",
               c_grey, c_reset);
    }
    printf("\n\n");
}

/* ── Main ─────────────────────────────────────── */

int main(
    int   argc,
    char **argv)
{
    (void)argc;
    (void)argv;

    init_colors();

    print_banner();
    print_overview();
    print_executables();
    print_workflow();
    print_options_reference();
    print_build();
    print_links();
    print_color_status();

    return 0;
}
