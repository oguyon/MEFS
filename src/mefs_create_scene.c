#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

/* Option character codes for long options */
#define OPT_TYPE  2000
#define OPT_INNER 2001
#define OPT_OUTER 2002

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const char *c_reset   = "";
const char *c_bold    = "";
const char *c_cyan    = "";
const char *c_b_green = "";
const char *c_green   = "";
const char *c_magenta = "";
const char *c_yellow  = "";
const char *c_b_red   = "";
const char *c_grey    = "";

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

void print_help(
    const char *prog_name)
{
    printf("%sNAME%s\n", c_cyan, c_reset);
    printf("  %s%s%s - Generate point source scene definition files.\n\n",
           c_b_green, prog_name, c_reset);

    printf("%sUSAGE%s\n", c_cyan, c_reset);
    printf("  %s%s%s [OPTIONS]\n\n", c_b_green, prog_name, c_reset);

    printf("%sDESCRIPTION%s\n", c_cyan, c_reset);
    printf("  Generates ASCII point source lists (x y flux per line) suitable for input\n");
    printf("  to %smefs_psf%s scene mode. Can generate single planets, uniform circular\n",
           c_b_green, c_reset);
    printf("  disks, or exozodi disks. If run with no options, generates 4 default scenes.\n\n");

    printf("%sOPTIONS%s\n", c_cyan, c_reset);
    printf("  %s--type%s %s<planet|disk|exozodi>%s  Type of component to generate.\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-x%s %s<val>%s                      X position in l/D (planet; default: 0.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-y%s %s<val>%s                      Y position in l/D (planet; default: 0.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-r%s %s<val>%s                      Radius in l/D (uniform disk; default: 5.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-f%s %s<val>%s                      Total integrated flux (default: 1.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-s%s %s<val>%s                      Sampling resolution in l/D (default: 0.1).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s--inner%s %s<val>%s                 Inner radius in l/D (exozodi; default: 1.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s--outer%s %s<val>%s                 "
           "Outer radius in l/D (exozodi; default: 10.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-p%s %s<val>%s                      Exozodi power law index (default: -2.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-i, --inclination%s %s<val>%s       "
           "Inclination in degrees (exozodi; default: 60.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-o%s %s<file>%s                     Output file path (default: scene.txt).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-a%s                           Append to output file instead of overwriting.\n",
           c_green, c_reset);
    printf("  %s-h, --help%s                   Show this help message.\n\n",
           c_green, c_reset);

    printf("%sEXAMPLES%s\n", c_cyan, c_reset);
    printf("  1. Generate off-axis planet:\n");
    printf("     %s$%s %s%s --type planet -x 2.5 -y 0.0 -f 1e-5 -o scene.txt%s\n",
           c_grey, c_reset, c_b_green, prog_name, c_reset);
    printf("  2. Append uniform disk to existing scene:\n");
    printf("     %s$%s %s%s --type disk -r 3.0 -f 1e-3 -s 0.2 -a -o scene.txt%s\n",
           c_grey, c_reset, c_b_green, prog_name, c_reset);
    printf("  3. Generate power-law exozodi:\n");
    printf("     %s$%s %s%s --type exozodi --inner 1.0 --outer 5.0 "
           "-p -2.0 -f 1e-4 -o zodi.txt%s\n\n",
           c_grey, c_reset, c_b_green, prog_name, c_reset);

    printf("%sCOLOR MODE%s\n", c_cyan, c_reset);
    printf("  Colors: %s%s%s. Set %sNO_COLOR%s env var to disable "
           "(see https://no-color.org).\n",
           c_yellow, (getenv("NO_COLOR") == NULL) ? "ENABLED" : "DISABLED", c_reset,
           c_bold, c_reset);
}

/* Helper function to generate a scene file */
int generate_scene_file(
    const char *filename,
    const char *type,
    double      x,
    double      y,
    double      r,
    double      flux,
    double      sampling,
    double      inner_rad,
    double      outer_rad,
    double      power_law,
    double      inclination_deg)
{
    FILE *f = fopen(filename, "w");
    if (f == NULL)
    {
        fprintf(stderr, "%sERROR: Could not open file %s for writing.%s\n",
                c_b_red, filename, c_reset);
        return 1;
    }

    if (strcmp(type, "planet") == 0)
    {
        fprintf(f, "%f %f %.15e\n", x, y, flux);
        printf("Created %s (planet: x=%.2f, y=%.2f, flux=%.2e)\n", filename, x, y, flux);
    }
    else if (strcmp(type, "disk") == 0)
    {
        int count = 0;
        for (double dy = -r; dy <= r + 1e-9; dy += sampling)
        {
            for (double dx = -r; dx <= r + 1e-9; dx += sampling)
            {
                if (dx * dx + dy * dy <= r * r)
                {
                    count++;
                }
            }
        }
        if (count == 0)
        {
            fprintf(stderr, "%sERROR: No grid points found inside disk for %s.%s\n",
                    c_b_red, filename, c_reset);
            fclose(f);
            return 1;
        }
        double pt_flux = flux / count;
        for (double dy = -r; dy <= r + 1e-9; dy += sampling)
        {
            for (double dx = -r; dx <= r + 1e-9; dx += sampling)
            {
                if (dx * dx + dy * dy <= r * r)
                {
                    fprintf(f, "%f %f %.15e\n", dx, dy, pt_flux);
                }
            }
        }
        printf("Created %s (disk: radius=%.4f, total_flux=%.2e, points=%d)\n",
               filename, r, flux, count);
    }
    else if (strcmp(type, "exozodi") == 0)
    {
        double total_weight = 0.0;
        int count = 0;
        double inc_rad = inclination_deg * M_PI / 180.0;
        double cos_inc = cos(inc_rad);
        if (cos_inc < 1e-9) cos_inc = 1e-9;

        for (double dy = -outer_rad; dy <= outer_rad + 1e-9; dy += sampling)
        {
            for (double dx = -outer_rad; dx <= outer_rad + 1e-9; dx += sampling)
            {
                double dist = sqrt(dx * dx + (dy / cos_inc) * (dy / cos_inc));
                if (dist >= inner_rad && dist <= outer_rad)
                {
                    total_weight += pow(dist, power_law);
                    count++;
                }
            }
        }
        if (count == 0 || total_weight <= 0.0)
        {
            fprintf(stderr, "%sERROR: No grid points found in exozodi bounds for %s.%s\n",
                    c_b_red, filename, c_reset);
            fclose(f);
            return 1;
        }
        for (double dy = -outer_rad; dy <= outer_rad + 1e-9; dy += sampling)
        {
            for (double dx = -outer_rad; dx <= outer_rad + 1e-9; dx += sampling)
            {
                double dist = sqrt(dx * dx + (dy / cos_inc) * (dy / cos_inc));
                if (dist >= inner_rad && dist <= outer_rad)
                {
                    double w = pow(dist, power_law);
                    double pt_flux = flux * (w / total_weight);
                    fprintf(f, "%f %f %.15e\n", dx, dy, pt_flux);
                }
            }
        }
        printf("Created %s (exozodi: inner=%.2f, outer=%.2f, total_flux=%.2e, "
               "inc=%.1f deg, points=%d)\n",
               filename, inner_rad, outer_rad, flux, inclination_deg, count);
    }
    fclose(f);
    return 0;
}

int main(
    int   argc,
    char **argv)
{
    init_colors();

    if (argc == 1)
    {
        printf("%sNo options specified. Generating default benchmark scenes...%s\n",
               c_yellow, c_reset);
        int st = 0;
        st |= generate_scene_file("scene.star.txt", "disk", 0.0, 0.0, 0.0075,
                                  1e10, 0.001, 0.0, 0.0, 0.0, 0.0);
        st |= generate_scene_file("scene.zodi.txt", "disk", 0.0, 0.0, 8.0,
                                  200.0, 0.1, 0.0, 0.0, 0.0, 0.0);
        st |= generate_scene_file("scene.exozodi.txt", "exozodi", 0.0, 0.0, 0.0,
                                  200.0, 0.1, 0.15, 3.0, -2.0, 60.0);
        st |= generate_scene_file("scene.planet.txt", "planet", 1.5, 0.0, 0.0,
                                  1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        if (st == 0)
        {
            printf("%sAll default scenes successfully generated.%s\n",
                   c_green, c_reset);
        }
        return st;
    }

    char *type = NULL;
    double x = 0.0;
    double y = 0.0;
    double r = 5.0;
    double flux = 1.0;
    double sampling = 0.1;
    double inner_rad = 1.0;
    double outer_rad = 10.0;
    double power_law = -2.0;
    double inclination_deg = -1.0;
    char *out_file = "scene.txt";
    int append_mode = 0;

    static struct option long_options[] = {
        {"type", required_argument, 0, OPT_TYPE},
        {"inner", required_argument, 0, OPT_INNER},
        {"outer", required_argument, 0, OPT_OUTER},
        {"inclination", required_argument, 0, 'i'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "x:y:r:f:s:p:o:i:ah",
                             long_options, &option_index)) != -1)
    {
        switch (opt)
        {
            case OPT_TYPE:
                type = optarg;
                break;
            case OPT_INNER:
                inner_rad = atof(optarg);
                break;
            case OPT_OUTER:
                outer_rad = atof(optarg);
                break;
            case 'x':
                x = atof(optarg);
                break;
            case 'y':
                y = atof(optarg);
                break;
            case 'r':
                r = atof(optarg);
                break;
            case 'f':
                flux = atof(optarg);
                break;
            case 's':
                sampling = atof(optarg);
                break;
            case 'p':
                power_law = atof(optarg);
                break;
            case 'o':
                out_file = optarg;
                break;
            case 'i':
                inclination_deg = atof(optarg);
                break;
            case 'a':
                append_mode = 1;
                break;
            case 'h':
                print_help(argv[0]);
                return 0;
            default:
                print_help(argv[0]);
                return 1;
        }
    }

    if (type == NULL)
    {
        fprintf(stderr, "%sERROR: Option --type is required (planet, disk, or exozodi).%s\n\n",
                c_b_red, c_reset);
        print_help(argv[0]);
        return 1;
    }

    if (sampling <= 0.0)
    {
        fprintf(stderr, "%sERROR: Sampling resolution must be positive.%s\n", c_b_red, c_reset);
        return 1;
    }

    if (flux <= 0.0)
    {
        fprintf(stderr, "%sERROR: Total flux must be positive.%s\n", c_b_red, c_reset);
        return 1;
    }

    FILE *f = fopen(out_file, append_mode ? "a" : "w");
    if (f == NULL)
    {
        fprintf(stderr, "%sERROR: Could not open output file %s%s\n", c_b_red, out_file, c_reset);
        return 1;
    }

    if (strcmp(type, "planet") == 0)
    {
        fprintf(f, "%f %f %.15e\n", x, y, flux);
        printf("Planet scene component added to %s: at (%f, %f) with flux %.3e\n",
               out_file, x, y, flux);
    }
    else if (strcmp(type, "disk") == 0)
    {
        if (r <= 0.0)
        {
            fprintf(stderr, "%sERROR: Disk radius must be positive.%s\n", c_b_red, c_reset);
            fclose(f);
            return 1;
        }

        /* First count points in grid */
        int count = 0;
        for (double dy = -r; dy <= r + 1e-9; dy += sampling)
        {
            for (double dx = -r; dx <= r + 1e-9; dx += sampling)
            {
                if (dx * dx + dy * dy <= r * r)
                {
                    count++;
                }
            }
        }

        if (count == 0)
        {
            fprintf(stderr, "%sERROR: No grid points found inside disk with sampling %f.%s\n",
                    c_b_red, sampling, c_reset);
            fclose(f);
            return 1;
        }

        /* Write out each point source with equal share of flux */
        double pt_flux = flux / count;
        for (double dy = -r; dy <= r + 1e-9; dy += sampling)
        {
            for (double dx = -r; dx <= r + 1e-9; dx += sampling)
            {
                if (dx * dx + dy * dy <= r * r)
                {
                    fprintf(f, "%f %f %.15e\n", dx, dy, pt_flux);
                }
            }
        }
        printf("Uniform disk scene component added to %s: radius %f, points %d, total flux %.3e\n",
               out_file, r, count, flux);
    }
    else if (strcmp(type, "exozodi") == 0)
    {
        if (inner_rad < 0.0 || outer_rad <= inner_rad)
        {
            fprintf(stderr, "%sERROR: Exozodi bounds must satisfy 0 <= inner < outer.%s\n",
                    c_b_red, c_reset);
            fclose(f);
            return 1;
        }

        if (inclination_deg < 0.0)
        {
            inclination_deg = 60.0;
        }

        double inc_rad = inclination_deg * M_PI / 180.0;
        double cos_inc = cos(inc_rad);
        if (cos_inc < 1e-9) cos_inc = 1e-9;

        /* First loop to calculate total weight */
        double total_weight = 0.0;
        int count = 0;
        for (double dy = -outer_rad; dy <= outer_rad + 1e-9; dy += sampling)
        {
            for (double dx = -outer_rad; dx <= outer_rad + 1e-9; dx += sampling)
            {
                double dist = sqrt(dx * dx + (dy / cos_inc) * (dy / cos_inc));
                if (dist >= inner_rad && dist <= outer_rad)
                {
                    total_weight += pow(dist, power_law);
                    count++;
                }
            }
        }

        if (count == 0 || total_weight <= 0.0)
        {
            fprintf(stderr, "%sERROR: No grid points in exozodi bounds with sampling %f.%s\n",
                    c_b_red, sampling, c_reset);
            fclose(f);
            return 1;
        }

        /* Write out each point source with weight normalized to total flux */
        for (double dy = -outer_rad; dy <= outer_rad + 1e-9; dy += sampling)
        {
            for (double dx = -outer_rad; dx <= outer_rad + 1e-9; dx += sampling)
            {
                double dist = sqrt(dx * dx + (dy / cos_inc) * (dy / cos_inc));
                if (dist >= inner_rad && dist <= outer_rad)
                {
                    double w = pow(dist, power_law);
                    double pt_flux = flux * (w / total_weight);
                    fprintf(f, "%f %f %.15e\n", dx, dy, pt_flux);
                }
            }
        }
        printf("Exozodi disk scene added to %s: inner %f, outer %f, inc %.1f deg, "
               "points %d, total flux %.3e\n",
               out_file, inner_rad, outer_rad, inclination_deg, count, flux);
    }
    else
    {
        fprintf(stderr, "%sERROR: Unknown component type: %s%s\n", c_b_red, type, c_reset);
        fclose(f);
        return 1;
    }

    fclose(f);
    return 0;
}
