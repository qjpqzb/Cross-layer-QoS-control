/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* Application that performs computation and communication to a user-given ratio.
   Then it performs output I/O with posix, parallel hdf5 and adios.
   It writes a global 2D array of N*NX x NY, where N is the number of processes and
   NX, NY are user-given sizes.
   
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
//#include <unistd.h>
//#include <fcntl.h>
#include <errno.h>
#include "mpi.h"
#include "output.h"
#include "timing.h"
#include "comp_comm.h"


// User parameters
int    nsteps;            // number of output steps
int    niterations;       // number of iterations between steps
int    ncomp;             // computation units in one iteration
int    ncomm;             // communication units in one iteration
int    nx;                // array size in x direction
int    ny;                // array size in y direction


void printUsage(char *prgname)
{
    printf("Usage: mpirun -np <N> %s <steps>  <it>  <ncomp>  <ncomm>  <nx>  <ny>\n"
           "  <steps>   number of steps (outputs)\n"
           "  <it>      number of iterations in each step (between outputs)\n"
           "\n"
           "  <ncomp>   computation units per iteration"
           "  <ncomm>   communication units per iteration"
           "\n"
           "  <nx> <ny> size of a 2D array generated by each process\n"
           "            The output 2D array will be <N>*<nx> times <ny>\n"
           "            e.g. for N=12 nx=4 ny=3, output will be 48x3\n"
        ,prgname);
}

int convert_arg_to_int (char **argv, int argpos, char *name, int *value)
{
    char *end;
    errno = 0;
    *value = strtol(argv[argpos], &end, 10);
    if (errno || (end != 0 && *end != '\0')) {
        printf ("ERROR: Invalid argument for %s: '%s'\n", name, argv[argpos]);
        printUsage(argv[0]);
        return 1;
    }
    return 0;
}

int processArgs(int argc, char ** argv)
{
    if (argc < 7) {
        printUsage (argv[0]);
        return 1;
    }

    if (convert_arg_to_int (argv, 1, "<steps>", &nsteps)) 
        return 1;
    if (convert_arg_to_int (argv, 2, "<it>", &niterations)) 
        return 1;
    if (convert_arg_to_int (argv, 3, "<ncomp>", &ncomp)) 
        return 1;
    if (convert_arg_to_int (argv, 4, "<ncomm>", &ncomm)) 
        return 1;
    if (convert_arg_to_int (argv, 5, "<nx>", &nx)) 
        return 1;
    if (convert_arg_to_int (argv, 6, "<ny>", &ny)) 
        return 1;

    return 0;
}

// other global variables
MPI_Comm  comm;
int    rank;
int    nproc;             // # of total procs
int    offs_x, offs_y;    // offset in x and y direction
int    gnx, gny;          // global array size
double * data;            // data array to do calc on it and output


void data_init()
{
    int i,j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            data[i*ny + j] = 1.0*rank;
        }
    }
}


#define MAX(a,b) (a<b ? b : a)

int main (int argc, char ** argv ) 
{
    int i;
    MPI_Init (&argc, &argv);
    MPI_Comm_dup (MPI_COMM_WORLD, &comm);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &nproc);
    comp_comm_init(comm);

    if (processArgs(argc, argv)) {
        return 1;
    }

    if (!rank) {
        printf ("Setup parameters:\n");
        printf ("  # of steps (outputs):     %d\n", nsteps);
        printf ("  # of iterations per step: %d\n", niterations);
        printf ("  # of computation units in each iteration:   %d\n", ncomp);
        printf ("  # of communication units in each iteration: %d\n", ncomm);
        printf ("  output size per process: %d x %d doubles = %lld bytes\n", 
                nx, ny, sizeof(double) * nx * (uint64_t) ny);
        printf ("  output size per step: %lld bytes\n", 
                nproc * sizeof(double) * nx * (uint64_t) ny);
    }

    //2D array with 1D decomposition
    offs_x = rank * nx;
    offs_y = 0;
    gnx = nproc * nx;
    gny = ny;

    data = (double*) malloc (sizeof(double) * nx * (size_t) ny);
    timing_alloc(nsteps);


    int bufsizeMB = nx*ny*sizeof(double)/1048576 + 1;
    output_init(comm, bufsizeMB);
    output_define(nx, ny, gnx, gny, offs_x, offs_y);

    int it, step, icomp, icomm;

    /* Warm up a bit */
    if (!rank) printf ("Warm up for 1 steps, %d iterations per step...\n", niterations);
    for (step=0; step < 1; step++) {
        for (it=0; it < niterations; it++) {
            for (icomp=0; icomp < ncomp; icomp++) {
                do_calc_unit (data, nx, ny);
            }
            for (icomm=0; icomm < ncomm; icomm++) {
                do_comm_unit (comm);
            }
        }
    }


    /* Do the steps with output now */
    data_init();
    if (!rank) printf ("Start running with I/O and measurements...\n");
    double Tcalc_it, Tcomm_it;
    double Truntime; //to print total time for the loop below (for overhead calculation)
    char filename[256];

    MPI_Barrier (comm);
    Truntime = MPI_Wtime();

    for (step=0; step < nsteps; step++) {
        if (!rank) printf ("Start step %d\n", step);
        Tcalc[step] = 0;
        Tcomm[step] = 0;
        for (it=0; it < niterations; it++) {
            // spend some time with computation
            Tcalc_it = MPI_Wtime();
            for (icomp=0; icomp < ncomp; icomp++) {
                do_calc_unit (data, nx, ny);
            }
            Tcalc_it = MPI_Wtime() - Tcalc_it;
            Tcalc[step] += Tcalc_it;

            // spend some time with communication
            Tcomm_it = MPI_Wtime();
            for (icomm=0; icomm < ncomm; icomm++) {
                do_comm_unit (comm);
            }
            Tcomm_it = MPI_Wtime() - Tcomm_it;
            Tcomm[step] += Tcomm_it;
        }
        // output per step
        snprintf (filename, sizeof(filename), "data%6.6d", step);
        if (!rank) printf ("Output to %s\n", filename);
        MPI_Barrier (comm);
        output_dump(filename, step, data);
    }

    MPI_Barrier (comm);
    Truntime = MPI_Wtime() - Truntime;

    if (!rank) printf ("Finalize...\n");
    MPI_Barrier (comm);
    output_finalize (rank);

    timing_report(nsteps, comm);
    double Truntime_max;
    MPI_Reduce (&Truntime, &Truntime_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    if (!rank) printf ("Total runtime of main loop: %9.3f\n", Truntime);
    free (data);
    timing_free();

    MPI_Barrier (comm);
    MPI_Finalize ();
    return 0;
}