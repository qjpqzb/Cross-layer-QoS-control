/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 *
 * This code is using the generic read API, which can read in
 * arbitrary slices of an array and thus we can read in an array
 * on arbitrary number of processes (provided our code is smart 
 * enough to do the domain decomposition).
 *
 * Run this example after adios_global, which generates 
 * adios_global.bp. Run this example on equal or less 
 * number of processes since we decompose only on one 
 * dimension of the global array here. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"
#include "adios.h"
#include <time.h>
#define TEST 

int main (int argc, char ** argv) 
{
    int         rank, size, i, j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;
    void * mesh = NULL;
    double * data = NULL, * grad = NULL, * R = NULL, * Z = NULL;
    uint64_t start[2], count[2];

    MPI_Init (&argc, &argv);
    //MPI_Init(NULL, NULL);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("From processor %s, rank %d out of %d processors\n",
           processor_name, rank, size);
    int step;
    double iotime = 0.0;
#ifdef TEST 
    for (step = 0; step < 1; step++)
    {
        char fname[64];

        //sprintf (fname, "xgc.3d.080%02d.bp", step);
        sprintf (fname, "larger_data/current_data/astro2D.bp");
        //sprintf (fname, "larger_data/data32X/xgc.3d.08000.bp");
        //sprintf (fname, "/nvme/sirius_storage/examples/C/global-array/larger_data/128X_data/xgc.3d.08000.bp");
        printf("dodosada=%f\n",1e10);
        adios_read_init_method (method, comm, "verbose=3");

        ADIOS_FILE * f = adios_read_open_file (fname, method, comm);
        printf("Finish read file\n"); 
        if (f == NULL)
        {
            printf ("%s\n", adios_errmsg());
            return -1;
        }

        ADIOS_VARINFO * v = adios_inq_var (f, "dpot");
        printf ("%lu, %lu\n", v->dims[0], v->dims[1]);
        //printf("size=%d\n",size);
        //int intv = v->dims[0] / size;
        //int remain = v->dims[0] % size;
	int intv = v->dims[1]/size;
	int remain = v->dims[1] % size;
       
        printf("intv=%d\n",intv);
        printf("rank=%d\n",rank);
/*
        if (rank == size -1) count[0] = intv + remain;
	else count[0] = intv;
        start[0] = rank * intv;
        //count[0] = intv;
        start[1] = 0;
        count[1] = 1;
*/
	start[0] = 0;
	count[0] = v->dims[0];
	start[1] = intv * rank;
	if (rank == size -1) count[1] = intv + remain;
	else count[1] = intv;	
        
        //start[0] = rank*intv;
        //count[0] = intv;
        //start[1] = 0;
        //count[1] = 1;
        //printf("v->dims[0]=%d\n",v->dims[0]);
        printf("count[0] * count[1]=%d\n",count[0] * count[1]); 
        data = malloc (count[0] * count[1] * sizeof (double));
        grad = malloc (count[0] * count[1] * sizeof (double));
        R = malloc (count[0] * count[1] * sizeof (double));
        Z = malloc (count[0] * count[1] * sizeof (double));
        if (data == NULL || grad == NULL || R == NULL || Z == NULL)
        {
            fprintf (stderr, "malloc failed.\n");
            return -1;
        }

        /* Read a subset of the temperature array */
        printf("start[0]=%d, count[0]=%d, start[1]=%d, count[1]=%d\n", start[0], count[0], start[1], count[1]);
        sel = adios_selection_boundingbox (v->ndim, start, count);
        adios_schedule_read (f, sel, "dpot", 0, 1, data);
        adios_perform_reads (f, 1);
        /*
        for (i = 0; i < count[0] * count[1]; i++) {
            printf (" %e", * ((double *)data + i ));
        }
        */
//    printf ("\n");

        adios_read_close (f);
        MPI_Barrier (comm);
        adios_read_finalize_method (method);

        // read mesh
        adios_read_init_method (method, comm, "verbose=3");
        //ADIOS_FILE * fmesh = adios_read_open_file ("/nvme/sirius_storage/examples/C/global-array/larger_data/128X_data/xgc.mesh.bp", method, comm);
        ADIOS_FILE * fmesh = adios_read_open_file ("larger_data/current_data/astro2D.mesh.bp", method, comm);
        //ADIOS_FILE * fmesh = adios_read_open_file ("larger_data/data32X/xgc.mesh.bp", method, comm);
        
//        ADIOS_FILE * fmesh = adios_read_open_file ("data1/dpot.3d.08000.bp", method, comm);
        
        if (fmesh == NULL)
        {
            printf ("%s\n", adios_errmsg());
            return -1;
        }

        ADIOS_VARINFO * conn = adios_inq_var (fmesh, "/cell_set[0]/node_connect_list");
//        ADIOS_VARINFO * conn = adios_inq_var (fmesh, "mesh/L1");
        printf ("conn: %lu, %lu\n", conn->dims[0], conn->dims[1]);


        start[0] = 0;
        count[0] = conn->dims[0];

        start[1] = 0;
        count[1] = conn->dims[1];


        mesh = malloc (count[0] * count[1] * sizeof (int));
        if (mesh == NULL)
        {
            fprintf (stderr, "malloc failed.\n");
            return -1;
        }



        /* Read a subset of the temperature array */
//        sel = adios_selection_boundingbox (conn->ndim, start, count);
        adios_schedule_read (fmesh, sel, "/cell_set[0]/node_connect_list", 0, 1, mesh);
//        adios_schedule_read (fmesh, sel, "mesh/L1", 0, 1, mesh);
        //start[0] = rank * intv;
        
        start[0] = 0;
        start[1] = 0;
        //if (rank == size -1) count[0] = intv + remain;
        //else count[0] = intv;
        count[0] = v->dims[0];
        count[1] = 1;
	printf("start[0]=%d, count[0]=%d, start[1]=%d, count[1]=%d\n", start[0], count[0], start[1], count[1]);
        sel = adios_selection_boundingbox (v->ndim, start, count);
        adios_schedule_read (fmesh, sel, "/coordinates/values", 0, 1, R);
//        adios_schedule_read (fmesh, sel, "R/L1", 0, 1, R);
printf("v->ndim = %d\n", v->ndim);
        //start[0] = rank * intv;
	    start[0] = 0; 
        start[1] = 1;
        //if (rank == size -1) count[0] = intv + remain;
        //else count[0] = intv;
        count[0] = v->dims[0];
        count[1] = 1;
        printf("start[0]=%d, count[0]=%d, start[1]=%d, count[1]=%d\n", start[0], count[0], start[1], count[1]);
        sel = adios_selection_boundingbox (v->ndim, start, count);
        adios_schedule_read (fmesh, sel, "/coordinates/values", 0, 1, Z);
//        adios_schedule_read (fmesh, sel, "Z/L1", 0, 1, Z);
        adios_perform_reads (fmesh, 1);
/*
    for (i = 0; i < count[0]; i++) {
        int n1 = * ((int *) mesh + i * 3);
        int n2 = * ((int *) mesh + i * 3 + 1);
        int n3 = * ((int *) mesh + i * 3 + 2);

        grad[n1] = grad[n2] = grad[n3] = 
              (double)(abs (data[n1] - data[n2]) + abs(data[n1] - data[n3]) + abs(data[n2] - data[n3])) / 3;
    }
*/      printf("read finish\n"); 
        adios_read_close (fmesh);
        MPI_Barrier (comm);
        adios_read_finalize_method (method);
        //double z_min=Z[0], z_max=Z[0]; 
        //for (int i=0; i< count[0] * count[1];i++ )
	//{
	//	printf("R[%d],Z[%d]=(%f,%f)\n",i,i,R[i],Z[i]);
        //        if (Z[i] > z_max) z_max = Z[i];
        //        else if (Z[i] < z_min) z_min = Z[i]; 
	//}
        //printf("z_max = %f, z_min = %f\n",z_max, z_min);

        int NX = count[0], GX = v->dims[0], OX = 0;
        int MY = conn->dims[0], MX = conn->dims[1];
        //printf("count[0]=%d\n",count[0]);
        
        uint64_t    adios_groupsize, adios_totalsize;
        int64_t     adios_handle;

        MPI_Barrier (MPI_COMM_WORLD);
        double start_io_time = MPI_Wtime ();
        adios_init ("test_xgc.xml", comm);
        
        sprintf (fname, "astro2D.bp");
        //sprintf (fname, "dpot.3d.080%02d_1.bp", step);

        adios_open (&adios_handle, "field", fname, "w", comm);
        adios_groupsize = 6 * 4 + 3 * 8 * NX + 4 * MX * MY;
        adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
        printf("start write in test_xgc\n");
        printf("NX\n"); 
        adios_write (adios_handle, "NX", &NX);
        printf("GX\n");
        adios_write (adios_handle, "GX", &GX);
        printf("OX\n");
        adios_write (adios_handle, "OX", &OX);
	printf("MX\n");
        adios_write (adios_handle, "MX", &MX);
	printf("MY\n");
        adios_write (adios_handle, "MY", &MY);
        
	printf("mesh\n");
        adios_write (adios_handle, "mesh", mesh);
        
	printf("R\n");
        adios_write (adios_handle, "R", R);
	printf("Z\n");
        adios_write (adios_handle, "Z", Z);
	printf("dpot\n");
        adios_write (adios_handle, "dpot", data);       
	

        adios_close (adios_handle);

        free (data);
        free (grad);
        free (mesh);
        adios_finalize (rank);

        MPI_Barrier (MPI_COMM_WORLD);
        double end_io_time = MPI_Wtime ();
    
        iotime += end_io_time - start_io_time;
        printf ("io time = %f\n", iotime);
    }
    //if (rank == 0) printf ("io time = %f\n", iotime);
#endif   
    MPI_Finalize ();
    return 0;
}
