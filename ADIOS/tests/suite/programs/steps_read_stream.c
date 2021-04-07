/* Staged write example code.
   Assumptions:
     - one output step fits into the memory of the staged writer.
       Actually, this means, even more memory is needed than the size of output.
       We need to read each variable while also buffering all of them for output.
     - output steps contain the same variable set (no changes in variables)
     - attributes are the same for all steps (will write only once here)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include "adios_read.h"
#include "adios_error.h"

static const int max_read_buffer_size  = 1024*1024*1024;
static const int max_write_buffer_size = 1024*1024*1024;

static int timeout_sec = 30; // will stop if no data found for this time (-1: never stop)
enum ADIOS_READ_METHOD read_method = ADIOS_READ_METHOD_BP;


// Global variables
int         rank, numproc;
MPI_Comm    comm; 
ADIOS_FILE *f;      // stream for reading
int64_t    fh;     // ADIOS output file handle
int64_t     gh;     // ADIOS group for output definitions
uint64_t    write_total; // data size written by one processor
uint64_t    largest_block; // the largest variable block one process reads
char     ** group_namelist; // name of ADIOS group
char       *readbuf; // read buffer
int         decomp_values[10];


int process_metadata();
int read_records();


int main (int argc, char ** argv) 
{
    int         err;
    int         steps = 0, curr_step;
    int         retval = 0;

    MPI_Init (&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &numproc);

    
    char infilename[] = "steps.bp";

    err = adios_read_init_method(read_method, comm, 
                                 "max_chunk_size=100; "
                                 "app_id =32767; \n"
                                 "verbose= 3;"
                                 "poll_interval  =  100;"
                                );

    if (!err) {
        printf ("%s\n", adios_errmsg());
    }


    printf ("Waiting to open stream %s...\n", infilename);
    f = adios_read_open_stream (infilename, read_method, comm,
            ADIOS_LOCKMODE_ALL, timeout_sec);
    if (adios_errno == err_file_not_found) 
    {
        printf ("rank %d: Stream not found after waiting %d seconds: %s\n", 
               rank, timeout_sec, adios_errmsg());
        retval = adios_errno;
    } 
    else if (adios_errno == err_end_of_stream) 
    {
        printf ("rank %d: Stream terminated before open. %s\n", rank, adios_errmsg());
        retval = adios_errno;
    } 
    else if (f == NULL) {
        printf ("rank %d: Error at opening stream: %s\n", rank, adios_errmsg());
        retval = adios_errno;
    } 
    else 
    {
        // read data here... 
        while(1) {
            steps++; // start counting from 1

            /*
            printf ("File info:\n");
            printf ("  current step:   %d\n", f->current_step);
            printf ("  last step:      %d\n", f->last_step);
            printf ("  # of variables: %d:\n", f->nvars);
            */

            if (steps==1)
                retval = process_metadata();

            retval = read_records();

            // advance to 1) next available step with 2) blocking wait 
            curr_step = f->current_step; // save for final bye print
            adios_advance_step (f, 0, 0.0);

            if (adios_errno == err_end_of_stream)
            {
                break; // quit while loop
            }
            else if (adios_errno == err_step_notready)
            {
                printf ("rank %d: No new step arrived within the timeout. Quit. %s\n",
                        rank, adios_errmsg());
                break; // quit while loop
            }
            else if (f->current_step != curr_step+1)
            {
                // we missed some steps
                printf ("rank %d: WARNING: steps %d..%d were missed when advancing.\n",
                        rank, curr_step+1, f->current_step-1);
            }


        }
        adios_read_close (f);
    } 

    adios_read_finalize_method (read_method);
    MPI_Finalize ();

    return retval;
}


typedef struct {
    ADIOS_VARINFO * v;
    uint64_t        start[10];
    uint64_t        count[10];
    uint64_t        writesize; // size of subset this process writes, 0: do not write
} VarInfo;

VarInfo * varinfo;

int NX, Width, nblocks;

int process_metadata()
{
    int retval = 0;
    int i, j;
    char gdims[256], ldims[256], offs[256];
    uint64_t sum_count;
    ADIOS_VARINFO *v; // shortcut pointer

    /* First step processing */

    printf ("Get info on variable Width\n"); 
    v = adios_inq_var (f, "Width");
    if (v == NULL) {
        printf ("rank %d: ERROR: Variable %s inquiry failed: %s\n", 
                rank, "Width", adios_errmsg());
        return 1;
    }
    Width = *(int *)v->value;
    adios_free_varinfo (v);
    printf ("rank %d: Width = %d\n", rank, Width); 

    printf ("Get info on variable NX\n"); 
    v = adios_inq_var (f, "NX");
    if (v == NULL) {
        printf ("rank %d: ERROR: Variable %s inquiry failed: %s\n", 
                rank, "NX", adios_errmsg());
        return 1;
    }
    NX = *(int *)v->value;
    adios_free_varinfo (v);
    printf ("rank %d: NX = %d\n", rank, NX); 

    printf ("Get info on variable record\n"); 
    v = adios_inq_var (f, "record");
    if (v == NULL) {
        printf ("rank %d: ERROR: Variable %s inquiry failed: %s\n", 
                rank, "record", adios_errmsg());
        return 1;
    }
    nblocks = v->nblocks[0];
    printf ("rank %d: record dims = %llu * %llu \n", rank, v->dims[0], v->dims[1]); 
    adios_free_varinfo (v);
    printf ("rank %d: nblocks = %d\n", rank, nblocks); 

    return retval;
}

int read_records()
{
    int retval = 0;
    int i,j;

    int N = nblocks/numproc;
    int startidx=N*rank;

    char *text;
    text = malloc (Width*NX+1);
    text[Width*NX] = 0;
    
    for (i=0; i<N; i++) 
    {
        memset (text, '+', Width*NX);
        // read one block of records
        //printf ("rank %d: Read block %d\n", rank, startidx+i); 
        ADIOS_SELECTION *sel = adios_selection_writeblock (startidx+i);
        adios_schedule_read (f, sel, "record", 0, 1, text);
        adios_perform_reads (f, 1);   
        adios_selection_delete (sel);

        printf ("block %2d = ", startidx+i);
        for (j=0; j<NX; j++) 
            printf ("[%s]", text+j*Width);
        printf ("\n");
    }

    return retval;
}



