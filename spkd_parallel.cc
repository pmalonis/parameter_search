#include <iostream>
#include <mpi.h>

int main(int argc, char** argv) {

    int status, tid, nthreads;
    MPI_File fh;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);

    char filename[] = "/lustre/beagle2/pmalonis/parameter_search/search.bin";
    status = MPI_File_open(MPI_COMM_WORLD, filename, 
                           MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    MPI_Offset filesize;
    status = MPI_File_get_size(fh, &filesize);
    
    if (tid == 0) {
        std::cout << "filesize: " << filesize << "\n";
    }
    int row_size = 9;
    int spike_idx = 6;
    int max_nspikes = 20;
    
    MPI_Finalize();

    return 0;
}
