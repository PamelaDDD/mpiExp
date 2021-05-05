#include <iostream>
#include <mpi.h>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <cmath>
#define INT_MAX 2147483647
using namespace std;

/**
 * @param binary_file: the binary file on Disk.
 * @param vector: the values from binary file will store in the vector.
 * @return EXIT_SUCCESS or EXIT_FAILURE
 * */
/*int fill_vector_from_binary_file(int **data,
                                 char *binary_file,
                                 long rank,
                                 int count_processes,
                                 unsigned long &data_size)
{
    ifstream bin_file(binary_file, ios::in | ios::binary);
    bin_file.seekg(0, ios::end);

    long long int count_all_bytes = bin_file.tellg();
    long long int total_size;
    total_size = count_all_bytes / sizeof(int);
    if (total_size >= pow(2,32)){
        total_size -= 4096;
        printf("total_size has reduce, now is %lld\n",total_size);
    }
    data_size                      = (total_size / count_processes);
    if(data_size < 1) {
        cout << "The amount of processes is higher than the amount of values. In this case not every process will become values." << endl;
        return EXIT_FAILURE;
    }
    else if( (count_all_bytes / sizeof(int)) % count_processes != 0  ) {
        printf("the data num : %lld\n",count_all_bytes / sizeof(int));
        cout << "The amount of values have to be even on every process. Otherwise the sorting would be incorrect." << endl;
        return EXIT_FAILURE;
    }
    
    *data = new int[data_size];

    bin_file.seekg(rank * sizeof(int) * data_size, ios::beg);
    bin_file.read(reinterpret_cast<char *>(*data), sizeof(int) * data_size);
    printf("read successed! data size =%ld\n",data_size);

    return EXIT_SUCCESS;
}
*/

int fill_vector_from_binary_file(float **data,
                                 char *binary_file,
                                 long rank,
                                 int count_processes,
                                 unsigned long &data_size)
{
    ifstream bin_file(binary_file, ios::in | ios::binary);
    bin_file.seekg(0, ios::end);

    const long int count_all_bytes = bin_file.tellg();
    long int total_size = count_all_bytes / sizeof(int);
    if ((total_size) > pow(2,32)|| total_size == pow(2,32)){
        total_size -= 4096;
    }
    data_size                      = total_size / count_processes;

    if(data_size < 1) {
        cout << "The amount of processes is higher than the amount of values. In this case not every process will become values." << endl;
        return EXIT_FAILURE;
    }
    else if( (count_all_bytes / sizeof(int)) % count_processes != 0  ) {
        cout << "The amount of values have to be even on every process. Otherwise the sorting would be incorrect." << endl;
        return EXIT_FAILURE;
    }

    *data = new float[data_size];

    bin_file.seekg(rank * sizeof(float) * data_size, ios::beg);
    bin_file.read(reinterpret_cast<char *>(*data), sizeof(float) * data_size);

    return EXIT_SUCCESS;
}


void print(float *data, int rank, unsigned long data_size)
{
    cout << "rank " << rank << " : ";
    for(unsigned long int i=0; i<data_size; i++)
        cout << data[i] << " ";
    cout << endl;
}

int findPartner(int phase, int rank) {
    int partner;

    /* if it's an even phase */
    if (phase % 2 == 0) {
        /* if we are an even process */
        if (rank % 2 == 0) {
            partner = rank + 1;
        } else {
            partner = rank - 1;
        }
    } else {
        /* it's an odd phase - do the opposite */
        if (rank % 2 == 0) {
            partner = rank - 1;
        } else {
            partner = rank + 1;
        }
    }
    return partner;
}

int compare (const void * a, const void * b) {
    return  *(float*)a > *(float*)b ? 1 : 0 ;
}

void parallel_sort(float *data, int rank, int count_processes, unsigned long data_size)
{
    const unsigned long concat_data_size = data_size * 2;
    //printf("rank:%d,concat_data_size %ld\n",rank,concat_data_size);
    auto *other      = new int[data_size];
    auto *concatData = new int[concat_data_size];

    for (int i=0; i<count_processes; i++)
    {
        int partner = findPartner(i, rank);
        if (partner < 0 || partner >= count_processes)
          continue;

        if (rank % 2 == 0) {

            MPI_Send(data, (int) data_size,  MPI_FLOAT, partner, 0, MPI_COMM_WORLD);
            MPI_Recv(other, (int) data_size,  MPI_FLOAT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(other, (int) data_size,  MPI_FLOAT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(data, (int) data_size,  MPI_FLOAT, partner, 0, MPI_COMM_WORLD);
        }

        merge(data,  data  + data_size,
              other, other + data_size,
              concatData);
        qsort(concatData, concat_data_size, sizeof(float), compare);

        auto posHalfConcatData = concatData + data_size;
        if (rank < partner)
            copy(concatData, posHalfConcatData, data);
        else
            copy(posHalfConcatData, concatData + concat_data_size, data);
    }
}

/**
 * Compile:      mpic++ OddEvenSort.cpp -o OddEvenSort -std=gnu++0x
 * Example-Call: mpirun -np 4 ./OddEvenSort "<numbers_file.bin>" <y>
 * <y> output on console
 * */
int main(int argCount, char** argValues)
{
    int rank, count_processes;
    double start, end;
    MPI_Init(&argCount, &argValues);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    ifstream myFile(argValues[1], ios::in | ios::binary);
    /*myFile.seekg(0, ios::end);
    long long int count_all_bytes = myFile.tellg();
    unsigned long array_size           = count_all_bytes / sizeof(int);
    int *array               = new int[array_size];
    printf("array new success!\n");

    myFile.seekg(0, ios::beg);
    myFile.read(reinterpret_cast<char *>(array), sizeof(int) * array_size);
    
    clock_t qstart,qend;
    qstart = clock();

    
    qsort(array,array_size,sizeof(int),compare);
    printf("qsort success!\n");
    qend = clock();
    double qruntime;
    qruntime = (double)(qend - qstart) / CLOCKS_PER_SEC;
    */
    float *data;
    unsigned long data_size;
    printf("start reading from binary file\n");
    int status = fill_vector_from_binary_file(&data, argValues[1], rank, count_processes, data_size);

    if (status == EXIT_FAILURE) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    start = MPI_Wtime();
    parallel_sort(data, rank, count_processes, data_size);

    if(argCount > 2 && strcmp(argValues[2], "y") == 0) {
        print(data, rank, data_size);
    }
    end = MPI_Wtime();
    MPI_Finalize();
    if (rank == 0) {
        printf("Runtime of mpi = %f\n",end-start);
    }
    return EXIT_SUCCESS;
}
