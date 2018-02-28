#include <mpi.h>
#include <iostream>
#include <array>

using namespace std;

int main(int argc, char** argv){
	
	MPI_Init(&argc, &argv);

	// MPI
	int my_rank, num_procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// parameter, matrix size m×m
	int matrix_size = stoi(argv[1]);
	
	// process per dimension
	int dims[2] = {0, 0};
	MPI_Dims_create(num_procs, 2, dims);
	int procs_y = dims[0];
	int procs_x = dims[1];

	// new communicator with cartesian topology
	MPI_Comm comm_cart;
	int periods[2] = {0, 0};
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_cart);

	// determines process coords in Cartesian topology given rank
	int coords[2] = {0, 0};
	MPI_Cart_coords(comm_cart, my_rank, 2, coords);

	// divide number of column/rows through number of process on x/y axis
	auto div_x = div(matrix_size, procs_x);
	auto div_y = div(matrix_size, procs_y);
	int cols = div_x.quot;
	int rows = div_y.quot;

	// add remainder to processes
	if(coords[1] < div_x.rem)
		++cols;
	if(coords[0] < div_y.rem)
		++rows;

	// get neighbours
	int neighbour_left;
	int neighbour_right;
	int neighbour_up;
	int neighbour_down;

	MPI_Cart_shift(comm_cart, 1, 1, &neighbour_left, &neighbour_right);
	MPI_Cart_shift(comm_cart, 0, 1, &neighbour_up, &neighbour_down);

	int offset_left, offset_right, offset_up, offset_down;
	offset_down=offset_up=offset_left=offset_right=0;
	
	// add space for HALO		
	if(neighbour_left != MPI_PROC_NULL){
		++cols;
		offset_left=1;
	}
	if(neighbour_right != MPI_PROC_NULL){
		++cols;
		offset_right=1;
	}
	if(neighbour_up != MPI_PROC_NULL){
		++rows; 
		offset_up=1;
	}
 	if(neighbour_down != MPI_PROC_NULL){
		++rows; 
		offset_down=1;
 	}

	// data container
	float data[rows][cols];
	for(int i=0; i<rows; ++i)
		for(int j=0; j<cols;++j)
			data[i][j] = 0;
	// fill data
	for(int i=0+offset_up; i<rows-offset_down; ++i)
		for(int j=0+offset_left; j<cols-offset_right;++j)
			data[i][j] = 0.2+my_rank;

	// send data to process 0
	// if(my_rank == 0){
	// 	MPI_Request handler;
	// 	for(int rank=1; rank<num_procs; ++rank){
	// 		MPI_Irecv(tmp, rows*cols, MPI_FLOAT, rank, 0, MPI_COMM_WORLD, &handler);
	// 	}
	// 	MPI_Wait(&handler, MPI_STATUS_IGNORE);
	// } else {
	// 	MPI_SEND(data, rows*cols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
	// }






	// Output eve< rank
	for(int rank=0; rank<num_procs; ++rank){
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == rank){
			cout << "From process: " << my_rank << " coords: " << "(" << coords[0] << "," << coords[1] << ")" <<endl;
			for(int i=0; i < rows; ++i){
				for(int j=0; j < cols; ++j)
					cout << data[i][j] << "\t";
				cout << "\n";
			}
			cout << "\n" << endl;
		}	
	}

	MPI_Finalize();
	return 0;
}