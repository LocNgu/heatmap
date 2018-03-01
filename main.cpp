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

	int offset_left=0;
	int offset_right=0;
	int offset_up=0;
	int offset_down=0;
	int cols_halo = cols;
	int rows_halo = rows;

	// add space for HALO		
	if(neighbour_left != MPI_PROC_NULL){
		++cols_halo;
		offset_left=1;
	}
	if(neighbour_right != MPI_PROC_NULL){
		++cols_halo;
		offset_right=1;
	}
	if(neighbour_up != MPI_PROC_NULL){
		++rows_halo; 
		offset_up=1;
	}
 	if(neighbour_down != MPI_PROC_NULL){
		++rows_halo; 
		offset_down=1;
 	}

	// data container
	float data[rows_halo][cols_halo];
	// initialize container with 0
	for(int i=0; i<rows_halo; ++i)
		for(int j=0; j<cols_halo;++j)
			data[i][j] = 0;

	// fill with data
	for(int i=0+offset_up; i<rows_halo-offset_down; ++i)
		for(int j=0+offset_left; j<cols_halo-offset_right;++j)
			data[i][j] = 0.2;

	// heat soruce
	// border
	if(neighbour_up == MPI_PROC_NULL)
		for(int i=0; i < cols; ++i)
			data[0][i] = 1;
	if(neighbour_down == MPI_PROC_NULL)
		for(int i=0; i < cols; ++i)
			data[rows_halo-1][i] = 1;
	if(neighbour_left == MPI_PROC_NULL)
		for(int i=0; i < rows; ++i)
			data[i][0] = 1;
	if(neighbour_right == MPI_PROC_NULL)
		for(int i=0; i < rows; ++i)
			data[i][cols_halo-1] = 1;


	
	/////////////////////////
	// send data to neighbour
	// halo communication
	/////////////////////////
	
	// send data west -> east	
	MPI_Datatype new_type;
		// int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
	MPI_Type_vector(rows_halo, 1, cols_halo, MPI_FLOAT, &new_type);
	MPI_Type_commit(&new_type);

	if(neighbour_left == MPI_PROC_NULL && neighbour_right != MPI_PROC_NULL)
		MPI_Send(&data[0][cols-1], 1, new_type, neighbour_right, 0, MPI_COMM_WORLD);

	if(neighbour_left != MPI_PROC_NULL && neighbour_right != MPI_PROC_NULL)
	//MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,recvcount, recvtype, source, recvtag, comm, status, ierror)
		MPI_Sendrecv(data[cols-1], 1, new_type, neighbour_right, 0, data, 1, new_type, neighbour_left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if(neighbour_left != MPI_PROC_NULL && neighbour_right == MPI_PROC_NULL)
		MPI_Recv(data, 1, new_type, neighbour_left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// send data north -> south	
	if(neighbour_up == MPI_PROC_NULL && neighbour_down != MPI_PROC_NULL)
		MPI_Send(data[1], cols_halo, MPI_FLOAT, neighbour_down, 0, MPI_COMM_WORLD);

	if(neighbour_up != MPI_PROC_NULL && neighbour_down != MPI_PROC_NULL)
		MPI_Sendrecv(data[1], cols_halo, MPI_FLOAT, neighbour_down, 0, data, cols_halo, MPI_FLOAT, neighbour_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if(neighbour_up != MPI_PROC_NULL && neighbour_down == MPI_PROC_NULL)
		MPI_Recv(data, cols_halo, MPI_FLOAT, neighbour_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


	// send data east -> west	
	if(neighbour_left != MPI_PROC_NULL && neighbour_right == MPI_PROC_NULL)
		MPI_Send(&data[0][1], 1, new_type, neighbour_left, 0, MPI_COMM_WORLD);

	if(neighbour_left != MPI_PROC_NULL && neighbour_right != MPI_PROC_NULL)
		MPI_Sendrecv(&data[0][1], 1, new_type, neighbour_left, 0, &data[0][cols], 1, new_type, neighbour_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if(neighbour_left == MPI_PROC_NULL && neighbour_right != MPI_PROC_NULL)
		MPI_Recv(&data[0][cols], 1, new_type, neighbour_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
	MPI_Type_free(&new_type);

	// send data south -> north	
	if(neighbour_up != MPI_PROC_NULL && neighbour_down == MPI_PROC_NULL)
		MPI_Send(data[1], cols_halo, MPI_FLOAT, neighbour_up, 0, MPI_COMM_WORLD);

	if(neighbour_up != MPI_PROC_NULL && neighbour_down != MPI_PROC_NULL)
		MPI_Sendrecv(data[1], cols_halo, MPI_FLOAT, neighbour_up, 0, data[rows+1], cols_halo, MPI_FLOAT, neighbour_down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if(neighbour_up == MPI_PROC_NULL && neighbour_down != MPI_PROC_NULL)
		MPI_Recv(data[rows], cols_halo, MPI_FLOAT, neighbour_down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


	// Output eve< rank
	for(int rank=0; rank<num_procs; ++rank){
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == rank){
			cout << "From process: " << my_rank << " coords: " << "(" << coords[1] << "," << coords[0] << ")" <<endl;
			for(int i=0; i < rows_halo; ++i){
				for(int j=0; j < cols_halo; ++j)
					cout << data[i][j] << "\t";
				cout << "\n";
			}
			cout << "\n" << endl;
		}	
	}

	MPI_Finalize();
	return 0;
}