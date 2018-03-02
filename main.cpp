#include <mpi.h>
#include <iostream>
#include <array>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

int main(int argc, char** argv){
	
	MPI_Init(&argc, &argv);

	// MPI
	int my_rank, num_procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// parameter, matrix size m×m
	int matrix_size = stoi(argv[1]);
	int iterations = stoi(argv[2]);
	
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

	
	// prepare space for HALO		
	int offset_left=0;
	int offset_right=0;
	int offset_up=0;
	int offset_down=0;
	int cols_halo = cols;
	int rows_halo = rows;
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
	// vector<vector<float, cols_halo>, rows_halo> data;
	// initialize container with 0
	for(int i=0; i<=rows_halo; ++i)
		for(int j=0; j<=cols_halo;++j)
			data[i][j] = 0;

	// fill with data
	for(int i=0+offset_up; i<rows+offset_up; ++i)
		for(int j=0+offset_left; j<cols+offset_left;++j)
			data[i][j] = 0.2;

	// heat soruce
	// border
	// upper border
	if(neighbour_up == MPI_PROC_NULL)
		for(int i=offset_left; i < cols + offset_left; ++i)
			data[0][i] = 1;
	// lower border 
	if(neighbour_down == MPI_PROC_NULL)
		for(int i=offset_left; i < cols + offset_left; ++i)
			data[rows-1+offset_up][i] = 1;
	// left border
	if(neighbour_left == MPI_PROC_NULL)
		for(int i=offset_up; i < rows + offset_up; ++i)
			data[i][0] = 1;
	// right border
	if(neighbour_right == MPI_PROC_NULL)
		for(int i=offset_up; i < rows + offset_up; ++i)
			data[i][cols-1+offset_left] = 1;

	if(my_rank == 0)
		data[cols/2][rows/2] = 20;


	// 	halo communication    //
	// send data to neighbour //

	for(int i=0; i <iterations; ++i){

		// Output eve< rank
		for(int rank=0; rank<num_procs; ++rank){
			MPI_Barrier(MPI_COMM_WORLD);
			if(my_rank == rank){
				cout << "From process: " << my_rank << " coords: " << "(" << coords[1] << "," << coords[0] << ")" <<endl;
				for(int i=0; i < rows_halo; ++i){
					for(int j=0; j < cols_halo; ++j)
						cout << setprecision(2)<< data[i][j] << "\t";
					cout << "\n";
				}
				cout << "\n" << endl;
			}	
		}

		// send data west -> east	
		MPI_Datatype new_type;
			// int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
		MPI_Type_vector(rows, 1, cols_halo, MPI_FLOAT, &new_type);
		MPI_Type_commit(&new_type);

		if(neighbour_left == MPI_PROC_NULL && neighbour_right != MPI_PROC_NULL)
			MPI_Send(&data[offset_up][cols-1], 1, new_type, neighbour_right, 0, MPI_COMM_WORLD);

		if(neighbour_left != MPI_PROC_NULL && neighbour_right != MPI_PROC_NULL)
		//MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,recvcount, recvtype, source, recvtag, comm, status, ierror)
			MPI_Sendrecv(&data[offset_up][cols-1], 1, new_type, neighbour_right, 0, &data[offset_up][0], 1, new_type, neighbour_left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(neighbour_left != MPI_PROC_NULL && neighbour_right == MPI_PROC_NULL)
			MPI_Recv(&data[offset_up][0], 1, new_type, neighbour_left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// // send data north -> south	
		if(neighbour_up == MPI_PROC_NULL && neighbour_down != MPI_PROC_NULL)
			MPI_Send(&data[rows-1][offset_left], cols, MPI_FLOAT, neighbour_down, 0, MPI_COMM_WORLD);

		if(neighbour_up != MPI_PROC_NULL && neighbour_down != MPI_PROC_NULL)
			MPI_Sendrecv(&data[rows-1+offset_up][offset_left], cols, MPI_FLOAT, neighbour_down, 0, &data[0][offset_left], cols, MPI_FLOAT, neighbour_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(neighbour_up != MPI_PROC_NULL && neighbour_down == MPI_PROC_NULL)
			MPI_Recv(&data[0][offset_left], cols, MPI_FLOAT, neighbour_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


		// send data east -> west	
		if(neighbour_left != MPI_PROC_NULL && neighbour_right == MPI_PROC_NULL)
			MPI_Send(&data[offset_up][offset_left], 1, new_type, neighbour_left, 0, MPI_COMM_WORLD);

		if(neighbour_left != MPI_PROC_NULL && neighbour_right != MPI_PROC_NULL)
			MPI_Sendrecv(&data[offset_up][offset_left], 1, new_type, neighbour_left, 0, &data[offset_up][cols], 1, new_type, neighbour_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(neighbour_left == MPI_PROC_NULL && neighbour_right != MPI_PROC_NULL)
			MPI_Recv(&data[offset_up][cols], 1, new_type, neighbour_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
		MPI_Type_free(&new_type);

		// send data south -> north	
		if(neighbour_up != MPI_PROC_NULL && neighbour_down == MPI_PROC_NULL)
			MPI_Send(&data[offset_up][offset_left], cols, MPI_FLOAT, neighbour_up, 0, MPI_COMM_WORLD);

		if(neighbour_up != MPI_PROC_NULL && neighbour_down != MPI_PROC_NULL)
			MPI_Sendrecv(&data[offset_up][offset_left], cols, MPI_FLOAT, neighbour_up, 0, &data[rows+offset_up][offset_left], cols, MPI_FLOAT, neighbour_down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(neighbour_up == MPI_PROC_NULL && neighbour_down != MPI_PROC_NULL)
			MPI_Recv(&data[rows][offset_left], cols, MPI_FLOAT, neighbour_down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Barrier(MPI_COMM_WORLD);

		////////////////////////////
		// 		calculation		  //
		////////////////////////////
		float data_old[rows_halo][cols_halo]; 
		for(int m=0; m < rows_halo;++m)
			for (int n=0; n < cols_halo; ++n)
				data_old[m][n] = data[m][n];

		for(int i = 0 + offset_up; i < rows+offset_up; ++i){
			for(int j = 0 + offset_left; j < cols + offset_left; ++j){
				
				int neighbours = 0;
				float north = 0;
				float south = 0;
				float east = 0;
				float west = 0;

				float core_factor = 1;
				float n_factor = 1;
				// north
				if(i-1 >= 0){
					north = n_factor * data_old[i-1][j];
					++neighbours;
				}

				// south
				if(i+1 < rows_halo){
					south = n_factor * data_old[i+1][j];
					++neighbours;				
				}

				// east
				if(j-1 >= 0){
					east = n_factor * data_old[i][j-1];
					++neighbours;
				}

				// west 
				if(j+1<cols_halo){
					west = n_factor * data_old[i][j+1];
					++neighbours;
				}
				float core = data[i][j];
				float neighbours_avg = 0.15*( north + south + east + west);

				data[i][j] = 0.4*core + neighbours_avg;
			}
		}
	} // Iteration end


	// if(my_rank == 0)
		// float matrix[matrix_size][matrix_size];

	// int displs[num_procs];
	// displs[my_rank] = cols*rows;

	// // send all data back to process 0 for output
	// /*MPI_Gatherv(
	// 	const void *sendbuf, 
	// 	int sendcount, 
	// 	MPI_Datatype sendtype,
	// 	void *recvbuf, 
	// 	const int recvcounts[],
	// 	const int displs[], 
	// 	MPI_Datatype recvtype, 
	// 	int root, MPI_Comm comm)*/
	// MPI_Gatherv(data[offset_up][offset_left], cols*rows, MPI_FLOAT, matrix, matrix_size*matrix_size, const int displs[], MPI_FLOAT, 0, MPI_COMM_WORLD);
	

	// if(my_rank == 0){
	// 	// Output eve< rank
	// 	cout << "\n\nFinal Output" <<endl;
	// 	for(int i=0; i < matrix_size; ++i){
	// 		for(int j=0; j < matrix_size; ++j)
	// 			cout << setprecision(2)<< matrix[i][j] << "\t";
	// 		cout << "\n" << endl;
	// 	}

	// }





	MPI_Finalize();
	return 0;
}