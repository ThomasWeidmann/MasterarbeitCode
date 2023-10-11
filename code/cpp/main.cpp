#include <math.h> 
#include <kagen.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

#include <iostream>
#include <numeric>
#include <vector>
#include <iostream>
#include <algorithm>
#include <mpi.h>


#include "kamping/checking_casts.hpp"
#include "kamping/collectives/alltoall.hpp"
#include "kamping/communicator.hpp"
#include "kamping/environment.hpp"
#include "kamping/collectives/allgather.hpp"
#include "kamping/collectives/bcast.hpp"

#include "timer.cpp"
#include "generator.cpp"
#include "test.cpp"

#include "list_ranking/regular_ruling_set.cpp"
#include "list_ranking/regular_pointer_doubling.cpp"
#include "list_ranking/sequential_list_ranking.cpp"
#include "list_ranking/regular_ruling_set2.cpp"



#include "tree_rooting/wood_regular_ruling_set2.cpp"

#include "tree_rooting/tree_euler_tour.cpp"

int mpi_rank, mpi_size;




std::vector<std::uint64_t> generate_regular_successor_vector(std::uint64_t num_local_vertices);
std::vector<std::uint64_t> generate_regular_wood_vector(std::uint64_t num_local_vertices, kamping::Communicator<>& comm);
std::vector<std::uint64_t> generate_regular_tree_vector(std::uint64_t num_local_vertices, kamping::Communicator<>& comm);


void error(std::string output)
{
	if (mpi_rank == 0)
		std::cout << output << std::endl;
}

//das erste argument ist der algorithmus aka "ruling_set", "sequential", "pointer_doubling"
int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	srand((unsigned) time(NULL) + mpi_rank);
	
	
	
	kamping::Environment e;
	kamping::Communicator<> comm;
	std::string ruling_set = "ruling_set";
	std::string ruling_set2 = "ruling_set2";
	std::string sequential = "sequential";
	std::string pointer_doubling = "pointer_doubling";
	std::string tree_rooting = "tree_rooting";
	std::string euler_tour = "euler_tour";
	std::string init_mpi = "init_mpi";
	
	if (argc < 2)
	{
		error("First argument is the type of algorithm: " + ruling_set + ", " + pointer_doubling + ", " + sequential);
	}
	else
	{
		if (ruling_set.compare(argv[1]) == 0)
		{
			
			std::int64_t num_local_vertices = atoi(argv[2]);

			std::int64_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);

			regular_ruling_set algorithm(s, dist_rulers, 2);
			std::vector<std::int64_t> d = algorithm.start(comm);
			
			test::regular_test(comm, s, d);
		}
		else if (ruling_set2.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			regular_ruling_set2 algorithm(s, dist_rulers);
			std::vector<std::int64_t> d = algorithm.start(comm);
			
			test::regular_test(comm, s, d);

		}
		else if (pointer_doubling.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);

			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			regular_pointer_doubling algorithm(s, comm);
			
			algorithm.start(comm);
		}
		else if (sequential.compare(argv[1]) == 0 && mpi_size == 1)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);

			sequential_list_ranking algorithm(s);
			algorithm.start(comm);
		}
		else if (tree_rooting.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::vector<std::uint64_t> tree_vector = generator::generate_regular_tree_vector(num_local_vertices, comm);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::int64_t> d = wood_regular_ruling_set2(tree_vector, dist_rulers, comm).result_dist;
			
			test::regular_test(comm, tree_vector, d);
		}
		else if (euler_tour.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::vector<std::uint64_t> tree_vector = generator::generate_regular_tree_vector(num_local_vertices, comm);

			std::vector<std::int64_t> d = tree_euler_tour(comm, tree_vector).start(comm, tree_vector);
			test::regular_test(comm, tree_vector, d);
		}
		else if (init_mpi.compare(argv[1]) == 0)
		{
			std::vector<std::uint64_t> send_buf(comm.size(),std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
			std::vector<std::int32_t> send_counts(comm.size(),1);
			auto recv = comm.alltoallv(kamping::send_buf(send_buf), kamping::send_counts(send_counts)).extract_recv_buffer();
		}	
		else 
		{
			error(std::string(argv[1]) + " is not a name of an algorithm or wrong parameters");
		}
	}
	
	
	MPI_Finalize();
	return 0;
	
};

