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

#include "regular_ruling_set.cpp"
#include "timer.cpp"
#include "sequential_list_ranking.cpp"
#include "regular_ruling_set2.cpp"

int mpi_rank, mpi_size;

struct unidirectional_path {
  std::vector<std::int32_t> s;
  std::int32_t num_global_vertices;
  std::int32_t num_local_vertices;
  std::vector<std::int32_t> num_vertices_per_pe;
  std::vector<std::int32_t> prefix_sum_num_vertices_per_pe;
};

unidirectional_path generate_unidirectional_path(std::int32_t num_global_vertices);
std::vector<std::uint64_t> generate_regular_successor_vector(std::uint64_t num_local_vertices);

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
	
	if (argc < 2)
	{
		error("First argument is the type of algorithm: " + ruling_set + ", " + pointer_doubling + ", " + sequential);
	}
	else
	{
		if (ruling_set.compare(argv[1]) == 0)
		{
			std::int32_t num_global_vertices = atoi(argv[2]) * mpi_size;

			std::int32_t dist_rulers = atoi(argv[3]);
			unidirectional_path unidirectional_path = generate_unidirectional_path(num_global_vertices);

			regular_ruling_set algorithm(unidirectional_path.s, dist_rulers, 2);
			algorithm.start(comm);
		}
		else if (ruling_set2.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generate_regular_successor_vector(num_local_vertices);
			regular_ruling_set2 algorithm(s, dist_rulers);
			algorithm.start(comm);
		}
		else if (pointer_doubling.compare(argv[1]) == 0)
		{
			std::int32_t num_global_vertices = atoi(argv[2]);
			num_global_vertices  = mpi_size * (num_global_vertices / mpi_size);
			unidirectional_path unidirectional_path = generate_unidirectional_path(num_global_vertices);

			regular_pointer_doubling algorithm(unidirectional_path.s, comm);
			algorithm.start(comm);
		}
		else if (sequential.compare(argv[1]) == 0 && mpi_size == 1)
		{
			std::int32_t num_global_vertices = atoi(argv[2]);
			unidirectional_path unidirectional_path = generate_unidirectional_path(num_global_vertices);

			sequential_list_ranking algorithm(unidirectional_path.s);
			algorithm.start(comm);
		}
		else 
		{
			error(std::string(argv[1]) + " is not a name of an algorithm or wrong parameters");
		}
	}
	
	MPI_Finalize();
	return 0;
	
}

//every PE has same number of vertices
std::vector<std::uint64_t> generate_regular_successor_vector(std::uint64_t num_local_vertices)
{
	
	kagen::KaGen gen(MPI_COMM_WORLD);
	std::vector<std::uint64_t> s(num_local_vertices);
	std::uint64_t num_global_vertices = mpi_size * num_local_vertices;
	auto path = gen.GenerateDirectedPath(num_global_vertices, true);
	
	for (std::uint64_t i = 0; i < num_local_vertices; i++)
		s[i] = i + mpi_rank * num_local_vertices;
	
	for (auto const& [src, dst]: path.edges)
		s[src - mpi_rank * num_local_vertices] = dst;
	

	return s;
}

unidirectional_path generate_unidirectional_path(std::int32_t num_global_vertices)
{
	kagen::KaGen gen(MPI_COMM_WORLD);
	
	auto path = gen.GenerateDirectedPath(num_global_vertices, true);
    std::int32_t num_local_vertices = path.vertex_range.second - path.vertex_range.first;
	

	std::vector<std::int32_t> num_vertices_per_pe(mpi_size);
	MPI_Allgather(&num_local_vertices, 1, MPI_INT, &*num_vertices_per_pe.begin(), 1, MPI_INT, MPI_COMM_WORLD); 
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe(mpi_size + 1);
	prefix_sum_num_vertices_per_pe[0]=0;
	for (std::int32_t i = 0; i < mpi_size; i++)
		prefix_sum_num_vertices_per_pe[i+1] = prefix_sum_num_vertices_per_pe[i] + num_vertices_per_pe[i];
	
	std::vector<std::int32_t> s(num_local_vertices);
	for (std::int32_t i = 0; i < num_local_vertices; i++)
		s[i] = i + prefix_sum_num_vertices_per_pe[mpi_rank]; //last edge has pointer to itself and not no pointer like kagen does

	for (auto const& [src, dst]: path.edges){
		s[src - prefix_sum_num_vertices_per_pe[mpi_rank]] = dst; 
	}

	unidirectional_path unidirectional_path;
	unidirectional_path.s = s;
	unidirectional_path.num_global_vertices = num_global_vertices;
	unidirectional_path.num_local_vertices = num_local_vertices;
	unidirectional_path.num_vertices_per_pe = num_vertices_per_pe;
	unidirectional_path.prefix_sum_num_vertices_per_pe = prefix_sum_num_vertices_per_pe;
	
	return unidirectional_path;
}