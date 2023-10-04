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

#include "tree_rooting/tree_regular_ruling_set2.cpp"

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
		else if (tree_rooting.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::vector<std::uint64_t> tree_vector = generate_regular_tree_vector(num_local_vertices, comm);
			std::int32_t dist_rulers = atoi(argv[3]);
			tree_regular_ruling_set2(tree_vector, dist_rulers, comm);
		}
		else 
		{
			error(std::string(argv[1]) + " is not a name of an algorithm or wrong parameters");
		}
	}
	
	MPI_Finalize();
	return 0;
	
}

std::uint64_t hash64(std::uint64_t x) {
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
}

std::vector<std::uint64_t> generate_regular_tree_vector(std::uint64_t num_local_vertices, kamping::Communicator<>& comm)
{
	std::vector<std::uint64_t> s(num_local_vertices);
	std::uint64_t node_offset = num_local_vertices * comm.rank();
	for (std::uint64_t i = 0; i < num_local_vertices; i++)
	{
		if (i == 0 && comm.rank() == 0)
			s[i] = 0;
		else
			s[i] = hash64(i + node_offset) % (i + node_offset);
	}
	return s;
	
	
}

std::vector<std::uint64_t> generate_regular_wood_vector(std::uint64_t num_local_vertices, kamping::Communicator<>& comm)
{
	std::int32_t max_edge_weight = 12345678; 
	std::uint64_t node_offset = mpi_rank * num_local_vertices;
	kagen::KaGen gen(MPI_COMM_WORLD);
	auto graph = gen.GenerateUndirectedGNP(num_local_vertices * mpi_size, 0.5, false);
	//now every edge finds edge with lowest edge weight
	std::vector<std::uint64_t> s(num_local_vertices); //s[i] will be the node j that minimizes c(i,j) with
	std::iota(s.begin(), s.end(), node_offset);
	std::vector<std::int32_t> w(num_local_vertices, max_edge_weight + 1);
	

	for (auto const& [src, dst]: graph.edges)
	{
		std::string edge_string = src < dst ? std::to_string(src) + "," + std::to_string(dst): std::to_string(dst) + "," + std::to_string(src);
		std::uint64_t edge_weight = std::hash<std::string>{}(edge_string)  % (max_edge_weight + 1);
		
		edge_weight = hash64(hash64(src) + hash64(dst)) % (max_edge_weight +1);
		
		if (edge_weight < w[src - node_offset])
		{
			s[src - node_offset] = dst;
			w[src - node_offset] = edge_weight;
		}
		
	}
	
	//for each node u on this PE we calculate the lightest edge (u,v) and send this information to v
	struct lightest_edge {
		std::uint64_t source;
		std::uint64_t destination;
	};
	
	std::vector<std::int32_t> num_packets_per_PE(mpi_size,0);
	std::vector<std::int32_t> send_displacements(mpi_size + 1,0);
	for (std::uint64_t i = 0; i < num_local_vertices; i++)
	{
		std::int32_t targetPE = s[i] / num_local_vertices;
		num_packets_per_PE[targetPE]++;
	}
	std::vector<lightest_edge> send(num_local_vertices);
	for (std::int32_t i = 1; i < mpi_size + 1; i++)
		send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
	std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
	for (std::uint64_t i = 0; i < num_local_vertices; i++)
	{
		std::int32_t targetPE = s[i] / num_local_vertices;
		std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
		send[packet_index].source = i + node_offset;
		send[packet_index].destination = s[i];
	}
	
	auto recv = comm.alltoallv(kamping::send_buf(send), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
	for (std::uint64_t i = 0; i < recv.size(); i++)
	{
		lightest_edge e = recv[i];
		if (s[e.destination - node_offset] == e.source && e.destination < e.source) // der node mit kleinerer id wird root
			s[e.destination - node_offset] = e.destination;
	}

	return s;
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