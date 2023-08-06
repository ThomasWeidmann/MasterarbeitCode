#include <math.h> 
#include <kagen.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

#include <iostream>
#include <numeric>
#include <vector>

#include <mpi.h>

#include "kamping/examples/usage/helpers_for_examples.hpp"
#include "kamping/checking_casts.hpp"
#include "kamping/collectives/alltoall.hpp"
#include "kamping/communicator.hpp"
#include "kamping/environment.hpp"


int mpi_rank, mpi_size;

struct unidirectional_path {
  std::vector<std::int32_t> s;
  std::int32_t num_global_vertices;
  std::int32_t num_local_vertices;
  std::vector<std::int32_t> num_vertices_per_pe;
  std::vector<std::int32_t> prefix_sum_num_vertices_per_pe;
};

unidirectional_path generate_unidirectional_path(std::int32_t num_global_vertices);
void sparse_ruling_set(unidirectional_path unidirectional_path, std::int32_t dist_rulers, kamping::Communicator<>& comm);



int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	srand((unsigned) time(NULL) + mpi_rank);
	if (argc != 3){
		if (mpi_rank == 0)
			std::cout << "Error\nParameter 1: Number of global vertices\nParameter 2: Distance of rulers" << std::endl;
		MPI_Finalize();
		return 0;
	}	


	kamping::Environment e;
	kamping::Communicator<>        comm;

	std::int32_t num_global_vertices = atoi(argv[1]);
	std::int32_t dist_rulers = atoi(argv[2]);
	
	unidirectional_path unidirectional_path = generate_unidirectional_path(num_global_vertices);
	sparse_ruling_set(unidirectional_path, dist_rulers, comm);


/*
    // Rank i sends i values to rank 0, i+1 values to rank 1, ...
    std::vector<int> counts_per_rank(comm.size());
    std::iota(counts_per_rank.begin(), counts_per_rank.end(), comm.rank_signed());

    int                 num_elements = std::reduce(counts_per_rank.begin(), counts_per_rank.end(), 0);
    std::vector<size_t> input(asserting_cast<size_t>(num_elements));
    // Rank i sends it own rank to all others
    std::fill(input.begin(), input.end(), comm.rank());

	auto recv = comm.alltoallv(send_buf(input), send_counts(counts_per_rank));

    std::vector<size_t> output = recv.extract_recv_buffer();

	std::vector<int> recvdislp  = recv.extract_recv_counts();

    print_result_on_root(output, comm);

	print_result_on_root(recvdislp, comm);*/
	MPI_Finalize();
    return 0;
}

std::int32_t calculate_targetPE(std::int32_t global_index, std::int32_t num_global_vertices)
{
	std::int32_t lower_bound_num_local_vertices = num_global_vertices / mpi_size;
	std::int32_t upper_bound_num_local_vertices = lower_bound_num_local_vertices + 1;
	std::int32_t rest = num_global_vertices - mpi_size * lower_bound_num_local_vertices;
	
	if (global_index < rest * upper_bound_num_local_vertices)
		return global_index / upper_bound_num_local_vertices;

	return rest + (global_index - (rest * upper_bound_num_local_vertices)) / lower_bound_num_local_vertices;
}

void sparse_ruling_set(unidirectional_path unidirectional_path, std::int32_t dist_rulers, kamping::Communicator<>& comm)	
{
	
	
	using namespace kamping;
	
	//an element is a ruler, iff its global_index % dist_rulers == 0, this way we assure every PE has same number of rulers
	std::vector<std::int32_t> s = unidirectional_path.s;
	std::int32_t num_global_vertices = unidirectional_path.num_global_vertices;
	std::int32_t num_local_vertices = unidirectional_path.num_local_vertices;
	std::vector<std::int32_t> num_vertices_per_pe = unidirectional_path.num_vertices_per_pe;
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe = unidirectional_path.prefix_sum_num_vertices_per_pe;
	
	std::vector<std::int32_t> counts_per_rank(mpi_size,0);
	std::int32_t first_local_ruler_index = ((prefix_sum_num_vertices_per_pe[mpi_rank] + dist_rulers - 1) / dist_rulers) * dist_rulers - prefix_sum_num_vertices_per_pe[mpi_rank]; //every PE must have a ruler, otherwise undefined behavior (probably array out of bounds later)


	std::vector<std::int32_t> is_reached(num_local_vertices, 0);

	for (std::int32_t i = first_local_ruler_index; i < num_local_vertices; i+= dist_rulers)
	{
		counts_per_rank[calculate_targetPE(s[i], num_global_vertices)]++;
		is_reached[i] = 1;
	}
	std::vector<std::int32_t> send_displacements(mpi_size + 1,0);
	for (std::int32_t i = 1; i < mpi_size + 1; i++)
		send_displacements[i] = send_displacements[i-1] + counts_per_rank[i-1];
	
	
	
	std::vector<std::int32_t> out_buffer(2 * send_displacements[mpi_size]);
	
	//packets in first round are defined like in practical_parallel_list_ranking aka (suc(r),r) --> suc(r) of ith packet = out_buffer[2*i], r of ith packet = out_buffer[2*i+1]
	std::fill(counts_per_rank.begin(), counts_per_rank.end(), 0);
	
	for (std::int32_t i = first_local_ruler_index; i < num_local_vertices; i+= dist_rulers)
	{
		
		std::int32_t targetPE = calculate_targetPE(s[i], num_global_vertices);
		
		std::int32_t packet_index = send_displacements[targetPE]+counts_per_rank[targetPE]++;
		out_buffer[2*packet_index] = s[i];
		out_buffer[2*packet_index+1] = i + prefix_sum_num_vertices_per_pe[mpi_rank];
	}
	
	for (std::int32_t i = 0; i < mpi_size; i++)
		counts_per_rank[i] *= 2;
	
	auto recv = comm.alltoallv(send_buf(out_buffer), send_counts(counts_per_rank));

	std::vector<std::int32_t> recv_buffer = recv.extract_recv_buffer();

	std::vector<std::int32_t> recv_displs  = recv.extract_recv_displs();
	
	std::vector<std::int32_t> recv_counts  = recv.extract_recv_counts();
	
	//print_result(s, comm);
	//print_result(recv_buffer, comm);
	//print_result(recv_displs, comm);
	
	//now packets are forwarded in a loop, iteration < dist_rulers ohne testen
	for (std::int32_t iteration = 0; iteration < 2 * dist_rulers; iteration++)
	{
		std::fill(counts_per_rank.begin(), counts_per_rank.end(), 0);
		
		out_buffer.resize(recv_buffer.size()); //er wird maximal so groß, da pakete die auf ruler zeigen, "verschlungen" werden
		
		for (std::int32_t i = 0; i < recv_buffer.size() / 2; i++)
		{
			if (recv_buffer[2*i] % dist_rulers != 0) //iff paket doesn't point auf ruler
			{
				std::int32_t succ_node_global_index = s[recv_buffer[2*i] - prefix_sum_num_vertices_per_pe[mpi_rank]];
				counts_per_rank[calculate_targetPE(succ_node_global_index, num_global_vertices)]++;
				
				is_reached[recv_buffer[2*i]- prefix_sum_num_vertices_per_pe[mpi_rank]] = 1;
			}	
		}		
		for (std::int32_t i = 1; i < mpi_size + 1; i++)
				send_displacements[i] = send_displacements[i-1] + counts_per_rank[i-1];
		std::fill(counts_per_rank.begin(), counts_per_rank.end(), 0);

		//jetzt wurden displacements für out_buffer berechnet und jetzt wird out_buffer befüllt
		
		for (std::int32_t i = 0; i < recv_buffer.size() / 2; i++)
		{
			if (recv_buffer[2*i] % dist_rulers != 0) //iff paket doesn't point auf ruler
			{
				std::int32_t succ_node_global_index = s[recv_buffer[2*i] - prefix_sum_num_vertices_per_pe[mpi_rank]];
				
				if (recv_buffer[2*i] == succ_node_global_index)//iff paket points auf list end
					continue; 
				
				std::int32_t targetPE = calculate_targetPE(succ_node_global_index, num_global_vertices);
		
				std::int32_t packet_index = send_displacements[targetPE]+counts_per_rank[targetPE]++;
				
				out_buffer[2*packet_index] = succ_node_global_index;
				out_buffer[2*packet_index+1] = recv_buffer[2*i+1];
				
			}	
		}
		
		for (std::int32_t i = 0; i < mpi_size; i++)
			counts_per_rank[i] *= 2;
		
		auto recv = comm.alltoallv(send_buf(out_buffer), send_counts(counts_per_rank));

		recv_buffer = recv.extract_recv_buffer();

		recv_displs  = recv.extract_recv_displs();
		
		recv_counts  = recv.extract_recv_counts();
		
		
		std::int32_t num_not_reached_nodes = 0;
		for (std::int32_t i = 0; i < num_local_vertices; i++)
			if (is_reached[i] == 0)
				num_not_reached_nodes++;
			
		if (mpi_rank == 0)
			std::cout << "iteration " << iteration << " num_not_reached_nodes=" << num_not_reached_nodes << std::endl;
		
	}
	
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