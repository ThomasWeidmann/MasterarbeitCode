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
#include "kamping/data_buffer.hpp"
#include "kamping/named_parameters.hpp"

#include "karam/mpi/grid_alltoall.hpp"

#include "timer.cpp"
#include "generator.cpp"
#include "test.cpp"
#include "interfaces.cpp"
#include "local_contraction.cpp"
#include "forest_local_contraction.cpp"
#include "forest_local_contraction2.cpp"

#include "helper_functions.cpp"
#include "analyze_instances.cpp"

#include "communicator.cpp"
//#include "karam/mpi/grid_alltoall.hpp"
//#include "list_ranking/example.cpp"
#include "list_ranking/asynchron_ruling_set2.cpp"
#include "list_ranking/grid_asynchron_ruling_set2.cpp"


#include "list_ranking/regular_ruling_set.cpp"
#include "list_ranking/regular_pointer_doubling.cpp"
#include "list_ranking/sequential_list_ranking.cpp"
#include "list_ranking/regular_ruling_set2.cpp"

#include "grid_list_ranking/grid_regular_ruling_set2.cpp"


#include "tree_rooting/forest_regular_ruling_set2.cpp"

#include "tree_rooting/tree_euler_tour.cpp"
#include "tree_rooting/forest_euler_tour.cpp"
#include "tree_rooting/forest_load_balance_regular_ruling_set2.cpp"


#include "tree_rooting/real_load_balance.cpp"
#include "tree_rooting/shuffle_load_balance.cpp"

int mpi_rank, mpi_size;






void error(std::string output)
{
	if (mpi_rank == 0)
		std::cout << output << std::endl;
}

int arg_index = 1;
char* get_next_arg(char* argv[])
{
	return argv[arg_index++];
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
	std::string ruling_set_rec = "ruling_set_rec";
	std::string ruling_set2 = "ruling_set2";
	std::string ruling_set2_rec = "ruling_set2_rec";
	std::string sequential = "sequential";
	std::string pointer_doubling = "pointer_doubling";
	std::string tree_rooting = "tree_rooting";
	std::string tree_rooting_rec = "tree_rooting_rec";
	std::string euler_tour = "euler_tour";
	std::string forest_rooting_euler = "forest_euler_tour";
	std::string lb_tree_rooting = "lb_tree_rooting";
	std::string asynchron = "asynchron";
	std::string asynchron2 = "asynchron2";
	std::string grid_asynchron = "grid_asynchron";
	std::string grid_ruling_set2 = "grid_ruling_set2";
	std::string grid_ruling_set2_rec = "grid_ruling_set2_rec";
	std::string local_contract = "local_contract";
	
	if (argc < 2)
	{
		error("First argument is the type of algorithm: " + ruling_set + ", " + pointer_doubling + ", " + sequential);
	}
	else
	{
		std::vector<std::uint64_t> send_buf(comm.size(),std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
		std::vector<std::int32_t> send_counts(comm.size(),1);
		auto recv = comm.alltoallv(kamping::send_buf(send_buf), kamping::send_counts(send_counts)).extract_recv_buffer();
		karam::mpi::GridCommunicator grid_comm;
		
		int communication_mode = 0; //aka direct
		
		
		//#######################################################################################################################################
		bool grid = atoi(get_next_arg(argv));
		
		char* algorithm_name = get_next_arg(argv);
		
		if (ruling_set.compare(algorithm_name) == 0)
		{
			std::int64_t num_local_vertices = atoi(get_next_arg(argv));
			std::int64_t dist_rulers = atoi(get_next_arg(argv));
			std::uint32_t num_iterations = atoi(get_next_arg(argv));
			
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			regular_ruling_set algorithm = regular_ruling_set(s, dist_rulers, num_iterations, grid);
			algorithm.start(s, comm, grid_comm);
			std::vector<std::int64_t> d = algorithm.get_ranks();
			test::regular_test_ranks(comm, s, d);

		}
		else if (ruling_set2.compare(algorithm_name) == 0)
		{
			
			std::int64_t num_local_vertices = atoi(get_next_arg(argv));
			//std::int64_t dist_rulers = atoi(get_next_arg(argv));
			//std::uint32_t num_iterations = atoi(get_next_arg(argv));
			
			
			
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			
			regular_ruling_set2 algorithm(s, grid, comm);
			std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
			
			
			
			
			test::regular_test_ranks(comm, s, d);

		}
		else if (grid_ruling_set2.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			int communication_mode = atoi(argv[4]);
			
			grid_regular_ruling_set2 algorithm(s, dist_rulers, 1, communication_mode);
			std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
			
			test::regular_test_ranks(comm, s, d);

		}
		else if (grid_ruling_set2_rec.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			int communication_mode = atoi(argv[4]);
			
			grid_regular_ruling_set2 algorithm(s, dist_rulers, 2, communication_mode);
			std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
			
			test::regular_test_ranks(comm, s, d);

		}
		else if (ruling_set2_rec.compare(argv[1]) == 0)
		{
			/*
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			regular_ruling_set2 algorithm(s, dist_rulers, 2);
			std::vector<std::int64_t> d = algorithm.start(comm);
			
			test::regular_test_ranks(comm, s, d);*/
		}
		else if (pointer_doubling.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);

			std::vector<std::uint64_t> s = generator::generate_regular_wood_vector(num_local_vertices, comm);
			regular_pointer_doubling algorithm(s, comm, communication_mode);
			
			std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
			
			test::regular_test_ranks(comm, s, d);
			analyze_instances::analyze_regular_instance(s, comm);
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
			std::vector<std::uint64_t> tree_vector = generator::generate_regular_wood_vector_from_irregular_graph(num_local_vertices, comm);
			std::int32_t dist_rulers = atoi(argv[3]);
			forest_regular_ruling_set2 algorithm(tree_vector, dist_rulers, comm,1);
			std::vector<std::int64_t> d = algorithm.result_dist;
			std::vector<std::uint64_t> roots = algorithm.result_root;
			test::regular_test_ranks_and_roots(comm, tree_vector, d, roots);
			//std::vector<std::int64_t> d = forest_regular_ruling_set2(tree_vector, dist_rulers, comm,1).result_dist;
			
			//analyze_instances::analyze_regular_instance(tree_vector, comm);
			//test::regular_test_ranks(comm, tree_vector, d);
		}
		else if (tree_rooting_rec.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::vector<std::uint64_t> tree_vector = generator::generate_regular_wood_vector(num_local_vertices, comm);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::int64_t> d = forest_regular_ruling_set2(tree_vector, dist_rulers, comm,2).result_dist;
			
			//analyze_instances::analyze_regular_instance(tree_vector, comm);
			test::regular_test_ranks(comm, tree_vector, d);
		}
		else if (euler_tour.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> tree_vector = generator::generate_regular_tree_vector(num_local_vertices, comm);

			std::vector<std::int64_t> d = tree_euler_tour(comm, tree_vector, dist_rulers).start(comm, tree_vector);
			test::regular_test_ranks(comm, tree_vector, d);
		}
		else if (lb_tree_rooting.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_tree_vector(num_local_vertices, comm);

			//forest_load_balance_regular_ruling_set2(dist_rulers).start2(s,comm);
			//real_load_balance(dist_rulers).start(s,comm, grid_comm);
			std::vector<std::int64_t> d = shuffle_load_balance().start(s,comm, grid_comm, dist_rulers);
			test::regular_test_ranks(comm, s, d);

			//analyze_instances::analyze_regular_instance(s, comm);
		}
		else if (local_contract.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_wood_vector(num_local_vertices, comm);
			std::cout << mpi_rank << " with graph korrekt insanziiert" << std::endl;	
			forest_local_contraction2 algorithm;
			std::vector<std::int64_t> d = algorithm.start(comm, s);
			
			test::regular_test_ranks(comm, s, d);
			analyze_instances::analyze_regular_instance(s, comm);
			/*
			kagen::KaGen gen(MPI_COMM_WORLD);
			//rgg, rmat(a=0.59,b=0.19,c=0.19, d=1-(a+b+c)) aus graph500 challenge
			kagen::Graph graph = gen.GenerateGrid2D_N(4, 1, false);
			for (auto const& [src, dst]: graph.edges)
				std::cout << "(" << src << "," << dst << ")" << std::endl;*/
		
			
		}
		else if (forest_rooting_euler.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_wood_vector(num_local_vertices, comm);

			std::vector<std::int64_t> d = forest_euler_tour(comm, s, dist_rulers).start(comm, s);
			test::regular_test_ranks(comm, s, d);
			
			//analyze_instances::analyze_regular_instance(s, comm);

			
		}
		else if (asynchron.compare(argv[1]) == 0)
		{
			/*
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			example example;
			
			std::vector<std::int64_t> d = example.test(s, dist_rulers, comm);
			
			test::regular_test_ranks(comm, s, d);*/
		}
		else if (asynchron2.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			asynchron_ruling_set2 algorithm;
			
			std::vector<std::int64_t> d = algorithm.test(s, dist_rulers, comm);
			
			test::regular_test_ranks(comm, s, d);
		}
		else if (grid_asynchron.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			grid_asynchron_ruling_set2 algorithm;
			
			std::vector<std::int64_t> d = algorithm.test(s, dist_rulers, comm);
			
			test::regular_test_ranks(comm, s, d);
		}
		else
		{		
			
			grid = true;
			
			std::int64_t num_local_vertices = 1000000;
			
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			int num_runs = 1;
			
			std::uint32_t num_iterations = 1;
			
			for (int num_iterations = 1; num_iterations < 5;  num_iterations++)
			{
				for (int dist_rulers = 50; dist_rulers <= 200; dist_rulers += 10)
				{
					for (int runs = 0; runs < num_runs; runs++)
					{
						regular_ruling_set2 algorithm(s, dist_rulers, num_iterations, grid);
						std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
					}
				}
				
			}
			
			//regular_ruling_set2 algorithm(s, 50, 3, grid);
			//std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
			
		
			
			
			
			
			
	
	/*
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_wood_vector_from_irregular_graph(num_local_vertices, comm);
			forest_local_contraction2 algorithm;
			std::vector<std::int64_t> d = algorithm.start(comm, s);
			
			test::regular_test_ranks(comm, s, d);
			//analyze_instances::analyze_regular_instance(s, comm);*/
			//error(std::string(argv[1]) + " is not a name of an algorithm or wrong parameters");
			
			
		}
	}
	
	
	MPI_Finalize();
	return 0;
	
};

