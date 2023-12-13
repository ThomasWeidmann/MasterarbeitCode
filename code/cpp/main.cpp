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

#include "timer.cpp"
#include "generator.cpp"
#include "test.cpp"
#include "interfaces.cpp"
#include "local_contraction.cpp"
#include "forest_local_contraction.cpp"
#include "forest_local_contraction2.cpp"


#include "grid_all_to_all.cpp"
#include "helper_functions.cpp"
#include "analyze_instances.cpp"

#include "communicator.cpp"
//#include "karam/mpi/grid_alltoall.hpp"
#include "list_ranking/example.cpp"
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
		
		
		if (ruling_set.compare(argv[1]) == 0)
		{
			
			std::int64_t num_local_vertices = atoi(argv[2]);

			std::int64_t dist_rulers = atoi(argv[3]);
			/*
			double n = num_local_vertices * mpi_size;
			double dist = dist_rulers;
			double exact = std::log(0.5 / n) / std::log((dist - 1)/dist);
			double approx = dist * std::log(2*n);
			if (mpi_rank == 0) std::cout << "exact = " << exact << std::endl;
			if (mpi_rank == 0) std::cout << "approx = " << approx << std::endl;*/
			
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			regular_ruling_set algorithm = regular_ruling_set(s, dist_rulers, 1);
			std::vector<std::int64_t> d = algorithm.start(comm, s, grid_comm);
			test::regular_test(comm, s, d);
			
			

			
		}
		else if (ruling_set_rec.compare(argv[1]) == 0)
		{
			
			std::int64_t num_local_vertices = atoi(argv[2]);

			std::int64_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);

			regular_ruling_set algorithm = regular_ruling_set(s, dist_rulers, 2);
			
			std::vector<std::int64_t> d = algorithm.start(comm, s, grid_comm);
			
			test::regular_test(comm, s, d);
		}
		else if (ruling_set2.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			
			regular_ruling_set2 algorithm(s, dist_rulers, 1);
			std::vector<std::int64_t> d = algorithm.start(comm);
			
			test::regular_test(comm, s, d);

		}
		else if (grid_ruling_set2.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			int communication_mode = atoi(argv[4]);
			
			grid_regular_ruling_set2 algorithm(s, dist_rulers, 1, communication_mode);
			std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
			
			test::regular_test(comm, s, d);

		}
		else if (grid_ruling_set2_rec.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			int communication_mode = atoi(argv[4]);
			
			grid_regular_ruling_set2 algorithm(s, dist_rulers, 2, communication_mode);
			std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
			
			test::regular_test(comm, s, d);

		}
		else if (ruling_set2_rec.compare(argv[1]) == 0)
		{
			
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			regular_ruling_set2 algorithm(s, dist_rulers, 2);
			std::vector<std::int64_t> d = algorithm.start(comm);
			
			test::regular_test(comm, s, d);
		}
		else if (pointer_doubling.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);

			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			regular_pointer_doubling algorithm(s, comm);
			
			std::vector<std::int64_t> d = algorithm.start(comm, grid_comm);
			
			test::regular_test(comm, s, d);
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
			std::vector<std::uint64_t> tree_vector = generator::generate_regular_successor_vector(num_local_vertices, comm);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::int64_t> d = forest_regular_ruling_set2(tree_vector, dist_rulers, comm,1).result_dist;
			
			//analyze_instances::analyze_regular_instance(tree_vector, comm);
			test::regular_test(comm, tree_vector, d);
		}
		else if (tree_rooting_rec.compare(argv[1]) == 0)
		{
			std::int32_t num_local_vertices = atoi(argv[2]);
			std::vector<std::uint64_t> tree_vector = generator::generate_regular_wood_vector(num_local_vertices, comm);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::int64_t> d = forest_regular_ruling_set2(tree_vector, dist_rulers, comm,2).result_dist;
			
			//analyze_instances::analyze_regular_instance(tree_vector, comm);
			test::regular_test(comm, tree_vector, d);
		}
		else if (euler_tour.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> tree_vector = generator::generate_regular_tree_vector(num_local_vertices, comm);

			std::vector<std::int64_t> d = tree_euler_tour(comm, tree_vector, dist_rulers).start(comm, tree_vector);
			test::regular_test(comm, tree_vector, d);
		}
		else if (lb_tree_rooting.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_tree_vector(num_local_vertices, comm);

			//forest_load_balance_regular_ruling_set2(dist_rulers).start2(s,comm);
			//real_load_balance(dist_rulers).start(s,comm, grid_comm);
			std::vector<std::int64_t> d = shuffle_load_balance().start(s,comm, grid_comm, dist_rulers);
			test::regular_test(comm, s, d);

			//analyze_instances::analyze_regular_instance(s, comm);
		}
		else if (local_contract.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_wood_vector(num_local_vertices, comm);
			forest_local_contraction2 algorithm;
			std::vector<std::int64_t> d = algorithm.start(comm, s);
			
			test::regular_test(comm, s, d);
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
			test::regular_test(comm, s, d);
			
			//analyze_instances::analyze_regular_instance(s, comm);

			
		}
		else if (asynchron.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			example example;
			
			std::vector<std::int64_t> d = example.test(s, dist_rulers, comm);
			
			test::regular_test(comm, s, d);
		}
		else if (asynchron2.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			asynchron_ruling_set2 algorithm;
			
			std::vector<std::int64_t> d = algorithm.test(s, dist_rulers, comm);
			
			test::regular_test(comm, s, d);
		}
		else if (grid_asynchron.compare(argv[1]) == 0)
		{
			std::uint64_t num_local_vertices = atoi(argv[2]);
			std::int32_t dist_rulers = atoi(argv[3]);
			std::vector<std::uint64_t> s = generator::generate_regular_successor_vector(num_local_vertices, comm);
			
			grid_asynchron_ruling_set2 algorithm;
			
			std::vector<std::int64_t> d = algorithm.test(s, dist_rulers, comm);
			
			test::regular_test(comm, s, d);
		}
		else
		{
			/*
			int rank, size;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			CLI::App app;
			app.option_defaults()->always_capture_default();
			int queue_version = 1;
			app.add_option("--queue_version", queue_version)->expected(1, 2);
			std::size_t number_of_messages = 10;
			app.add_option("--number_of_messages", number_of_messages, "The number of messages to send from each PE");
			std::size_t iterations = 1;
			app.add_option("--iterations", iterations);
			bool use_test_any = false;
			app.add_flag("--use_test_any", use_test_any);
			bool use_custom_implementations = false;
			app.add_flag("--use_custom_implementations", use_custom_implementations);

			//CLI11_PARSE(app, argc, argv);

			DEBUG_BARRIER(rank);
			const int message_size = 10;

			auto merger = [](std::vector<int>& buffer, std::vector<int> msg, int) {
				for (auto elem : msg) {
					buffer.emplace_back(elem);
				}
				return msg.size();
			};
			auto splitter = [](std::vector<int>& buffer, auto on_message, message_queue::PEID sender) {
				for (size_t i = 0; i < buffer.size(); i += message_size) {
					on_message(buffer.cbegin() + i, buffer.cbegin() + i + message_size, sender);
				}
			};
			for (size_t i = 0; i < iterations; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				double start = MPI_Wtime();
				int local_max_test_size = 0;
				size_t local_max_active_requests = 0;

				using MessageContainer = std::vector<int>;
				// using MessageContainer = message_queue::testing::OwnContainer<int>;
				// using MessageContainer = std::vector<int, message_queue::testing::CustomAllocator<int>>;

				auto queue = message_queue::MessageQueue<int, MessageContainer>{};
				if (use_test_any) {
					queue.use_test_any();
				}
				std::default_random_engine eng;
				eng.seed(rank);
				std::bernoulli_distribution bernoulli_dist(0.1);
				std::uniform_int_distribution<size_t> rank_dist(1, size - 1);
				// auto queue = message_queue::make_mesqueue<int>(std::move(merger), std::move(splitter));
				// queue.set_threshold(200);
				message_queue::PEID receiver = rank_dist(eng);
				MessageContainer message(message_size);
				message[0] = rank;
				message[1] = 0;
				for (size_t i = 0; i < number_of_messages; ++i) {
					message[2] = i;
					queue.post_message(MessageContainer(message), (rank + rank_dist(eng)) % size);
				}
				auto on_message = [&](message_queue::Envelope<int> auto envelope) {
					if (bernoulli_dist(eng)) {
						auto begin = envelope.message.begin();
						std::stringstream ss;
						ss << "Message " << *(begin + 2) << " from " << *begin << " arrived after " << *(begin + 1) << " hops.";
						message_queue::atomic_debug(ss.str());
					} else {
						KASSERT(envelope.message.size() > 1);
						envelope.message[1]++;
						queue.post_message(std::move(envelope.message), (rank + rank_dist(eng)) % size);
					}
				};
				queue.poll(on_message);
				queue.terminate(on_message);
				using namespace std::chrono;
				local_max_active_requests = queue.max_active_requests();
				MPI_Barrier(MPI_COMM_WORLD);
				double end = MPI_Wtime();
				int global_max_test_size;
				MPI_Reduce(&local_max_test_size, &global_max_test_size, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

				size_t global_max_active_requests;
				MPI_Reduce(&local_max_active_requests, &global_max_active_requests, 1, MPI_UINT64_T, MPI_MAX, 0,
						   MPI_COMM_WORLD);
				// print CLI options
				std::unordered_map<std::string, std::string> stats;
				for (const auto& option : app.get_options()) {
					if (option->get_single_name() == "help") {
						continue;
					}
					stats[option->get_single_name()] = option->as<std::string>();
					if (stats[option->get_single_name()].empty()) {
						stats[option->get_single_name()] = "false";
					}
				}
				stats["ranks"] = fmt::format("{}", size);
				stats["time"] = fmt::format("{}", end - start);
				stats["iteration"] = fmt::format("{}", i);
				stats["max_test_size"] = fmt::format("{}", global_max_test_size);
				stats["max_active_requests"] = fmt::format("{}", global_max_active_requests);

				if (rank == 0) {
					std::cout << "RESULT";
					for (const auto& [key, value] : stats) {
						std::cout << " " << key << "=" << value;
					}
					std::cout << "\n";
				}
			}*/
			
			error(std::string(argv[1]) + " is not a name of an algorithm or wrong parameters");
		}
	}
	
	
	MPI_Finalize();
	return 0;
	
};

