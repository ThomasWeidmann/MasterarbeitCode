class generator
{
	public:
	

	static std::vector<std::uint64_t> generate_regular_tree_vector(std::uint64_t num_local_vertices, kamping::Communicator<>& comm)
	{
		std::uint64_t mpi_rank = comm.rank();
		std::uint64_t mpi_size = comm.size();
		
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

	static std::vector<std::uint64_t> generate_regular_wood_vector(std::uint64_t num_local_vertices, kamping::Communicator<>& comm)
	{
		std::uint64_t mpi_rank = comm.rank();
		std::uint64_t mpi_size = comm.size();
		
		std::int32_t max_edge_weight = 150; 
		std::uint64_t node_offset = mpi_rank * num_local_vertices;
		kagen::KaGen gen(MPI_COMM_WORLD);
		
		double prob = std::log(num_local_vertices) / ((double) num_local_vertices * mpi_size);
		prob = prob > 1 ? 1 : prob;
		prob = 30 / ((double) num_local_vertices * mpi_size);
		prob = prob > 1 ? 1 : prob;
		//GNM
		//auto graph = gen.GenerateUndirectedGNP(num_local_vertices * mpi_size, prob, false);
		auto graph = gen.GenerateUndirectedGNM(num_local_vertices * mpi_size, 20 * num_local_vertices * mpi_size, false);
		
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


	
	static std::vector<std::uint64_t> generate_regular_successor_vector(std::uint64_t num_local_vertices, kamping::Communicator<>& comm)
	{
		
		kagen::KaGen gen(MPI_COMM_WORLD);
		std::vector<std::uint64_t> s(num_local_vertices);
		std::uint64_t num_global_vertices = comm.size() * num_local_vertices;
		auto path = gen.GenerateDirectedPath(num_global_vertices, true);
		
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
			s[i] = i + comm.rank() * num_local_vertices;
		
		for (auto const& [src, dst]: path.edges)
			s[src - comm.rank() * num_local_vertices] = dst;
		

		return s;
	}
	
	private:
	
	static std::uint64_t hash64(std::uint64_t x) {
		x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
		x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
		x = x ^ (x >> 31);
		return x;
	}
	
	std::uint64_t rank;
	std::uint64_t size;
	
	
};