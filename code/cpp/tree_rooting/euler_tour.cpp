class euler_tour2
{

	public:
	
	euler_tour2(kamping::Communicator<>& comm, std::uint64_t num_local_vertices_param)
	{
		
		num_local_vertices = num_local_vertices_param;
		size = comm.size();
		rank = comm.rank();
		node_offset = rank * num_local_vertices;
		
		generate_regular_tree_vector(comm);
		
		calculate_euler_tour(comm);
	}
	
	
	void calculate_euler_tour(kamping::Communicator<>& comm)
	{
		num_local_edges = all_edges.size();
		std::vector<std::int32_t> send(1,num_local_edges);
		std::vector<std::int32_t> recv_sizes;
		comm.allgather(kamping::send_buf(send), kamping::recv_buf(recv_sizes));
		
		std::vector<std::uint64_t> prefix_sum_num_edges_per_PE(size + 1,0);
		for (std::uint64_t i = 1; i < size + 1; i++)
			prefix_sum_num_edges_per_PE[i] = recv_sizes[i-1] + prefix_sum_num_edges_per_PE[i-1];
		
		
		struct edge{
			std::uint64_t source;
			std::uint64_t destination;
		};
		
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		for (std::uint64_t j = bounds[i]; j < bounds[i+1]; j++)
		{
			std::int32_t targetPE = all_edges[j] / num_local_vertices;
			num_packets_per_PE[targetPE]++;	
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		std::vector<edge> request(send_displacements[size]);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		for (std::uint64_t j = bounds[i]; j < bounds[i+1]; j++)
		{
			std::int32_t targetPE = all_edges[j] / num_local_vertices;
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			
			request[packet_index].source = i + node_offset;
			request[packet_index].destination = all_edges[j];
			
		}
			
		auto recv_request = comm.alltoallv(kamping::send_buf(request), kamping::send_counts(num_packets_per_PE));
		std::vector<edge> recv_buffer = recv_request.extract_recv_buffer();
		std::vector<std::uint64_t> answers(recv_buffer.size());
		
		for (std::uint64_t i = 0; i < recv_buffer.size(); i++)
		{
			std::uint64_t lower_bound = bounds[recv_buffer[i].destination - node_offset];
			std::uint64_t upper_bound = bounds[recv_buffer[i].destination - node_offset + 1];
			
			for (std::uint64_t j = lower_bound; j < upper_bound; j++)
				if (all_edges[j] == recv_buffer[i].source)
				{
					answers[i] = prefix_sum_num_edges_per_PE[rank] + lower_bound + ((j - lower_bound + 1) % (upper_bound - lower_bound));
					
					break;
				}
			
		}
		
		auto recv_answers = comm.alltoallv(kamping::send_buf(answers), kamping::send_counts(recv_request.extract_recv_counts())).extract_recv_buffer();
		std::vector<std::uint64_t> s(num_local_edges);
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);

		
		
			
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		for (std::uint64_t j = bounds[i]; j < bounds[i+1]; j++)
		{
			std::int32_t targetPE = all_edges[j] / num_local_vertices;
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			
			s[j] = recv_answers[packet_index];
			
			
		}
		
		
		
		std::vector<std::uint64_t> r(num_local_edges,1);
		if (rank == 0)
		{	
			s[0] = 0;
			r[0] = 0;
		}
		std::vector<std::uint32_t> targetPEs(num_local_edges);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		for (std::uint64_t j = bounds[i]; j < bounds[i+1]; j++)
			targetPEs[j] = all_edges[j] / num_local_vertices; 
		
		
		irregular_pointer_doubling algorithm(s,r,targetPEs,prefix_sum_num_edges_per_PE);
		std::vector<std::uint64_t> ranks = algorithm.start(comm);
		
		std::cout << "PE " << rank << " with rank array:";
		for (int i = 0; i < num_local_edges; i++)
			std::cout << ranks[i] << ",";
		std::cout << std::endl;
		
		
	}
	
	std::uint64_t hash64(std::uint64_t x) {
		x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
		x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
		x = x ^ (x >> 31);
		return x;
	}

	void generate_regular_tree_vector(kamping::Communicator<>& comm)
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
		
	
		struct edge{
			std::uint64_t source;
			std::uint64_t destination;
		};
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		for (std::uint64_t i = 0; i < s.size(); i++)
		{
			if (i + node_offset == s[i])
				continue;
			
			std::int32_t targetPE = s[i] / s.size();
			num_packets_per_PE[targetPE]++;		
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		std::vector<edge> edges(s.size());
		
		for (std::uint64_t i = 0; i < s.size(); i++)
		{
			if (i + node_offset == s[i])
				continue;
			
			std::int32_t targetPE = s[i] / s.size();
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			
			edges[packet_index].source = i + node_offset;
			edges[packet_index].destination = s[i];
		}
		
	

		
		auto recv = comm.alltoallv(kamping::send_buf(edges), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
	
		all_edges = std::vector<std::uint64_t>(recv.size() + s.size());
		std::vector<std::uint64_t> edges_per_node(s.size(),0);
		
		for (std::uint64_t i = 0; i < recv.size(); i++)
		{
			edges_per_node[recv[i].destination - node_offset]++;
		}
		for (std::uint64_t i = 0; i < s.size(); i++)
			edges_per_node[i]++; //for only not turned edge, since it was originally a tree
		
		bounds = std::vector<std::uint64_t>(s.size() + 1, 0);
		for (std::uint64_t i = 1; i <= s.size(); i++)
			bounds[i] = bounds[i-1] + edges_per_node[i-1];
		std::fill(edges_per_node.begin(), edges_per_node.end(), 0);
		for (std::uint64_t i = 0; i < s.size(); i++)
		{
			std::uint64_t packet_index = bounds[i] + edges_per_node[i]++;
			all_edges[packet_index] = s[i];
		}
		for (std::uint64_t i = 0; i < recv.size(); i++)
		{
			std::uint64_t target_node = recv[i].destination - node_offset;
			std::uint64_t packet_index = bounds[target_node] + edges_per_node[target_node]++;
			all_edges[packet_index] = recv[i].source;
		}
		
		
		/*
		//print for testing
		std::cout << "s: ";
		for (int i = 0; i < s.size(); i++)
			std::cout << "(" << i + node_offset << "," << s[i] << "),";
		std::cout << "\nbounds: ";
		for (int i = 0; i < bounds.size(); i++)
			std::cout << bounds[i] << ",";
		std::cout << "\nturned edges: ";
		for (int i = 0; i < all_edges.size(); i++)
			std::cout << all_edges[i] << ",";
		std::cout << std::endl;
		*/
	
	}
	
	void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t	>& num_packets_per_PE)
	{
		send_displacements[0]=0;
		for (std::int32_t i = 1; i < size + 1; i++)
			send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
	}
	
	
	private: 
	
	
	
	std::int32_t size;
	std::int32_t rank;
	std::int32_t num_local_vertices;
	std::int32_t node_offset;
	
	std::uint64_t num_local_edges;
	
	std::vector<std::uint64_t> all_edges;
	std::vector<std::uint64_t> bounds;
	
};