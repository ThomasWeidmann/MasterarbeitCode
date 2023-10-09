/*
This class calculates from a given sucessor array of a regular tree with the euler tour tree rooting
*/
class tree_euler_tour
{

	public:
	
	tree_euler_tour(kamping::Communicator<>& comm, std::vector<std::uint64_t>& s)
	{

		num_local_vertices = s.size();
		size = comm.size();
		rank = comm.rank();
		node_offset = rank * num_local_vertices;
		
		
		tree_rooting(comm, s);
	}
	
	std::vector<std::uint64_t> tree_rooting(kamping::Communicator<>& comm, std::vector<std::uint64_t>& s)
	{
		std::int64_t local_root = -1;
		
		//first calculate adjacency array by also turning edges
		struct edge{
			std::uint64_t source;
			std::uint64_t destination;
		};
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		for (std::uint64_t i = 0; i < s.size(); i++)
		{
			if (i + node_offset == s[i])
			{
				local_root = i + node_offset;
				continue;
			}
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
		all_edges_weights = std::vector<std::int64_t>(recv.size() + s.size());
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
			all_edges_weights[packet_index] = 1;
		}
		for (std::uint64_t i = 0; i < recv.size(); i++)
		{
			std::uint64_t target_node = recv[i].destination - node_offset;
			std::uint64_t packet_index = bounds[target_node] + edges_per_node[target_node]++;
			all_edges[packet_index] = recv[i].source;
			all_edges_weights[packet_index] = -1;
		}
		
		// graph sucessfully added all back edges to adjacency array
		num_local_edges = all_edges.size();
		std::vector<std::uint64_t> send(1,num_local_edges);
		std::vector<std::uint64_t> recv_sizes;
		comm.allgather(kamping::send_buf(send), kamping::recv_buf(recv_sizes));
		
		std::vector<std::uint64_t> prefix_sum_num_edges_per_PE(size + 1,0);
		for (std::uint64_t i = 1; i < size + 1; i++)
			prefix_sum_num_edges_per_PE[i] = recv_sizes[i-1] + prefix_sum_num_edges_per_PE[i-1];
		
		//idea, if we want the following edge of (a,b) in the euler tour, we need to fing the edge (b,a) in adjacency array to answer this
		//so we sort adjacency array to perfom binary sorts, so for every node i the values all_edges[bounds[i]], ..., all_edges[bounds[i+1]-1] are in ascending order
		
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
			std::sort(&all_edges[bounds[i]], &all_edges[bounds[i+1]]);
		
	
		
		
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);

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
			
			std::uint64_t dynamic_lower_bound = lower_bound;
			std::uint64_t dynamic_upper_bound = upper_bound;
			
			
			while (dynamic_upper_bound - dynamic_lower_bound > 1)
			{
				std::uint64_t middle = (dynamic_lower_bound + dynamic_upper_bound) / 2;
				
				if (all_edges[middle] > recv_buffer[i].source)
					dynamic_upper_bound = middle;
				else
					dynamic_lower_bound = middle;
				
			}
			std::uint64_t result = dynamic_lower_bound; //all_edges[result] == recv_buffer[i].source
			
			
			
			answers[i] = prefix_sum_num_edges_per_PE[rank] + lower_bound + ((result - lower_bound + 1) % (upper_bound - lower_bound));
			
			
		}
		
		
		
		auto recv_answers = comm.alltoallv(kamping::send_buf(answers), kamping::send_counts(recv_request.extract_recv_counts())).extract_recv_buffer();
		std::vector<std::uint64_t> s_edges(num_local_edges);
	
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		for (std::uint64_t j = bounds[i]; j < bounds[i+1]; j++)
		{
			std::int32_t targetPE = all_edges[j] / num_local_vertices;
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			
			s_edges[j] = recv_answers[packet_index];
			
			
		}
		
		
		if (local_root != -1)
		{
			std::uint64_t dynamic_lower_bound = bounds[local_root];
			std::uint64_t dynamic_upper_bound = bounds[local_root+1];
			
			
			while (dynamic_upper_bound - dynamic_lower_bound > 1)
			{
				std::uint64_t middle = (dynamic_lower_bound + dynamic_upper_bound) / 2;
				
				if (all_edges[middle] > local_root + prefix_sum_num_edges_per_PE[rank])
					dynamic_upper_bound = middle;
				else
					dynamic_lower_bound = middle;
				
			}
			std::uint64_t result = dynamic_lower_bound;
			
			s_edges[result] = local_root + prefix_sum_num_edges_per_PE[rank];
			all_edges_weights[result] = 0;
		}
		
		
		std::cout << "PE " << rank << " with s arr: ";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << s[i] << ",";
		std::cout << std::endl;
		std::cout << "PE " << rank << " with bounds arr: ";
		for (int i = 0; i < bounds.size(); i++)
			std::cout << bounds[i] << ",";
		std::cout << std::endl;
		std::cout << "PE " << rank << " with all_edges arr: ";
		for (int i = 0; i < all_edges.size(); i++)
			std::cout << all_edges[i] << ",";
		std::cout << std::endl;
		std::cout << "PE " << rank << " with s_edges arr: ";
		for (int i = 0; i < s_edges.size(); i++)
			std::cout << s_edges[i] << ",";
		std::cout << std::endl;
		std::cout << "PE " << rank << " with all_edges_weights arr: ";
		for (int i = 0; i < all_edges_weights.size(); i++)
			std::cout << all_edges_weights[i] << ",";
		std::cout << std::endl;
		
		std::vector<std::uint32_t> targetPEs(num_local_edges);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		for (std::uint64_t j = bounds[i]; j < bounds[i+1]; j++)
			targetPEs[j] = all_edges[j] / num_local_vertices; 
		
		/*
		irregular_pointer_doubling algorithm(s_edges,all_edges_weights,targetPEs,prefix_sum_num_edges_per_PE);
		std::vector<std::uint64_t> ranks = algorithm.start(comm);
		
		std::cout << "PE " << rank << " with rank array:";
		for (int i = 0; i < num_local_edges; i++)
			std::cout << ranks[i] << ",";
		std::cout << std::endl;*/
		
		return all_edges;
		
		
		
	}
	
	
	void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t>& num_packets_per_PE)
	{
		send_displacements[0]=0;
		for (std::int32_t i = 1; i < size + 1; i++)
			send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
	}
	
	
	private: 
	
	
	
	std::uint64_t size;
	std::uint64_t rank;
	std::uint64_t num_local_vertices;
	std::uint64_t node_offset;
	
	std::uint64_t num_local_edges;
	
	std::vector<std::uint64_t> all_edges;
	std::vector<std::int64_t> all_edges_weights;
	std::vector<std::uint64_t> bounds;
	
};