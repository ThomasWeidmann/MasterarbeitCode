class forest_load_balance_regular_ruling_set2 
{
	
	
	public:
	
	forest_load_balance_regular_ruling_set2(std::uint64_t comm_rounds)
	{
		this-> comm_rounds = comm_rounds;
	}
	
	void start(std::vector<std::uint64_t>& s, kamping::Communicator<>& comm)
	{
		size = comm.size();
		rank = comm.rank();
		num_local_vertices = s.size();
		node_offset = rank * num_local_vertices;
		
		
		//als allererstes kanten umdrehen, aber so load balancen, sodass jeder PE gleich viele kanten hat
		//dabei können auch knoten auf verschiedene PEs aufgeteilt werden
	
		//targetPE wird struct mit targetPE_lower_bound <= targetPE < targetPE_upper_bound
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		std::unordered_map<std::uint64_t, std::int64_t> local_node_indegrees;

		std::vector<std::uint64_t> num_edges_per_PE(size,0);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{	
			if (local_node_indegrees.contains(s[i]))
				local_node_indegrees[s[i]] = local_node_indegrees[s[i]] + 1;
			else
				local_node_indegrees[s[i]] = 1;
		}
		
		
		for (const auto& [key, value] : local_node_indegrees)
		{
			std::int32_t targetPE = key / num_local_vertices;
			num_packets_per_PE[targetPE]++;
		}
		
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		struct packet {
			std::uint64_t node;
			std::uint64_t indegree;
		};
		std::vector<packet> send_packets(send_displacements[size]);
		for (const auto& [key, value] : local_node_indegrees)
		{
			std::int32_t targetPE = key / num_local_vertices;
			std::int64_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			send_packets[packet_index].node = key;
			send_packets[packet_index].indegree = value;
		}
		
		std::vector<packet> recv_packets = comm.alltoallv(kamping::send_buf(send_packets), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
	
		std::vector<std::uint64_t> indegrees(num_local_vertices,0);
		std::uint64_t weight = num_local_vertices;
		for (int i = 0; i < recv_packets.size(); i++)
		{
			weight += recv_packets[i].indegree;
			indegrees[recv_packets[i].node - node_offset] += recv_packets[i].indegree;
		}
		
		std::vector<std::uint64_t> all_weights;
		comm.allgather(kamping::send_buf(weight), kamping::recv_buf<kamping::resize_to_fit>(all_weights));
		std::vector<std::uint64_t> prefix_sum_all_weights(size+1,0);
		for (std::uint32_t i = 1; i < size + 1; i++)
		{
			prefix_sum_all_weights[i] = prefix_sum_all_weights[i-1] + all_weights[i-1];
		}
		
		
/*
		//every PE needs to know the PEs from which we split, so 2 binary searches on all_weights
		std::uint64_t start_weight = 2 * rank * num_local_vertices;
		std::uint64_t end_weight = 2 * (rank + 1) * num_local_vertices;
		
		std::uint32_t lower_bound_PE=0;
		std::uint32_t upper_bound_PE =size;
		
		while (upper_bound_PE - lower_bound_PE > 1)
		{
			std::uint32_t middle = (lower_bound_PE + upper_bound_PE) / 2;
			if (all_weights[middle]  > start_weight)
				upper_bound_PE = middle;
			else
				lower_bound_PE = middle;
		}
		std::uint32_t start_weight_PE = lower_bound_PE;
		
		lower_bound_PE=0;
		upper_bound_PE =size;
		
		while (upper_bound_PE - lower_bound_PE > 1)
		{
			std::uint32_t middle = (lower_bound_PE + upper_bound_PE) / 2;
			if (all_weights[middle]  > end_weight)
				upper_bound_PE = middle;
			else
				lower_bound_PE = middle;
		}
		std::uint32_t end_weight_PE = lower_bound_PE;
		
		if (rank == 0)
		{
			std::cout << "prefix sum weight array:[";
			for (int i = 0; i < size+1; i++)
				std::cout << prefix_sum_all_weights[i] << ",";
			std::cout << "]" << std::endl;
		}*/
		
		struct node_range {
			std::uint64_t start_index;
			std::uint64_t end_index;
		};
		
		std::vector<node_range> local_node_ranges(0);
		
		std::uint64_t dynamic_start_weight = prefix_sum_all_weights[rank];
		while (dynamic_start_weight < prefix_sum_all_weights[rank+1])
		{
			if (dynamic_start_weight % (2* num_local_vertices) == 0)
			{
				if (prefix_sum_all_weights[rank+1]-dynamic_start_weight > 2 *  num_local_vertices)
				{
					local_node_ranges.push_back({dynamic_start_weight, dynamic_start_weight + 2*num_local_vertices});
					dynamic_start_weight += 2*num_local_vertices;
				}
				else
				{
					local_node_ranges.push_back({dynamic_start_weight, prefix_sum_all_weights[rank+1]});
					break;
				}
			}
			else
			{
				std::uint64_t dynamic_end_weight = ((dynamic_start_weight + (2*num_local_vertices-1))/(2*num_local_vertices))*2*num_local_vertices;
				if (dynamic_end_weight < prefix_sum_all_weights[rank+1])
				{
					local_node_ranges.push_back({dynamic_start_weight, dynamic_end_weight});
				}
				else
				{
					local_node_ranges.push_back({dynamic_start_weight, prefix_sum_all_weights[rank+1]});
				}
				dynamic_start_weight = dynamic_end_weight;
			}
			
		}
		
		std::cout << rank << " with node_ranges:";
		for (int i = 0; i < local_node_ranges.size(); i++)
			std::cout << "[" << local_node_ranges[i].start_index << "," << local_node_ranges[i].end_index << "],";
		std::cout << std::endl;
		
		struct node_info_for_node_range {
			std::uint64_t start_node;
			std::uint64_t start_node_start_degree;
			std::uint64_t end_node;
			std::uint64_t end_node_end_degree; //this means end_node_end_degree = 1 then we only want of 1 edge of end_node 
		};
		
		std::vector<node_info_for_node_range> node_info_for_node_ranges(local_node_ranges.size());
		dynamic_start_weight = prefix_sum_all_weights[rank];
		std::uint64_t dynamic_start_node = node_offset;
		std::uint64_t dynamic_node_degree = 0;
		for (int i = 0; i < local_node_ranges.size(); i++)
		{	
			node_info_for_node_range info;
			info.start_node = dynamic_start_node;
			info.start_node_start_degree = dynamic_node_degree;
			
			
		}
		
		
		
	
	
	}
	
	void start2(std::vector<std::uint64_t>& s, kamping::Communicator<>& comm)
	{
		size = comm.size();
		rank = comm.rank();
		num_local_vertices = s.size();
		node_offset = rank * num_local_vertices;
		
		
		//als allererstes kanten umdrehen, aber so load balancen, sodass jeder PE gleich viele kanten hat
		//dabei können auch knoten auf verschiedene PEs aufgeteilt werden
	
		//targetPE wird struct mit targetPE_lower_bound <= targetPE < targetPE_upper_bound
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		std::vector<std::uint64_t> num_edges_per_PE(size,0);
		for (std::uint64_t i = 0; i < s.size(); i++)
		{	
			std::int32_t targetPE = s[i] / num_local_vertices;
			num_edges_per_PE[targetPE]++;		
		}
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 1);

		auto recv_edges_per_PE = comm.alltoallv(kamping::send_buf(num_edges_per_PE), kamping::send_counts(num_packets_per_PE));
		std::vector<std::uint64_t> recv_num_edges_per_PE = recv_edges_per_PE.extract_recv_buffer();
		std::uint64_t num_local_edges = 0;
		for (std::uint32_t i = 0; i < recv_num_edges_per_PE.size(); i++)
			num_local_edges += recv_num_edges_per_PE[i];

		
		comm.allgather(kamping::send_buf(num_local_edges), kamping::recv_buf<kamping::resize_to_fit>(num_edges_per_PE));
		std::vector<std::uint64_t> prefix_sum_weight_per_PE(size+1,0);
		for (std::uint32_t i = 1; i < size +1; i++)
		{
			prefix_sum_weight_per_PE[i] = prefix_sum_weight_per_PE[i-1] + num_edges_per_PE[i-1] + num_local_vertices;
		}
		if (rank == 0)
		{
			std::cout << "weights per PE: [";
			for (int i = 0; i < size+1; i++)
				std::cout << prefix_sum_weight_per_PE[i] << ",";
			std::cout << "]" << std::endl;
		}
		//an edge has weight 1 and an node has weight 1, so the whole list has weight 2n and we want every PE to have 2n/p of weight
		std::uint64_t start_weight = 2 * rank * num_local_vertices;
		std::uint64_t end_weight = 2 * (rank + 1) * num_local_vertices;
		//so every PE has nodes and edges that have start_weight<= weight < end_weight
		
		//first find start_weight
		std::uint32_t lower_bound_PE=0;
		std::uint32_t upper_bound_PE =size;
		
		while (upper_bound_PE - lower_bound_PE > 1)
		{
			std::uint32_t middle = (lower_bound_PE + upper_bound_PE) / 2;
			if (prefix_sum_weight_per_PE[middle]  > start_weight)
				upper_bound_PE = middle;
			else
				lower_bound_PE = middle;
		}
		std::uint32_t start_weight_PE = lower_bound_PE;


		lower_bound_PE=0;
		upper_bound_PE =size;
		
		while (upper_bound_PE - lower_bound_PE > 1)
		{
			std::uint32_t middle = (lower_bound_PE + upper_bound_PE) / 2;
			if (prefix_sum_weight_per_PE[middle]  > end_weight)
				upper_bound_PE = middle;
			else
				lower_bound_PE = middle;
		}
		std::uint32_t end_weight_PE = lower_bound_PE;
		std::cout << rank << " with PE arr [" << start_weight_PE << "," << lower_bound_PE << "]" << std::endl;
		
		for (std::uint32_t i = start_weight_PE; i <= end_weight_PE; i++)
		{
			if (i == start_weight_PE && i == end_weight_PE)
				std::cout << rank << " will von " << i << " alle nodes mit " << start_weight << " <= weight < " << end_weight << std::endl; 
			else if (i == start_weight_PE)
				std::cout << rank << " will von " << start_weight_PE << " alle nodes mit weight >= " << start_weight << std::endl;
			else if (i == end_weight_PE)
				std::cout << rank << " will von " << end_weight_PE << " alle nodes mid weight < " <<  end_weight << std::endl;
			else
				std::cout << rank << " will von " << i << " alle nodes " << std::endl;
				
			
		}
		
		struct nodes_request {
			std::uint64_t lower_bound_weights;
			std::uint64_t upper_bound_weights;
		}; //ein node requested alle nodes die lower_bound_weights <= weight < upper_bound_weights;
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);

		std::vector<nodes_request> nodes_request(0);
		
		for (std::uint32_t targetPE = start_weight_PE; targetPE <= end_weight_PE; targetPE++)
		{
			num_packets_per_PE[targetPE]++;
			
			if (targetPE == start_weight_PE && targetPE == end_weight_PE)
				nodes_request.push_back({start_weight, end_weight});
			else if (targetPE == start_weight_PE)
				nodes_request.push_back({start_weight, prefix_sum_weight_per_PE[targetPE+1]});
			else if (targetPE == end_weight_PE)
				nodes_request.push_back({prefix_sum_weight_per_PE[targetPE], end_weight});
			else
				nodes_request.push_back({prefix_sum_weight_per_PE[targetPE], prefix_sum_weight_per_PE[targetPE+1]});
			
		}
		//this is just sparse all to all, one PE sends to at most 3 other PEs
		auto recv_nodes_request = comm.alltoallv(kamping::send_buf(nodes_request), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
		
		std::cout << rank << " with recv:[";
		for (int i = 0; i < recv_nodes_request.size(); i++)
			std::cout << "(" << recv_nodes_request[i].lower_bound_weights << "," << recv_nodes_request[i].upper_bound_weights << "),";
		std::cout << "]" << std::endl;
		
		std::vector<std::uint64_t> bounds(1,rank);
		std::vector<std::uint64_t> edges(1,rank);
		
		std::uint64_t current_weight = prefix_sum_weight_per_PE[rank];
		for (int i = 0; i < recv_nodes_request.size(); i++)
		{
			
		}
		
	}
	
	void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t>& num_packets_per_PE)
	{
		send_displacements[0]=0;
		for (std::int32_t i = 1; i < size + 1; i++)
			send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
	}
	
	
	private:
	std::uint64_t comm_rounds;
	std::uint64_t node_offset;
	std::uint64_t num_local_vertices;
	std::uint64_t rank, size;
};