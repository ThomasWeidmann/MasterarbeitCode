class forest_load_balance_regular_ruling_set2 
{
	
	
	public:
	
	forest_load_balance_regular_ruling_set2(std::uint64_t comm_rounds)
	{
		this-> comm_rounds = comm_rounds;
	}
	
	//hierbei werden, wenn ein node einen zu hohen indegree hat, pseudonodes hinzugefügt
	void start2(std::vector<std::uint64_t>& s, kamping::Communicator<>& comm)
	{
		size = comm.size();
		rank = comm.rank();
		num_local_vertices = s.size();
		node_offset = rank * num_local_vertices;
		
		
		//wir wollen kanten mit hohen eingangsgraden in mehrere nodes aufsplitten, sodass kein PE mehr als doppelt so viele Kanten hat, wie bei einer komlett random eingabe
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
		std::uint64_t sum_of_indegrees = 0;
		for (int i = 0; i < recv_packets.size(); i++)
		{
			sum_of_indegrees += recv_packets[i].indegree;
			indegrees[recv_packets[i].node - node_offset] += recv_packets[i].indegree;
		}
		
		struct node {
			std::uint64_t global_index;
			std::uint64_t indegree;
		};
		std::vector<node> nodes(num_local_vertices);
		for (int i = 0; i < num_local_vertices; i++)
		{
			nodes[i].global_index = i + node_offset;
			nodes[i].indegree = indegrees[i];
		}
		
		
		std::sort(nodes.begin(), nodes.end(), [](node  a, node b) { return a.indegree > b.indegree; });
		//so we want sum_of_indegrees <= 2*num_local_vertices
		//so we iterate over nodes with highest indegrees and outsource them until we reached out goal for sum_of_indegrees

		std::vector<node> nodes_to_outsource(0);
		std::int64_t sum_of_indegrees_to_outsource = sum_of_indegrees - 2*num_local_vertices;
		if (sum_of_indegrees_to_outsource > 0)
		{
			std::uint64_t outsourced_indegrees = 0;
			std::uint64_t index_of_nodes_to_outsource = 0;
			while (sum_of_indegrees_to_outsource > outsourced_indegrees)
			{
				nodes_to_outsource.push_back(nodes[index_of_nodes_to_outsource]);
				outsourced_indegrees += nodes[index_of_nodes_to_outsource].indegree;
				index_of_nodes_to_outsource++;
			}
			
		}
		
	}
	
	
	
	//hierbei werden die weight vectoren aufgeteilt und keine neuen vertices dazugemacht
	//worst case O(p^2 + n/p) space and time
	void start(std::vector<std::uint64_t>& s, kamping::Communicator<>& comm)
	{
		size = comm.size();
		rank = comm.rank();
		num_local_vertices = s.size();
		node_offset = rank * num_local_vertices;
		
		/*
		std::cout << rank << " with s array:";
		for (int i = 0; i < s.size(); i++)
			std::cout << s[i] << " ";
		std::cout << std::endl;*/
		
		
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
		
		auto recv = comm.alltoallv(kamping::send_buf(send_packets), kamping::send_counts(num_packets_per_PE));
		std::vector<packet> recv_packets = recv.extract_recv_buffer();
		std::vector<std::int32_t> recv_counts = recv.extract_recv_counts();
		std::vector<std::int32_t> recv_displs = recv.extract_recv_displs();
			
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
		std::cout << rank << " with indegrees:";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << indegrees[i] << " ";
		std::cout << std::endl;*/
		
		
		/*
		struct node_range {
			std::uint64_t start_index;
			std::uint64_t end_index;
		};
		
		struct node_info_for_node_range { //if any degree is -1 then that means every edge of the node
			std::uint64_t start_node;
			std::int64_t start_node_start_degree;
			std::uint64_t end_node;
			std::int64_t end_node_end_degree; //this means when end_node_end_degree = 1 then we only want of 1 edge of end_node 
		};
		
		std::vector<node_range> local_node_ranges(0);
		std::vector<node_info_for_node_range> local_node_info_for_node_range(0);
		
		std::uint64_t dynamic_start_weight = prefix_sum_all_weights[rank];
		while (dynamic_start_weight < prefix_sum_all_weights[rank+1])
		{
			node_info_for_node_range info;
			if (dynamic_start_weight == prefix_sum_all_weights[rank])
			
			
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
		
		
		
		
		if (rank >= 0)
		{
			std::vector<node_info_for_node_range> node_info_for_node_ranges(local_node_ranges.size());
			dynamic_start_weight = prefix_sum_all_weights[rank];
			std::uint64_t dynamic_start_node = 0; //this is local, not global index!
			std::uint64_t dynamic_node_degree = 0;
			for (int i = 0; i < local_node_ranges.size(); i++)
			{	
				dynamic_start_weight = local_node_ranges[i].start_index;
				node_info_for_node_range info;
				info.start_node = dynamic_start_node + node_offset;
				info.start_node_start_degree = dynamic_node_degree;
				
				if (dynamic_node_degree > 0) //iff last node got splittet
				{
					dynamic_start_weight += indegrees[dynamic_start_node++] - dynamic_node_degree;
				}
				
				while (dynamic_start_weight < local_node_ranges[i].end_index && dynamic_start_node < num_local_vertices)
				{
					dynamic_start_weight += indegrees[dynamic_start_node++] + 1;
				}
				
				
				if (dynamic_start_node == num_local_vertices)
				{
					dynamic_node_degree = 0;
					info.end_node = dynamic_start_node + node_offset - 1;
					info.end_node_end_degree = indegrees[dynamic_start_node-1];
				}
				else if (dynamic_start_weight == local_node_ranges[i].end_index)
				{
					dynamic_node_degree = 0;
					info.end_node = dynamic_start_node + node_offset - 1;
					info.end_node_end_degree = indegrees[dynamic_start_node-1];
				}
				else
				{
					info.end_node_end_degree = dynamic_start_weight - local_node_ranges[i].end_index;
					dynamic_start_node--;
					info.end_node = dynamic_start_node + node_offset;
					dynamic_node_degree = info.end_node_end_degree +1;
				}
				
				
				
				
				std::cout << rank << " with [(start_node, start_degree),(end_node, end_node_end_degree)and node :[(" << info.start_node
				<< "," << info.start_node_start_degree << "),(" << info.end_node << "," << info.end_node_end_degree << ")]" << std::endl;
			}
		}*/
		

	
		struct nodes_to_PE_assignment {
			std::uint32_t targetPE;
			std::int64_t start_part_node;
			std::uint64_t start_part_node_start_degree;
			std::uint64_t node_start_index;
			std::uint64_t node_end_index; //this index is inclusive. So node_end_index = 10 means every node <= 10 is meant
			std::int64_t end_part_node;
			std::uint64_t end_part_node_end_degree;
		};
		
	
		
		std::vector<nodes_to_PE_assignment> nodes_to_PE_assignments(1);
		nodes_to_PE_assignments[0].start_part_node = -1;
		nodes_to_PE_assignments[0].node_start_index = node_offset;
		
		struct cut_node {
			std::uint64_t local_index;
			std::uint64_t cut_degree;
		};
		std::vector<cut_node> cut_nodes(0);
		std::uint64_t dynamic_start_weight = prefix_sum_all_weights[rank];
		std::uint64_t next_weight_cut = ((dynamic_start_weight + 2*num_local_vertices) / (2*num_local_vertices)) * (2*num_local_vertices) - 2*num_local_vertices;
		while (dynamic_start_weight < prefix_sum_all_weights[rank+1])
		{			
			next_weight_cut+= 2*num_local_vertices;
			
			nodes_to_PE_assignment& nodes_to_PE_assignment = nodes_to_PE_assignments.back();
			nodes_to_PE_assignment.targetPE = (next_weight_cut - 1) / (2*num_local_vertices);
			
			if (next_weight_cut >= prefix_sum_all_weights[rank+1])
			{
				nodes_to_PE_assignment.node_end_index = node_offset + num_local_vertices;
				nodes_to_PE_assignment.end_part_node = -1;
				
				//iff there is no weight cut
				
				break;
			}
			else
			{
				std::uint64_t start_node = nodes_to_PE_assignments.back().node_start_index - node_offset;
				while (dynamic_start_weight < next_weight_cut)
					dynamic_start_weight += indegrees[start_node++] +1;
				
				if (dynamic_start_weight == next_weight_cut)
				{
					nodes_to_PE_assignment.node_end_index = start_node + node_offset;
					nodes_to_PE_assignment.end_part_node = -1;
					nodes_to_PE_assignments.resize(nodes_to_PE_assignments.size() + 1);
					nodes_to_PE_assignments.back().node_start_index = start_node + node_offset;
					nodes_to_PE_assignments.back().start_part_node = -1;
					
					
					//iff there is a clean weight cut between following nodes
					

				}
				else
				{
					nodes_to_PE_assignment.node_end_index = start_node - 1 + node_offset;
					nodes_to_PE_assignment.end_part_node = start_node - 1 + node_offset;
					nodes_to_PE_assignment.end_part_node_end_degree = dynamic_start_weight - next_weight_cut + 1;
					nodes_to_PE_assignments.resize(nodes_to_PE_assignments.size() + 1);
					nodes_to_PE_assignments.back().start_part_node = start_node - 1 + node_offset;
					nodes_to_PE_assignments.back().start_part_node_start_degree = dynamic_start_weight - next_weight_cut + 1;
					nodes_to_PE_assignments.back().node_start_index = start_node + node_offset;
					
					
					//iff there is is a cut between edges of the same node

					cut_nodes.push_back({start_node - 1,  dynamic_start_weight - next_weight_cut + 1});
					
					start_node--;
					dynamic_start_weight -= indegrees[start_node] +1;
				}
				
			}
			
			
			
		}
		
		
		/*
		for (int i = 0; i < local_node_info_for_node_range.size(); i++)
		{
			node_info_for_node_range info = local_node_info_for_node_range[i];
			std::cout << rank << " with [(start_node, start_degree),(end_node, end_node_end_degree)and node :[(" << info.start_node
				<< "," << info.start_node_start_degree << "),(" << info.end_node << "," << info.end_node_end_degree << ")]" << std::endl;
		}*/
		
		
		for (int i = 0; i < nodes_to_PE_assignments.size(); i++)
		{
			nodes_to_PE_assignment info = nodes_to_PE_assignments[i];
			std::cout << rank << " hat info für PE " << info.targetPE << ": (" << info.start_part_node << "," << info.start_part_node_start_degree << "),(" <<info.node_start_index << "," << info.node_end_index << "),(" << info.end_part_node << "," << info.end_part_node_end_degree << ")" << std::endl;
		}
		
		
		for (int i = 0; i < cut_nodes.size(); i++)
		{
			cut_node node = cut_nodes[i];

		/* //this is just so i know how the variables are called....
		std::vector<packet> recv_packets = recv.extract_recv_buffer();
		std::vector<std::int32_t> recv_counts = recv.extract_recv_counts();*/


			if (node.cut_degree > 0)
			{
				std::uint64_t local_index = node.local_index;
				std::uint64_t current_degree = 0;
				for (std::uint32_t p = 0; p < size; p++)
				{
					std::uint64_t current_degree_from_PE_with_lower_index = current_degree;
					for (std::uint64_t i = 0; i < recv_counts[p]; i++)
					{
						if (recv_packets[recv_displs[p]+i].node == local_index + node_offset)
						{
							current_degree += recv_packets[recv_displs[p]+i].indegree;
							if (current_degree > node.cut_degree)
							{
								std::cout << "cut node (" << node.local_index + node_offset << "," << node.cut_degree << ") has to be cut at PE " << p << " at " << current_degree - current_degree_from_PE_with_lower_index << "th index" <<  std::endl;
								i = recv_counts[p];
								p=size;
							}
						}
					}
				}
				
			}
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