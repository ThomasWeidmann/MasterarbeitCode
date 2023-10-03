
class tree_irregular_pointer_doubling
{
	struct node_request {
		std::uint64_t node;
		std::uint64_t mst;
	};

	struct answer {
		std::uint64_t node;
		std::uint64_t r_of_mst;
		std::uint64_t mst_of_mst;
		std::uint32_t targetPE_of_mst;
		std::int64_t passive_of_mst;
	};
	

	
	public:
	
	tree_irregular_pointer_doubling(std::vector<std::uint64_t>& s, std::vector<std::uint64_t>& r, std::vector<std::uint32_t>& targetPEs, std::vector<std::uint64_t>& prefix_sum_num_vertices_per_PE, kamping::Communicator<>& comm)
	{
		rank = comm.rank();
		size = comm.size();
		node_offset = prefix_sum_num_vertices_per_PE[rank];
		
		this->s = s;
		this->r = r;
		this->targetPEs = targetPEs;
		this->prefix_sum_num_vertices_per_PE = prefix_sum_num_vertices_per_PE;
		
		
		
		start(comm);
	}
	
	/*
	//all edges will be turn around, therefore we have indegree 1 and outdegree can be any integer. Additionally adj_arr will have no self edges
	void calculate_adj_arr(std::vector<std::uint64_t>& s, std::vector<std::uint64_t> r, std::vector<std::uint32_t> targetPEs, std::vector<std::uint64_t> prefix_sum_num_vertices_per_PE, kamping::Communicator<>& comm)
	{		
		struct edge{
			std::uint64_t source;
			std::uint64_t destination;
			std::uint64_t weight; //aka the number of edges it represents in the original graph
			std::uint32_t source_PE;
		};
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		for (std::uint64_t i = 0; i < s.size(); i++)
		{
			if (i + node_offset == s[i])
				continue;
			
			std::int32_t targetPE = targetPEs[i];
			num_packets_per_PE[targetPE]++;		
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		std::vector<edge> edges(s.size());
		
		for (std::uint64_t i = 0; i < s.size(); i++)
		{
			if (i + node_offset == s[i])
				continue;
			
			std::int32_t targetPE = targetPEs[i];
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			
			edges[packet_index].source = i + node_offset;
			edges[packet_index].destination = s[i];
			edges[packet_index].weight = r[i];
			edges[packet_index].source_PE = rank;
		}
		
		auto recv = comm.alltoallv(kamping::send_buf(edges), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
	
		std::vector<std::uint64_t> turned_target_PEs(recv.size());
		std::vector<std::uint64_t> turned_edges(recv.size());
		std::vector<std::uint64_t> weights(recv.size());
		std::vector<std::uint64_t> edges_per_node(s.size(),0);
		
		for (std::uint64_t i = 0; i < recv.size(); i++)
		{
			edges_per_node[recv[i].destination - node_offset]++;
		}
		std::vector<std::uint64_t> bounds(s.size() + 1, 0);
		for (std::uint64_t i = 1; i <= s.size(); i++)
			bounds[i] = bounds[i-1] + edges_per_node[i-1];
		std::fill(edges_per_node.begin(), edges_per_node.end(), 0);
		for (std::uint64_t i = 0; i < recv.size(); i++)
		{
			std::uint64_t target_node = recv[i].destination - node_offset;
			std::uint64_t packet_index = bounds[target_node] + edges_per_node[target_node]++;
			turned_edges[packet_index] = recv[i].source;
			weights[packet_index] = recv[i].weight;
			turned_target_PEs[packet_index] = recv[i].source_PE;
		}
	}*/
	
	
	void start(kamping::Communicator<>& comm)
	{
		rank = comm.rank();
		size = comm.size();
		num_local_vertices = s.size();
		num_global_vertices = prefix_sum_num_vertices_per_PE[size];
		node_offset = prefix_sum_num_vertices_per_PE[rank];
		
		q = s;
		std::vector<bool> passive(num_local_vertices, false);
		
		std::uint64_t active_nodes = num_local_vertices;
		for (std::int32_t local_index = 0; local_index < num_local_vertices; local_index++)
			if (q[local_index] == local_index + node_offset)
			{
				passive[local_index] = true;
				active_nodes--;
			}
		
		
	
		
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		std::vector<node_request> requests(num_local_vertices);
		std::vector<answer> answers(num_local_vertices);

		
		
	
		
		while (any_PE_has_work(comm, active_nodes > 0))
		{
			
			

			std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			for (std::int32_t local_index = 0; local_index < num_local_vertices;local_index++)
			{
				if (!passive[local_index])
				{
					std::int32_t targetPE = targetPEs[local_index];
					num_packets_per_PE[targetPE]++;
				}
			}
			calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
			requests.resize(send_displacements[size]);
			
			//dann requests gefülllt
			for (std::int32_t local_index = 0; local_index < num_local_vertices;local_index++)
			{
				if (!passive[local_index])
				{
					std::int32_t targetPE = targetPEs[local_index];
					std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
					requests[packet_index].node = local_index + node_offset;
					requests[packet_index].mst = q[local_index];
					
				}
				
			}
			
			auto recv = comm.alltoallv(kamping::send_buf(requests), kamping::send_counts(num_packets_per_PE));
			std::vector<node_request> recv_requests = recv.extract_recv_buffer();
			num_packets_per_PE = recv.extract_recv_counts();
			answers.resize(recv_requests.size());
			
		
			
			for (std::int32_t i = 0; i < recv_requests.size(); i++)
			{
				

				
				std::int32_t local_index = recv_requests[i].mst - node_offset;
				answers[i].node = recv_requests[i].node;
				answers[i].r_of_mst = r[local_index];
				answers[i].mst_of_mst = q[local_index];
				answers[i].targetPE_of_mst = targetPEs[local_index];
				answers[i].passive_of_mst = passive[local_index];
			}
		
			std::vector<answer> recv_answers = comm.alltoallv(kamping::send_buf(answers), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
			//dann answers eingetragen
			
			
			for (std::int32_t i = 0; i < recv_answers.size(); i++)
			{				
				std::int32_t local_index = recv_answers[i].node - node_offset;
				
				targetPEs[local_index] = recv_answers[i].targetPE_of_mst;
				q[local_index] = recv_answers[i].mst_of_mst;
				r[local_index] = r[local_index] + recv_answers[i].r_of_mst;
				passive[local_index] = recv_answers[i].passive_of_mst;
				active_nodes -= passive[local_index];
			}
			
		}
		
		/*
		for (int i = 0; i< num_local_vertices; i++)
			std::cout << i + node_offset << "," << q[i] << "," << r[i] << std::endl;
		*/
		
		
	}
	
	
		
	
	void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t>& num_packets_per_PE)
	{
		send_displacements[0]=0;
		for (std::int32_t i = 1; i < size + 1; i++)
			send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
	}
	
	bool any_PE_has_work(kamping::Communicator<>& comm, bool this_PE_has_work)
	{
		std::int32_t work = this_PE_has_work;
		std::vector<std::int32_t> send(1,work);
		std::vector<std::int32_t> recv;
		comm.allgather(kamping::send_buf(send), kamping::recv_buf(recv));
		
		for (std::int32_t i = 0; i < size; i++)
			work += recv[i];
		return work > 0;
	}
	
	public:
	std::uint64_t node_offset;
	std::uint64_t num_local_vertices;
	std::uint64_t num_global_vertices;
	std::uint64_t rank, size;
	std::vector<std::uint64_t> s;
	std::vector<std::uint64_t> q;
	std::vector<std::uint64_t> r;
	std::vector<std::uint32_t> targetPEs;
	std::vector<std::uint64_t> prefix_sum_num_vertices_per_PE;
};