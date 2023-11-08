class real_load_balance 
{
	
	
	public:
	
	real_load_balance(std::uint64_t comm_rounds)
	{
		this-> comm_rounds = comm_rounds;
	}
	

	void start(std::vector<std::uint64_t>& s, kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm)
	{
		std::vector<std::string> categories = {"calc, other"};
		timer timer("indegrees", categories, "calc", "grid_indegrees_test");
		
		size = comm.size();
		rank = comm.rank();
		num_local_vertices = s.size();
		node_offset = rank * num_local_vertices;
		
		
		std::vector<std::uint64_t> indegrees = calculate_indegrees(s, comm, grid_comm);
		
		for (int i = 0; i < 10; i++)
			std::cout << i + node_offset << " with indegree " << indegrees[i] << std::endl;
		
		calculate_weight_cuts(indegrees, comm);
		
		timer.finalize(comm, "grid_indegrees_test");
	}
	
	void calculate_weight_cuts(std::vector<std::uint64_t>& indegrees, kamping::Communicator<>& comm)
	{
		std::uint64_t weight = num_local_vertices; //because every node has weight 1, the additional 1 for every edge in counted in the loop below
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
			weight += indegrees[i];
		
		std::vector<std::uint64_t> all_weights;
		comm.allgather(kamping::send_buf(weight), kamping::recv_buf<kamping::resize_to_fit>(all_weights));

		std::vector<std::uint64_t> prefix_sum_all_weights(size+1,0);
		for (std::uint32_t i = 1; i < size + 1; i++)
		{
			prefix_sum_all_weights[i] = prefix_sum_all_weights[i-1] + all_weights[i-1];
		}
		
		
		//end_degree means only the first #end_degree
		//start_degree means not the first #start_degree
		struct start_part_node {
			std::uint64_t node;
			std::uint64_t end_degree;
			std::uint32_t targetPE;
		};
		struct middle_part_node {
			std::uint64_t node;
			std::uint64_t start_degree;
			std::uint64_t end_degree;
			std::uint32_t targetPE;
		};
		struct end_part_node {
			std::uint64_t node;
			std::uint64_t start_degree;
			std::uint32_t targetPE;
		};
		
		std::vector<start_part_node> start_part_nodes(0);
		std::vector<middle_part_node> middle_part_nodes(0);
		std::vector<end_part_node> end_part_nodes(0);

		std::uint64_t dynamic_start_weight = prefix_sum_all_weights[rank];
		std::int64_t next_weight_cut = ((dynamic_start_weight + 2*num_local_vertices) / (2*num_local_vertices)) * (2*num_local_vertices) - 2*num_local_vertices;
		std::uint64_t dynamic_start_node = 0;
		while (dynamic_start_weight < prefix_sum_all_weights[rank+1])
		{			
			next_weight_cut+= 2*num_local_vertices;
			
			std::uint32_t targetPE = (next_weight_cut - 1) / (2*num_local_vertices);
		}
	}
	
	std::vector<std::uint64_t> calculate_indegrees(std::vector<std::uint64_t>& s, kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm)
	{
		struct packet {
			std::uint64_t node;
			std::uint64_t indegree;
		};
		
		std::vector<std::int32_t> send_counts(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		std::unordered_map<std::uint64_t, std::uint64_t> local_node_indegrees;

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
			send_counts[targetPE]++;
		}
		
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, send_counts);
		
		std::vector<packet> send_buf(send_displacements[size]);
		for (const auto& [key, value] : local_node_indegrees)
		{
			std::int32_t targetPE = key / num_local_vertices;
			std::int64_t packet_index = send_displacements[targetPE] + send_counts[targetPE]++;
			send_buf[packet_index].node = key;
			send_buf[packet_index].indegree = value;
		}

		std::vector<std::int32_t> send_counts_row = std::vector<std::int32_t>(grid_comm.row_comm().size(),0);
		for (std::int32_t p = 0; p < comm.size(); p++)
		{
			std::int32_t targetPE = grid_comm.proxy_col_index(static_cast<std::size_t>(p));
			send_counts_row[targetPE] += send_counts[p];
		}
		
		auto send_displacements3 = send_counts_row;
		std::exclusive_scan(send_counts_row.begin(), send_counts_row.end(), send_displacements3.begin(), 0ull);
		auto       index_displacements = send_displacements3;
		
		karam::utils::default_init_vector<karam::mpi::IndirectMessage<packet>> contiguous_send_buf(send_buf.size()); 
		
		std::uint64_t index = 0;
		for (std::int32_t p = 0; p < comm.size(); p++)
		{

			for (std::uint64_t i = 0; i < send_counts[p]; i++)
			{
				auto const final_destination = p;
				auto const destination_in_row = grid_comm.proxy_col_index(static_cast<std::size_t>(final_destination));
				auto const idx = index_displacements[destination_in_row]++;

				contiguous_send_buf[static_cast<std::size_t>(idx)] = karam::mpi::IndirectMessage<packet>(
					static_cast<std::uint32_t>(comm.rank()),
					static_cast<std::uint32_t>(final_destination),
					send_buf[index++]
				  );
				
				
			}
			
		}
		
		auto mpi_result_rowwise = grid_comm.row_comm().alltoallv(
		kamping::send_buf(contiguous_send_buf),
		kamping::send_counts(send_counts_row));
		
		auto rowwise_recv_buf    = mpi_result_rowwise.extract_recv_buffer();
		
		local_node_indegrees = std::unordered_map<std::uint64_t, std::uint64_t>();
		std::unordered_map<std::uint64_t, std::uint32_t> node_targetPE;
		
		for (std::uint64_t i = 0; i < rowwise_recv_buf.size(); i++)
		{
			if (local_node_indegrees.contains(rowwise_recv_buf[i].payload().node))
				local_node_indegrees[rowwise_recv_buf[i].payload().node] += rowwise_recv_buf[i].payload().indegree;
			else
				local_node_indegrees[rowwise_recv_buf[i].payload().node] = rowwise_recv_buf[i].payload().indegree;
			node_targetPE[rowwise_recv_buf[i].payload().node] = rowwise_recv_buf[i].get_destination();
		}
		
		std::vector<std::int32_t> num_packets_per_PE2(size,0);
		std::vector<std::int32_t> send_displacements2(size + 1,0);
		
		for (const auto& [key, value] : local_node_indegrees)
		{
			std::int32_t targetPE = node_targetPE[key];
			num_packets_per_PE2[targetPE]++;
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements2, num_packets_per_PE2);
		std::vector<karam::mpi::IndirectMessage<packet>> send(send_displacements2[size]);
		for (const auto& [key, value] : local_node_indegrees)
		{
			std::int32_t targetPE = node_targetPE[key];
			std::int64_t packet_index = send_displacements2[targetPE] + num_packets_per_PE2[targetPE]++;
		
			send[packet_index] = karam::mpi::IndirectMessage<packet>(
					static_cast<std::uint32_t>(comm.rank()),
					static_cast<std::uint32_t>(targetPE),
					{key,value}
				  );
		}
		
		std::vector<karam::mpi::IndirectMessage<packet>> recv = columnwise_exchange(send, grid_comm).extract_recv_buffer();
		//auto rowwise_recv_buf    = mpi_result_rowwise.extract_recv_buffer();
		
		/*
		std::cout << rank << " with: ";
		for (int i = 0; i < recv.size(); i++)
			std::cout << recv[i].payload().node << " hat indegree " << recv[i].payload().indegree << std::endl;
		*/
		local_node_indegrees = std::unordered_map<std::uint64_t, std::uint64_t>();
		for (int i = 0; i < recv.size(); i++)
		{
			if (local_node_indegrees.contains(recv[i].payload().node))
				local_node_indegrees[recv[i].payload().node] += recv[i].payload().indegree;
			else
				local_node_indegrees[recv[i].payload().node] = recv[i].payload().indegree;
		}
		
		std::vector<std::uint64_t> indegrees(num_local_vertices,0);
		
		for (const auto& [key, value] : local_node_indegrees)
		{
			indegrees[key - node_offset] = value;
		}
		return indegrees;
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