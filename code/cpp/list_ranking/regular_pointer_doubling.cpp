#pragma once



class regular_pointer_doubling
{
	struct node_request {
		std::uint64_t node;
		std::uint64_t mst;
	};

	struct answer {
		std::uint64_t node;
		std::int64_t r_of_mst;
		std::uint64_t mst_of_mst;
		bool passive_of_mst;
	};
	
	public:
	
	//if this PE has final node, then final node is set to a valid value, otherweise it is -1
	regular_pointer_doubling(std::vector<std::uint64_t>& successors, std::vector<std::int64_t>& ranks, std::uint64_t local_index_final_node)
	{
		s = successors;
		num_local_vertices = s.size();
		q = s;
		r = ranks;

		passive = std::vector<bool>(num_local_vertices, false);
		
		if (local_index_final_node != -1)
			passive[local_index_final_node] = true;
	}
	
	
	regular_pointer_doubling(std::vector<std::uint64_t>& successors, kamping::Communicator<>& comm)
	{
		s = successors;
		num_local_vertices = s.size();
		q = s;
		node_offset = num_local_vertices * comm.rank();
		
		r = std::vector<std::int64_t>(num_local_vertices, 1);
		passive = std::vector<bool>(num_local_vertices, false);
		
		for (std::int32_t local_index = 0; local_index < num_local_vertices; local_index++)
			if (q[local_index] == local_index + node_offset)
			{
				r[local_index] = 0;
				passive[local_index] = true;
			}
	}
	
	
	std::vector<std::int64_t> start(kamping::Communicator<>& comm)
	{
		std::vector<std::string> categories = {"local_work", "communication"};
		timer timer("start", categories, "local_work", "regular_pointer_doubling");
		
		timer.add_info(std::string("num_local_vertices"), std::to_string(num_local_vertices));

		
		size = comm.size();
		rank = comm.rank();
		
		/*
		std::cout << rank << " mit successor array:\n";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << s[i] << " ";
		std::cout <<std::endl;
		*/
		num_global_vertices = num_local_vertices * size;
		node_offset = num_local_vertices * rank;
		
		std::vector<node_request> requests(num_local_vertices);
		std::vector<answer> answers(num_local_vertices);

		
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		std::int32_t max_iteration = std::log2(num_global_vertices) + 2;
		for (std::int32_t iteration = 0; iteration < max_iteration; iteration++)
		{
			
			timer.add_checkpoint("iteration " + std::to_string(iteration));
			
			//zuerst request packets gezählt
			std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			for (std::int32_t local_index = 0; local_index < num_local_vertices;local_index++)
			{
				if (!passive[local_index])
				{
					std::int32_t targetPE = calculate_targetPE(q[local_index]);
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
					std::int32_t targetPE = calculate_targetPE(q[local_index]);
					std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
					requests[packet_index].node = local_index + node_offset;
					requests[packet_index].mst = q[local_index];
				}
				
			}
			timer.switch_category("communication");
		
			std::vector<node_request> recv_requests = comm.alltoallv(kamping::send_buf(requests), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
			timer.switch_category("local_work");

			//dann answers gezählt
			
			answers.resize(recv_requests.size());
			std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			for (std::int32_t i = 0; i < recv_requests.size(); i++)
			{
				std::int32_t targetPE = calculate_targetPE(recv_requests[i].node);
				num_packets_per_PE[targetPE]++;
			}
			//dann answers ausgefüllt
			calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
			for (std::int32_t i = 0; i < recv_requests.size(); i++)
			{
				std::int32_t targetPE = calculate_targetPE(recv_requests[i].node);
				std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
				
				std::int32_t local_index = recv_requests[i].mst - node_offset;
				answers[packet_index].node = recv_requests[i].node;
				answers[packet_index].r_of_mst = r[local_index];
				answers[packet_index].mst_of_mst = q[local_index];
				answers[packet_index].passive_of_mst = passive[local_index];
			}
			timer.switch_category("communication");
	
			std::vector<answer> recv_answers = comm.alltoallv(kamping::send_buf(answers), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
			timer.switch_category("local_work");
				//dann answers eingetragen
			
			for (std::int32_t i = 0; i < recv_answers.size(); i++)
			{
				std::int32_t local_index = recv_answers[i].node - node_offset;
				
				q[local_index] = recv_answers[i].mst_of_mst;
				r[local_index] = r[local_index] + recv_answers[i].r_of_mst;
				passive[local_index] = recv_answers[i].passive_of_mst;
			}
	
		}
		//timer.finalize(comm, "regular_pointer_doubling");

		/*
		std::cout << rank << " mit rank array:\n";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << r[i] << " ";
		std::cout <<std::endl;
		*/
		return r;
		
	}
	
	void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t>& num_packets_per_PE)
	{
		send_displacements[0]=0;
		for (std::int32_t i = 1; i < size + 1; i++)
			send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
	}
	
	std::int32_t calculate_targetPE(std::uint64_t global_index)
	{
		return global_index / num_local_vertices;
	}
	
	
	private:
	std::uint64_t num_local_vertices;
	std::uint64_t final_node;
	std::vector<std::uint64_t> s;
	std::vector<std::uint64_t> q;
	std::vector<std::int64_t> r;
	std::vector<bool> passive;
	
	std::int32_t rank;
	std::int32_t size;
	std::int32_t num_global_vertices;
	std::int32_t node_offset;

};