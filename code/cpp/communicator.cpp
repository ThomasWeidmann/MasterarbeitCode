//i will introduce three different communication modes
//communication_mode = 0
//this means all communication is direct but request and replies are optimized
//communication_mode = 1
//this menas request and replies are direct and optimized and the rest is grid communication
//communication_mode = 2
//this menas all communication is grid communication
	
	template<typename packet>
	static std::vector<packet> alltoall(timer& timer, std::vector<packet>& send_buf, std::vector<std::int32_t>& send_counts, kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm, int communication_mode)
	{
		switch (communication_mode)
		{
			case 0:
				return alltoall_normal(timer, send_buf, send_counts, comm);
			case 1:
				return alltoall_grid(timer, send_buf, send_counts, comm, grid_comm);
			case 2:
				return alltoall_grid(timer, send_buf, send_counts, comm, grid_comm);
			default:
				return std::vector<packet>();
		}
	}
	
	template<typename request, typename answer>
	static std::vector<answer> request_reply(timer& timer, std::vector<request>& requests, std::vector<std::int32_t>& send_counts, std::function<answer(const request)> lambda, kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm, int communication_mode)
	{
		switch (communication_mode)
		{
			case 0:
				return request_reply_normal<request,answer>(timer, requests, send_counts, lambda, comm);
			case 1:
				return request_reply_normal<request,answer>(timer, requests, send_counts, lambda, comm);
			case 2:
				return request_reply_grid<request,answer>(timer, requests, send_counts, lambda, comm, grid_comm);
			default:
				return std::vector<answer>();
		}
	}


	template<typename packet>
	static std::vector<packet> alltoall_normal(timer& timer, std::vector<packet>& send_buf, std::vector<std::int32_t>& send_counts, kamping::Communicator<>& comm)
	{
		timer.switch_category("communication");
		std::vector<packet> recv =  comm.alltoallv(kamping::send_buf(send_buf), kamping::send_counts(send_counts)).extract_recv_buffer();
		timer.switch_category("local_work");
		return recv;
	}
	
	template<typename packet>
	static std::vector<packet> alltoall_grid(timer& timer, std::vector<packet> send_buf, std::vector<std::int32_t> send_counts, kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm)
	{
		timer.switch_category("communication");
		auto recv = my_grid_all_to_all(send_buf, send_counts, grid_comm, comm).extract_recv_buffer();
		timer.switch_category("local_work");

		std::vector<packet> recv_vector(recv.size());
		for (std::uint64_t i = 0; i < recv.size(); i++)
			recv_vector[i] = recv[i].payload();
		
		return recv_vector;
	}

	
	template<typename request, typename answer>
	static std::vector<answer> request_reply_normal(timer& timer, std::vector<request>& requests, std::vector<std::int32_t>& send_counts, std::function<answer(const request)> lambda, kamping::Communicator<>& comm)
	{
		timer.switch_category("communication");
		auto recv = comm.alltoallv(kamping::send_buf(requests), kamping::send_counts(send_counts));
		timer.switch_category("local_work");
	
		std::vector<request> recv_request = recv.extract_recv_buffer();
		
		std::uint64_t size = recv_request.size();
		std::vector<answer> answers(size);
		
		for (std::uint64_t i = 0; i < size; i++)
		{
			answers[i] = lambda(recv_request[i]);

		}
		
		timer.switch_category("communication");
		std::vector<answer> recv_answers = comm.alltoallv(kamping::send_buf(answers), kamping::send_counts(recv.extract_recv_counts())).extract_recv_buffer();
		timer.switch_category("local_work");

		return recv_answers;
	}
	
	template<typename request, typename answer>
	static std::vector<answer> request_reply_grid(timer& timer, std::vector<request>& requests, std::vector<std::int32_t>& send_counts, std::function<answer(const request)> lambda, kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm)
	{
		timer.switch_category("communication");
		auto recv_request = my_grid_all_to_all(requests, send_counts, grid_comm, comm).extract_recv_buffer();
		timer.switch_category("local_work");

		std::uint64_t size = recv_request.size();
		std::vector<answer> answers(recv_request.size());
		
		std::vector<std::int32_t> recv_counts(comm.size(), 0);
		for (std::uint64_t i = 0; i < size; i++)
		{
			recv_counts[recv_request[i].get_source()]++;
		}
		std::vector<std::uint64_t> recv_displacements(comm.size(), 0);
		for (std::int32_t p = 1; p < comm.size(); p++)
			recv_displacements[p] = recv_displacements[p-1] + recv_counts[p-1];
		std::fill(recv_counts.begin(), recv_counts.end(), 0);

		for (std::uint64_t i = 0; i < size; i++)
		{
			std::int32_t targetPE = recv_request[i].get_source();
			std::uint64_t packet_index = recv_displacements[targetPE] + recv_counts[targetPE]++;
			answers[packet_index] = lambda(recv_request[i].payload());
			
		}
		timer.switch_category("communication");
		auto recv_indirect_answers = my_grid_all_to_all(answers, recv_counts, grid_comm, comm).extract_recv_buffer();
		timer.switch_category("local_work");
		std::vector<std::uint64_t> send_displacements(comm.size(),0);
		for (std::int32_t p = 1; p < comm.size(); p++)
			send_displacements[p] = send_displacements[p-1] + send_counts[p-1];
		
		std::fill(send_counts.begin(), send_counts.end(), 0);

		size = recv_indirect_answers.size();
		std::vector<answer> recv_answers(size);
		for (std::uint64_t i = 0; i < size; i++)
		{	
			std::int32_t sourcePE = recv_indirect_answers[i].get_source();
			std::uint64_t packet_index = send_displacements[sourcePE] + send_counts[sourcePE]++;
			recv_answers[packet_index] = recv_indirect_answers[i].payload();
		}
		return  recv_answers;
	}
	
	
	
	static void test(kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm)
	{
		int rank = comm.rank();
		int size = comm.size();
		
		std::vector<std::uint64_t> requests(size, rank);
		std::vector<std::int32_t> send_counts(size,1);
		
		auto lambda = [&](std::uint64_t i ) {return i + rank;}; 
		/*
		std::vector<std::uint64_t> reply = request_reply_normal<std::uint64_t,std::uint64_t>(requests, send_counts, lambda, comm);
		
		
		std::cout << rank << " with: ";
		for (int i = 0; i < reply.size(); i++)
			std::cout << reply[i] << " ";
		std::cout << std::endl;
		
		return;
		int i = 0;
		
		auto test = [&](float a) {
            return a + i;
        };
		
		std::cout << test(0) << std::endl;
		
		i = 1;
		
		std::cout << test(0) << std::endl;
		*/
		
	}
		
	
