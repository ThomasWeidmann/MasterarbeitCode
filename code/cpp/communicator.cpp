class my_communicator
{
	
	public:
	
	
	template<typename request, typename answer>
	std::vector<answer> request_reply(std::vector<request> requests, std::vector<std::int32_t> send_counts, std::function<answer(const request)> lambda, kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm)
	{
		
		auto recv = comm.alltoallv(kamping::send_buf(requests), kamping::send_counts(send_counts));
		std::vector<request> recv_request = recv.extract_recv_buffer();
		
		int size = recv_request.size();
		std::vector<answer> answers(size);
		
		for (int i = 0; i < size; i++)
			answers[i] = lambda(recv_request[i]);
		
		return comm.alltoallv(kamping::send_buf(answers), kamping::send_counts(recv.extract_recv_counts())).extract_recv_buffer();
	}
	
	template<typename request, typename answer>
	std::vector<answer> request_reply_test(std::vector<request> requests, std::vector<std::int32_t> send_counts, std::function<answer(const request)> lambda, kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm)
	{
		int size = requests.size();
		std::vector<answer> answers(size);
		
		for (int i = 0; i < size; i++)
			answers[i] = lambda(requests[i]);
		
		return answers;
	}
	
	void test(kamping::Communicator<>& comm, karam::mpi::GridCommunicator grid_comm)
	{
		int rank = comm.rank();
		int size = comm.size();
		
		std::vector<std::uint64_t> requests(size, rank);
		std::vector<std::int32_t> send_counts(size,1);
		
		auto lambda = [&](std::uint64_t i ) {return i + rank;}; 
		
		std::vector<std::uint64_t> reply = request_reply<std::uint64_t,std::uint64_t>(requests, send_counts, lambda, comm, grid_comm);
		
		
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
		
		
	}
		
	
};