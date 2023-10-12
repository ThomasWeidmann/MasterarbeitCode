template <typename T>
class grid_all_to_all
{
	public:
	
	
	std::vector<T> comm(kamping::Communicator<>& comm, std::vector<std::int32_t> num_packets_per_PE, std::vector<T> send_buffer)
	{
		std::vector<std::int32_t> ranks(2);
		for (int i = 0; i < 2; i++)
			ranks[i] = (comm.rank() % 2) + 2*i;
		
		std::cout << "PE " << comm.rank() << ":[";
		for (int i = 0; i < 2; i++)
			std::cout << ranks[i] << ",";
		std::cout << std::endl;
		
		//kamping::RankRange range = {1,2,3};
		//kamping::RankRanges ranges(range);
		
		kamping::Communicator<> comm2 = comm.create_subcommunicators(ranks);
		
		std::vector<std::int32_t> send_counts = {2,2};
		
		
		return comm2.alltoallv(kamping::send_buf(send_buffer), kamping::send_counts(send_counts)).extract_recv_buffer();
		mpirun -np 4 code d
	}
	
	private:


};