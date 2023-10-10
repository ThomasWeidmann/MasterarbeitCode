class mpi_init
{
	public:

	mpi_init(kamping::Communicator<>& comm)
	{
		std::vector<std::int32_t> send_buf(comm.size(),comm.rank());
		std::vector<std::int32_t> send_counts(comm.size(),1);
		
		recv = comm.alltoallv(kamping::send_buf(send_buf), kamping::send_counts(send_counts)).extract_recv_buffer();
	}
	
	
	private:
	
	std::vector<std::int32_t> recv;
};