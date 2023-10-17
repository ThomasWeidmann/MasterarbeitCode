#include "karam/mpi/grid_alltoall.hpp"

template <typename T> 
class my_grid
{
	public:
	
	void grid_send(std::vector<std::int32_t>& num_packets_per_PE, std::vector<T>& send, karam::mpi::GridCommunicator& grid_com, kamping::Communicator<>& comm)
	{
		struct send_element
		{
			T value;
			std::int32_t target_PE;
		};
		
		std::vector<send_element> real_send(send.size());
		for (std::uint32_t i = 0; i < send.size(); i++)
		{
			real_send[i].value = send[i];
		}
		
		std::uint32_t index = 0;
		for (std::uint32_t i = 0; i < comm.size(); i++)
			for (std::uint32_t j = 0; j < num_packets_per_PE[i]; j++)
				real_send[index++].target_PE = i;
		
		
		
		auto get_destination = [](const send_element& e) {
			return e.target_PE;
		};

		//std::vector<karam::mpi::IndirectMessage<send_element>> result = grid_mpi_all_to_all(real_send, get_destination, grid_com).extract_recv_buffer();
		std::vector<karam::mpi::IndirectMessage<send_element>> result = grid_mpi_all_to_all(real_send, get_destination, grid_com).extract_recv_buffer();

		
		recv_counts = std::vector<std::int32_t>(comm.size());
		recv_buffer = std::vector<T>(result.size());
		
		for (std::uint32_t i = 0; i < result.size(); i++)
		{
			recv_buffer[i] = result[i].payload().value;
			recv_counts[result[i].get_source()]++;
		}
	
		
	
	}
	
	
	std::vector<std::int32_t>& extract_recv_counts()
	{
		return recv_counts;
	}
	
	std::vector<T>& extract_recv_buffer()
	{
		return recv_buffer;
	}
	
	private:
	
	std::vector<std::int32_t> recv_counts;
	std::vector<T> recv_buffer;
	
	
	
};