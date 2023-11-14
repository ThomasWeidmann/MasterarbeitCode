#pragma once

template <typename T>
inline auto alltoall(std::vector<std::int32_t>& num_packets_per_PE, std::vector<T>& send, kamping::Communicator<>& comm, karam::mpi::GridCommunicator& grid_comm)
{
		//return grid_alltoall(num_packets_per_PE, send, comm, grid_comm);
		return comm.alltoallv(kamping::send_buf(send), kamping::send_counts(num_packets_per_PE));
}


template <typename T>
auto normal_alltoall(std::vector<std::int32_t>& num_packets_per_PE, std::vector<T>& send, kamping::Communicator<>& comm, karam::mpi::GridCommunicator& grid_comm)
{
	return comm.alltoallv(kamping::send_buf(send), kamping::send_counts(num_packets_per_PE));
}

template <typename T>
auto grid_alltoall(std::vector<std::int32_t>& num_packets_per_PE, std::vector<T>& send, kamping::Communicator<>& comm, karam::mpi::GridCommunicator& grid_comm)
{
	my_grid<T>* grid;
	grid = new my_grid<T>();
	
	grid->grid_send(num_packets_per_PE, send, grid_comm, comm);
	
	return *grid;
}

std::string get_output_string(std::vector<std::uint64_t> to_output)
{
	bool everything = false; //output values from all PEs iff everything otherwise output [min, lower_quartil, median, upper_quartil, max]
	std::string output = "[";
	if (everything)
	{
		for (int i = 0; i < to_output.size(); i++)
			output += std::to_string(to_output[i]) + ",";
		output.pop_back();
		output += "]";	
	}
	else
	{
		std::sort(to_output.begin(), to_output.end());
		int parts = 5;
		for (int i = 0; i < parts; i++)
		{
			int index = (i * (to_output.size() - 1)) / (parts - 1);
			output += std::to_string(to_output[index]) + ",";
		}
		output.pop_back();
		output += "]";
	}
	
	return output;
	
}

void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t	>& num_packets_per_PE, kamping::Communicator<>& comm)
{
	send_displacements[0]=0;
	for (std::int32_t i = 1; i < comm.size() + 1; i++)
		send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
	std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
}