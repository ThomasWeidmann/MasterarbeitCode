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