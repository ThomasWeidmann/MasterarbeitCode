#include "communicator.cpp"

#pragma once

void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t	>& num_packets_per_PE, kamping::Communicator<>& comm)
{
	send_displacements[0]=0;
	for (std::int32_t i = 1; i < comm.size() + 1; i++)
		send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
	std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
}


bool is_pow_of_two(std::uint32_t p)
{
	return std::lround(std::pow(2, std::lround(std::log2(p)))) == p;
}

bool any_PE_has_work(kamping::Communicator<>& comm, karam::mpi::GridCommunicator& grid_comm, timer& timer, bool this_PE_has_work, bool grid)
{
	std::vector<std::int32_t> work_vec = allgatherv(timer, (std::int32_t) this_PE_has_work, comm, grid_comm, grid);
	
	std::int32_t size = comm.size();
	std::int32_t work = 0;
	for (std::int32_t i = 0; i < size; i++)
		work += work_vec[i];
	return work > 0;
}