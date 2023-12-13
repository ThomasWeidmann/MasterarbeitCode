#pragma once

void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t	>& num_packets_per_PE, kamping::Communicator<>& comm)
{
	send_displacements[0]=0;
	for (std::int32_t i = 1; i < comm.size() + 1; i++)
		send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
	std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
}