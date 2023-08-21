#include "kamping/checking_casts.hpp"
#include "kamping/collectives/allgather.hpp"
#include "kamping/communicator.hpp"
#include "kamping/data_buffer.hpp"
#include "kamping/environment.hpp"
#include "kamping/named_parameters.hpp"

#include "timer.cpp"
#include "regular_pointer_doubling.cpp"
struct packet {
	std::int32_t ruler_source;
	std::int32_t destination;
};



/*
here every PE must have the same number of nodes aka the length of successors is the same
also dist_rulers >= 3
*/
class regular_ruling_set
{
	public:
	
	regular_ruling_set(std::vector<std::int32_t>& successors, int32_t dist_rulers)
	{
		s = successors;
		num_local_vertices = s.size();
		num_local_rulers = num_local_vertices / dist_rulers;
		distance_rulers = dist_rulers;
	}
	
	
	void start(kamping::Communicator<>& comm)
	{
		timer timer("algorithmus");
		
		
		
		
		
		
		size = comm.size();
		rank = comm.rank();
		num_global_vertices = num_local_vertices * size;
		node_offset = num_local_vertices * rank;
		
		/*
		std::cout << rank << " mit successor array:\n";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << s[i] << " ";
		std::cout <<", rulers sind die ersten " << num_local_rulers << " nodes" << std::endl;
		*/
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		
		for (std::int32_t local_index = 0; local_index < num_local_rulers; local_index++)
		{
			if (!is_final(local_index))
			{
				std::int32_t targetPE = calculate_targetPE(s[local_index]);
				num_packets_per_PE[targetPE]++;
			}
		}
	
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		//now packets are written
		std::vector<packet> out_buffer(send_displacements[size]);
		for (std::int32_t local_index = 0; local_index < num_local_rulers; local_index++)
		{
			if (!is_final(local_index))
			{
				std::int32_t targetPE = calculate_targetPE(s[local_index]);
				std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
				
				out_buffer[packet_index].ruler_source = local_index + node_offset;
				out_buffer[packet_index].destination = s[local_index];
			}
		}

		auto recv = comm.alltoallv(kamping::send_buf(out_buffer), kamping::send_counts(num_packets_per_PE));
	
		std::vector<packet> recv_buffer = recv.extract_recv_buffer();

		std::vector<std::int32_t> mst(num_local_vertices, -1); //previous ruler
		std::vector<std::int32_t> del(num_local_vertices, -1); //dist to previous ruler
		
		std::int32_t num_reached_nodes = 0;
		bool more_nodes_reached = false;
		
		std::int32_t dist_rulers = num_local_vertices / num_local_rulers;
		std::int32_t max_iteration = dist_rulers * std::log(num_global_vertices);
		std::int32_t iteration;
		for (iteration = 1; iteration < max_iteration; iteration++)
		{
			//if (rank ==  0) std::cout << "iteration " << iteration << " mit " << num_local_vertices - num_reached_nodes << std::endl;
			
			
			std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			more_nodes_reached = false;
			
			for (packet& packet: recv_buffer)
			{
				
				if (packet_will_be_forwarded(packet))
				{
					std::int32_t local_index = packet.destination - node_offset;
					std::int32_t target_node = s[local_index];
					std::int32_t targetPE = calculate_targetPE(target_node);
					num_packets_per_PE[targetPE]++;
				}
			}
			calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
	
			out_buffer.resize(send_displacements[size]);
		
			for (packet& packet: recv_buffer)
			{
				std::int32_t local_index = packet.destination - node_offset;
				
				mst[local_index] = packet.ruler_source;
				del[local_index] = iteration;
				
				num_reached_nodes++;
				more_nodes_reached = true;
				
				if (packet_will_be_forwarded(packet))
				{
					std::int32_t target_node = s[local_index];
					std::int32_t targetPE = calculate_targetPE(target_node);
					std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
					
					out_buffer[packet_index].ruler_source = packet.ruler_source;
					out_buffer[packet_index].destination = target_node;
				
				}
			}

			auto recv = comm.alltoallv(kamping::send_buf(out_buffer), kamping::send_counts(num_packets_per_PE));

			recv_buffer = recv.extract_recv_buffer(); //wird der alte recv_buffer eigentlich gefreed?
			
		}
		
		//falls statistische abschätzung nicht gereicht hat, wird schleife weiter ausgeführt
		while (any_PE_has_work(comm, more_nodes_reached))
		{
			//if (rank == 0) std::cout << rank << " extra runde " << iteration - max_iteration + 1 << std::endl;
			
			
			iteration++;
			
			std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			more_nodes_reached = false;
			
			for (packet& packet: recv_buffer)
			{
				
				if (packet_will_be_forwarded(packet))
				{
					std::int32_t local_index = packet.destination - node_offset;
					std::int32_t target_node = s[local_index];
					std::int32_t targetPE = calculate_targetPE(target_node);
					num_packets_per_PE[targetPE]++;
				}
			}
			calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
	
			out_buffer.resize(send_displacements[size]);
		
			for (packet& packet: recv_buffer)
			{
				std::int32_t local_index = packet.destination - node_offset;
				
				mst[local_index] = packet.ruler_source;
				del[local_index] = iteration;
				
				num_reached_nodes++;
				more_nodes_reached = true;
				
				if (packet_will_be_forwarded(packet))
				{
					std::int32_t target_node = s[local_index];
					std::int32_t targetPE = calculate_targetPE(target_node);
					std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
					
					out_buffer[packet_index].ruler_source = packet.ruler_source;
					out_buffer[packet_index].destination = target_node;
				
				}
			}

			auto recv = comm.alltoallv(kamping::send_buf(out_buffer), kamping::send_counts(num_packets_per_PE));

			recv_buffer = recv.extract_recv_buffer(); //wird der alte recv_buffer eigentlich gefreed?
			
		}
		
		//wir müssen noch anfangsknoten zählen und dann die gesamtzahl als rank des final rulers setzten
		
		std::vector<std::int32_t> send_num_not_reached_nodes(1, num_local_vertices - num_reached_nodes);
		std::vector<std::int32_t> recv_num_not_reached_nodes;
		comm.allgather(kamping::send_buf(send_num_not_reached_nodes), kamping::recv_buf(recv_num_not_reached_nodes));
		std::int32_t sum = 0;
		for (std::int32_t i = 0; i < size; i++)
			sum+= recv_num_not_reached_nodes[i];
		
		std::int32_t local_index_final_node = -1;
		for (std::int32_t local_index = 0; local_index < num_local_rulers; local_index++)
			if (mst[local_index] == -1)
			{
				local_index_final_node = local_index;
				mst[local_index] = local_index + node_offset;
				del[local_index] = 0;
			}
			
			
		std::vector<std::int32_t> s_rec(num_local_rulers);
		std::vector<std::int32_t> r_rec(num_local_rulers);
		for (std::int32_t local_index = 0; local_index < num_local_rulers; local_index++)
		{
			std::int32_t next_ruler = mst[local_index];
			std::int32_t next_ruler_PE = calculate_targetPE(next_ruler);
			s_rec[local_index] = next_ruler - next_ruler_PE * num_local_vertices + next_ruler_PE * num_local_rulers;
			r_rec[local_index] = del[local_index];
		}
		timer.add_checkpoint("recursion");
		
		regular_pointer_doubling algorithm(s_rec, r_rec, local_index_final_node);
		std::vector<std::int32_t> result = algorithm.start(comm);
		
		timer.add_checkpoint("restore values");
		std::vector<std::int32_t> all_results;
		
		comm.allgather(kamping::send_buf(result), kamping::recv_buf(all_results));
		//jetzt müssen werte wiederhergestellt werden
		//dafür müssen alle ruler auf alle PE verteilt werden
		
		result.resize(num_local_vertices);
		for (std::int32_t local_index = num_local_rulers; local_index < num_local_vertices; local_index++)
		{
			if (mst[local_index] != -1)
			{
				std::int32_t prev_ruler = mst[local_index];
				std::int32_t prev_ruler_PE = calculate_targetPE(prev_ruler);
				std::int32_t rank_prev_ruler = all_results[prev_ruler - num_local_vertices * prev_ruler_PE + num_local_rulers * prev_ruler_PE];
				result[local_index] = del[local_index] + rank_prev_ruler;
			}
			else 
			{
				result[local_index] = -1;
			}
		}
		
		timer.finalize(comm, num_global_vertices, distance_rulers);
	/*
		std::cout << rank << " mit result array:\n";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << result[i] << " ";
		std::cout <<std::endl;*/
		//timer.finalize(comm);
		//jetzt alle nodes gereached
		
	}
	
	bool is_global_ruler(std::int32_t global_index)
	{
		std::int32_t targetPE = calculate_targetPE(global_index);
		return is_ruler(global_index - targetPE * num_local_vertices);
	}
	
	bool any_PE_has_work(kamping::Communicator<>& comm, bool this_PE_has_work)
	{
		std::int32_t work = this_PE_has_work;
		std::vector<std::int32_t> send(1,work);
		std::vector<std::int32_t> recv;
		comm.allgather(kamping::send_buf(send), kamping::recv_buf(recv));
		
		for (std::int32_t i = 0; i < size; i++)
			work += recv[i];
		return work > 0;
	}
	
	bool is_final(std::int32_t local_index)
	{
		return local_index + node_offset == s[local_index];
	}
	
	void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t>& num_packets_per_PE)
	{
		send_displacements[0]=0;
		for (std::int32_t i = 1; i < size + 1; i++)
			send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
	}
	
	std::int32_t calculate_targetPE(std::int32_t global_index)
	{
		return global_index / num_local_vertices;
	}
	
	bool is_ruler(std::int32_t local_index)
	{
		return local_index < num_local_rulers;
	}
	
	//a packet will be forwardef iff it doesn't point to ruler and doesn't point to final node
	bool packet_will_be_forwarded(packet packet)
	{
		std::int32_t local_index = packet.destination - node_offset;
		return !is_ruler(local_index) && !is_final(local_index);
	}


	
	private: 
	
	std::vector<std::int32_t> s;
	std::int32_t num_local_vertices;
	std::int32_t num_local_rulers;
	std::int32_t distance_rulers;
	
	std::int32_t size;
	std::int32_t rank;
	std::int32_t num_global_vertices;
	std::int32_t node_offset;
	
	
};