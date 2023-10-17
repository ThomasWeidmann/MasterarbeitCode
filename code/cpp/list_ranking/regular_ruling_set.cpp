#include "kamping/checking_casts.hpp"
#include "kamping/collectives/allgather.hpp"
#include "kamping/communicator.hpp"
#include "kamping/data_buffer.hpp"
#include "kamping/environment.hpp"
#include "kamping/named_parameters.hpp"

#include "../timer.cpp"
#include "regular_pointer_doubling.cpp"
#include "regular_ruling_set_rec.cpp"
#include "../interfaces.cpp"


/*
here every PE must have the same number of nodes aka the length of successors is the same
also dist_rulers >= 3
*/
class regular_ruling_set
{
	struct packet {
		std::int64_t ruler_source;
		std::int64_t destination;
	};

	struct node_packet {
		std::int64_t source;
		std::int64_t destination;
	};
	
	public:
	
	regular_ruling_set(std::vector<std::uint64_t>& successors, int64_t dist_rulers, int64_t iterations)
	{
		s = successors;
		num_local_vertices = s.size();
		num_local_rulers = num_local_vertices / dist_rulers;
		distance_rulers = dist_rulers;
		num_iterations = iterations;
	}
	
	
	std::vector<std::int64_t> start(kamping::Communicator<>& comm, std::vector<std::uint64_t>& successors, karam::mpi::GridCommunicator grid_comm)
	{
		
		
		
		std::vector<std::string> categories = {"local_work", "communication", "other"};
		timer timer("ruler_pakete_senden", categories, "local_work", "regular_ruling_set");
		
		timer.add_info("num_local_vertices", std::to_string(num_local_vertices));
		timer.add_info("dist_rulers", std::to_string(distance_rulers));
		timer.add_info("iterations", std::to_string(num_iterations));
		
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
		

		
		for (std::int64_t local_index = 0; local_index < num_local_rulers; local_index++)
		{
			if (!is_final(local_index))
			{
				std::int64_t targetPE = calculate_targetPE(s[local_index]);
				num_packets_per_PE[targetPE]++;
			}
		}
	
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		//now packets are written
		std::vector<packet> out_buffer(send_displacements[size]);
		for (std::int64_t local_index = 0; local_index < num_local_rulers; local_index++)
		{
			if (!is_final(local_index))
			{
				std::int64_t targetPE = calculate_targetPE(s[local_index]);
				std::int64_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
				
				out_buffer[packet_index].ruler_source = local_index + node_offset;
				out_buffer[packet_index].destination = s[local_index];
			}
		}
		timer.switch_category("communication");
		
		

		auto recv = comm.alltoallv(kamping::send_buf(out_buffer), kamping::send_counts(num_packets_per_PE));
		
		std::vector<packet> recv_buffer = recv.extract_recv_buffer();
		timer.switch_category("local_work");

		std::vector<std::int64_t> mst(num_local_vertices, -1); //previous ruler
		std::vector<std::int64_t> del(num_local_vertices, -1); //dist to previous ruler
		
		std::int64_t num_reached_nodes = 0;
		bool more_nodes_reached = false;
		
		timer.add_checkpoint("pakete_verfolgen");

		std::int64_t max_iteration = distance_rulers * std::log(num_global_vertices);
		std::int64_t iteration=0;
		

		while (iteration++ < max_iteration  || any_PE_has_work(comm, more_nodes_reached))
		{
			//if (rank ==  0) std::cout << "iteration " << iteration << " mit " << num_local_vertices - num_reached_nodes << std::endl;
			
			//timer.add_checkpoint("iteration " + std::to_string(iteration));

			
			std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			more_nodes_reached = false;
			
			for (packet& packet: recv_buffer)
			{
				//std::cout << "iteration " << iteration << ", PE " << rank << " recv packet (src,dst): (" << packet.ruler_source << "," << packet.destination << ")" << std::endl; 
				
				
				if (packet_will_be_forwarded(packet))
				{
					std::int64_t local_index = packet.destination - node_offset;
					std::int64_t target_node = s[local_index];
					std::int64_t targetPE = calculate_targetPE(target_node);
					num_packets_per_PE[targetPE]++;
				}
			}
			calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
	
			out_buffer.resize(send_displacements[size]);
		
			for (packet& packet: recv_buffer)
			{
				std::int64_t local_index = packet.destination - node_offset;
				
				mst[local_index] = packet.ruler_source;
				del[local_index] = iteration;
				
				num_reached_nodes++;
				more_nodes_reached = true;
				
				if (packet_will_be_forwarded(packet))
				{
					std::int64_t target_node = s[local_index];
					std::int64_t targetPE = calculate_targetPE(target_node);
					std::int64_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
					
					out_buffer[packet_index].ruler_source = packet.ruler_source;
					out_buffer[packet_index].destination = target_node;
				
				}
			}
			timer.switch_category("communication");
		
			auto recv = comm.alltoallv(kamping::send_buf(out_buffer), kamping::send_counts(num_packets_per_PE));
			recv_buffer = recv.extract_recv_buffer(); //wird der alte recv_buffer eigentlich gefreed?
			timer.switch_category("local_work");
			
			

		}
		timer.add_checkpoint("rekursion_vorbereiten");

		
		//wir müssen noch anfangsknoten zählen und dann die gesamtzahl als rank des final rulers setzten
		
		std::vector<std::int64_t> send_num_not_reached_nodes(1, num_local_vertices - num_reached_nodes);
		std::vector<std::int64_t> recv_num_not_reached_nodes;
		timer.switch_category("communication");
		
		comm.allgather(kamping::send_buf(send_num_not_reached_nodes), kamping::recv_buf<kamping::resize_to_fit>(recv_num_not_reached_nodes));
		timer.switch_category("local_work");

		std::int64_t sum = -1; //weil der erste ruler auch not reached ist, aber nicht mitgezählt werden soll. sonst nur nicht ruler unreached
		for (std::int64_t i = 0; i < size; i++)
			sum+= recv_num_not_reached_nodes[i];
		
		
		
		std::int64_t local_index_final_node = -1; //das hier ist der erste ruler, aber  im rekursiven aufruf der letzte node insgesamt
		for (std::int64_t local_index = 0; local_index < num_local_rulers; local_index++)
			if (mst[local_index] == -1)
			{
				
				local_index_final_node = local_index;
				mst[local_index] = local_index + node_offset;
				del[local_index] = sum;
			}
			
			
		std::vector<std::uint64_t> s_rec(num_local_rulers);
		std::vector<std::int64_t> r_rec(num_local_rulers);
		for (std::int64_t local_index = 0; local_index < num_local_rulers; local_index++)
		{
			std::int64_t next_ruler = mst[local_index];
			std::int64_t next_ruler_PE = calculate_targetPE(next_ruler);
			s_rec[local_index] = next_ruler - next_ruler_PE * num_local_vertices + next_ruler_PE * num_local_rulers;
			r_rec[local_index] = del[local_index];
		}
	
		timer.add_checkpoint("rekursion");
		timer.switch_category("other");
		std::vector<std::int64_t> result;
		if (num_iterations == 1)
		{
			regular_pointer_doubling algorithm(s_rec, r_rec, local_index_final_node);
			result = algorithm.start(comm, grid_comm);
		}
		else
		{
			regular_ruling_set_rec algorithm(s_rec, r_rec, local_index_final_node, distance_rulers);
			result = algorithm.start(comm, grid_comm);
		}
		timer.add_checkpoint("finalen_ranks_berechnen");
		timer.switch_category("local_work");
		
		for (std::int64_t local_index = 0; local_index < num_local_rulers; local_index++)
		{
			result[local_index] = num_global_vertices - 1 - result[local_index];
		}
		
	
		std::vector<std::int64_t> all_results;
		//falls dist_rulers << p, dann sollten results die benötigt werden requested werden und dann in eine lokale hasmap für effizienten zugriff geschrieben werden
		timer.switch_category("communication");
		
		comm.allgather(kamping::send_buf(result), kamping::recv_buf<kamping::resize_to_fit>(all_results));
		timer.switch_category("local_work");

		//jetzt müssen werte wiederhergestellt werden
		//dafür müssen alle ruler auf alle PE verteilt werden
		
		result.resize(num_local_vertices);
		
	
		std::vector<node_packet> local_unreached_nodes(num_local_vertices - num_reached_nodes + 1);
		std::int64_t local_unreached_nodes_index = 0;
		
		for (std::int64_t local_index = num_local_rulers; local_index < num_local_vertices; local_index++)
		{
			if (mst[local_index] != -1)
			{
				std::int64_t prev_ruler = mst[local_index];
				std::int64_t prev_ruler_PE = calculate_targetPE(prev_ruler);
				std::int64_t rank_prev_ruler = all_results[prev_ruler - num_local_vertices * prev_ruler_PE + num_local_rulers * prev_ruler_PE];
				result[local_index] =  rank_prev_ruler - del[local_index];
			}
			else 
			{
				local_unreached_nodes[local_unreached_nodes_index].source = local_index + node_offset;
				local_unreached_nodes[local_unreached_nodes_index].destination = s[local_index];
				local_unreached_nodes_index++;
			}
		}
		
		local_unreached_nodes.resize(local_unreached_nodes_index);
		std::vector<node_packet> global_unreached_nodes; //das hier sind jetzt genau die nodes, die vor dem ersten ruler sind
		timer.switch_category("communication");

		comm.allgatherv(kamping::send_buf(local_unreached_nodes), kamping::recv_buf<kamping::resize_to_fit>(global_unreached_nodes));
		timer.switch_category("local_work");

		std::unordered_map<std::int64_t, std::int64_t> node_map; //node_map[source] = destination, für jeden unreached node (source,destination)
		std::unordered_map<std::int64_t, std::int64_t> has_pred_map; //has_pred_map[source] = true, if any node source has any pred
		std::int64_t start_node;
		
		for (std::int64_t i = 0; i < global_unreached_nodes.size(); i++)
		{
			node_packet node = global_unreached_nodes[i];
			node_map[node.source] = node.destination;
			has_pred_map[node.destination] = true;
		}
		
		for (std::int64_t i = 0; i < global_unreached_nodes.size(); i++)
		{
			if (!has_pred_map.contains(global_unreached_nodes[i].source))
			{
				start_node = global_unreached_nodes[i].source;
				break;
			}
		}
		std::int64_t node = start_node;
		std::int64_t node_rank = num_global_vertices - 1;
		while (node_map.contains(node))
		{
			if (calculate_targetPE(node) == rank)
				result[node - node_offset] = node_rank;
			node_rank--;
			node = node_map[node];
		}
		
		
		std::string save_dir = "regular_ruling_set";
		if (num_iterations == 2)
			save_dir = "regular_ruling_set_rec";
		timer.finalize(comm, save_dir);

	
	/*
		std::cout << rank << " mit result array:\n";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << result[i] << " ";
		std::cout <<std::endl;
*/
		return result;
		
	}
	
	bool is_global_ruler(std::int64_t global_index)
	{
		std::int64_t targetPE = calculate_targetPE(global_index);
		return is_ruler(global_index - targetPE * num_local_vertices);
	}
	
	//ich habe die vermutung dass allgather kaputt ist...
	bool any_PE_has_work(kamping::Communicator<>& comm, bool this_PE_has_work)
	{
		std::int32_t work = this_PE_has_work;
		std::vector<std::int32_t> send(1,work);
		std::vector<std::int32_t> recv;
		comm.allgather(kamping::send_buf(send), kamping::recv_buf<kamping::resize_to_fit>(recv));
		
		for (std::int32_t i = 0; i < size; i++)
			work += recv[i];
		return work > 0;
	}

	
	bool is_final(std::int64_t local_index)
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
	
	std::int64_t calculate_targetPE(std::int64_t global_index)
	{
		return global_index / num_local_vertices;
	}
	
	bool is_ruler(std::int64_t local_index)
	{
		return local_index < num_local_rulers;
	}
	
	//a packet will be forwardef iff it doesn't point to ruler and doesn't point to final node
	bool packet_will_be_forwarded(packet packet)
	{
		std::int64_t local_index = packet.destination - node_offset;
		return !is_ruler(local_index) && !is_final(local_index);
	}


	
	private: 
	
	std::vector<std::uint64_t> s;
	std::int64_t num_local_vertices;
	std::int64_t num_local_rulers;
	std::int64_t distance_rulers;
	std::int64_t num_iterations;
	
	std::int64_t size;
	std::int64_t rank;
	std::int64_t num_global_vertices;
	std::int64_t node_offset;
	
	
};