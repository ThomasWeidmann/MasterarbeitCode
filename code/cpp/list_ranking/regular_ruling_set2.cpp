//this class is implemented like in the paper "ultimate parallel list ranking" to minimize number of communication steps in trade of local work

#include "irregular_pointer_doubling.cpp"

class regular_ruling_set2
{
	struct packet{
		std::uint64_t ruler_source;
		std::uint64_t destination;
		std::uint32_t distance;
	};

	

	
	public:
	
	regular_ruling_set2(std::vector<std::uint64_t>& s, std::uint64_t dist_rulers)
	{
		this->s = s;
		this->dist_rulers = dist_rulers;
	}
	
	
	std::vector<std::int64_t> start(kamping::Communicator<>& comm)
	{
		rank = comm.rank();
		size = comm.size();
		num_local_vertices = s.size();
		node_offset = rank * num_local_vertices;
		
		std::vector<std::string> categories = {"local_work", "communication"};
		timer timer("ruler_pakete_senden", categories, "local_work");

		
		//man kann ja wieder die ersten n/dist vielen nodes als ruler setzten. den ruler index speichern. wenn eine packet iteration durch ist, werden erreichte ruler gezählt und genau so viele neue ruler gemacht, in dem rulerindex erhöhrt wird. Dadruch wird nur ein einziges mal extra iteriert
		std::uint64_t out_buffer_size = num_local_vertices/dist_rulers;
		std::vector<packet> out_buffer(out_buffer_size);
		
		/*
		std::cout << rank << " mit successor array:\n";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << s[i] << " ";
		std::cout <<", es werden " << out_buffer_size << " pakete gesendet" << std::endl;*/

		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		std::vector<uint64_t> local_rulers(0);
		std::vector<uint64_t> rulers_to_send_packages(out_buffer_size);
		std::uint64_t ruler_index = 0; //this means that first free rulers has an index >= ruler_index
	
		
		for (std::uint64_t i = 0; i < out_buffer_size; i++)
		{
			while (is_final(ruler_index)) ruler_index++;
			
			local_rulers.push_back(ruler_index);
			rulers_to_send_packages[i] = ruler_index;
			std::int32_t targetPE = calculate_targetPE(s[ruler_index]);
			num_packets_per_PE[targetPE]++;
			ruler_index++;
		}
		
		
		
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		
		
		for (std::uint64_t i = 0; i < out_buffer_size; i++)
		{
			std::uint64_t local_index = rulers_to_send_packages[i];
			
			std::int32_t targetPE = calculate_targetPE(s[local_index]);
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			
			
			out_buffer[packet_index].ruler_source = local_index + node_offset;
			out_buffer[packet_index].destination = s[local_index];
			out_buffer[packet_index].distance = 1;
			
			mark_as_ruler(local_index);
		}

		
		std::vector<std::uint64_t> mst(num_local_vertices);
		std::iota(mst.begin(),mst.end(),node_offset); 
		std::vector<std::int64_t> del(num_local_vertices,0);
		
		timer.add_checkpoint("pakete_verfolgen");

		for (std::uint64_t iteration = 0; iteration < dist_rulers + 3; iteration++)
		{
			/*
			std::cout << rank << " in iteration " << iteration << " with following packages:\n";
			for (packet& packet: out_buffer)
				std::cout << "(" << packet.ruler_source << "," << packet.destination << "," << packet.distance << "),";
			std::cout << std::endl;*/
			
			
			timer.switch_category("communication");
			std::vector<packet> recv_buffer = comm.alltoallv(kamping::send_buf(out_buffer), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
			timer.switch_category("local_work");
			std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			
			
			std::uint64_t num_forwarded_packages = 0;
			for (packet& packet: recv_buffer)
			{
				std::uint64_t local_index = packet.destination - node_offset;
				
				
				mark_as_reached(local_index);
				
				if (!is_final(local_index) && !is_ruler(local_index))
				{
					std::int32_t targetPE = calculate_targetPE(s[local_index]);
					num_packets_per_PE[targetPE]++;
					num_forwarded_packages++;
				}

			}
			
			
			
			//select num_rulers_to_send_packages new rulers if possible
			std::uint64_t num_rulers_to_send_packages = out_buffer_size > num_forwarded_packages ? out_buffer_size - num_forwarded_packages : 0;
			
			rulers_to_send_packages.resize(num_rulers_to_send_packages);
			for (std::uint64_t i = 0; i < num_rulers_to_send_packages; i++)
			{
				while (ruler_index < num_local_vertices && (is_final(ruler_index) || is_reached(ruler_index) || is_ruler(ruler_index))) ruler_index++;
				
				if (ruler_index == num_local_vertices)
				{
					rulers_to_send_packages.resize(i);
					break;
				}
				
				local_rulers.push_back(ruler_index);
				mark_as_ruler(ruler_index);
				rulers_to_send_packages[i] = ruler_index;
				std::int32_t targetPE = calculate_targetPE(s[ruler_index]);
				num_packets_per_PE[targetPE]++;
				
			}
			

			calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
			out_buffer.resize(send_displacements[size]); 
			//now packets are written
			
	 
			for (packet& packet: recv_buffer) //vlt compileranweisung, dass schleife oft ausgeführt wird?
			{
				std::uint64_t local_index = packet.destination - node_offset;
				mst[local_index] = packet.ruler_source;
				del[local_index] = packet.distance;//ich könnte beide arrays zusammenlegen
				
				
				//das könnte auch in eine bit-operation umgewandelt werden
				if (!is_final(local_index) && !is_ruler(local_index)) //if kosten auch viel zeit, aber da führt glaube kein weg vorbei...
				{
					std::int32_t targetPE = calculate_targetPE(s[local_index]); //targetPE mit shift berechnen, division kostet viel zeit
					std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
					
					
					
					out_buffer[packet_index].ruler_source = packet.ruler_source; 
					out_buffer[packet_index].destination = unmask(s[local_index]); //hier kann man den local index schicken, das spart platz
					out_buffer[packet_index].distance = packet.distance + 1;
					
				}	

			}
			timer.switch_category("local_work");

			
			for (std::uint64_t i = 0; i < rulers_to_send_packages.size(); i++)
			{
				std::uint64_t local_index = rulers_to_send_packages[i];
				
				std::int32_t targetPE = calculate_targetPE(s[local_index]);
				std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
				
				out_buffer[packet_index].ruler_source = local_index + node_offset;
				out_buffer[packet_index].destination = unmask(s[local_index]);
				out_buffer[packet_index].distance = 1;
				
				
			}
			
			/*
			//count how many nodes reached
			std::uint64_t num_nodes_reached = 0;
			for (std::uint64_t i = 0; i < num_local_vertices; i++)
				if (is_reached(i))
					num_nodes_reached++;
				
			if (rank == 0)  std::cout << "iteration " << iteration << " with left " << num_local_vertices - num_nodes_reached << std::endl;
			*/
		}
		
		timer.add_checkpoint("rekursion_vorbereiten");

		//now just the global starting node is unreached, this node is also always a ruler
		std::vector<std::uint64_t> num_local_vertices_per_PE;
		timer.switch_category("communication");
		comm.allgather(kamping::send_buf(local_rulers.size()), kamping::recv_buf(num_local_vertices_per_PE));
		timer.switch_category("local_work");
		std::vector<std::uint64_t> prefix_sum_num_vertices_per_PE(size + 1,0);
		for (std::uint32_t i = 1; i < size + 1; i++)
		{
			prefix_sum_num_vertices_per_PE[i] = prefix_sum_num_vertices_per_PE[i-1] + num_local_vertices_per_PE[i-1];
		}
		
		
		std::vector<std::uint64_t> map_ruler_to_its_index(num_local_vertices);
		std::vector<std::uint64_t> s_rec(local_rulers.size());
		std::vector<std::int64_t> r_rec(local_rulers.size());
		std::vector<std::uint32_t> targetPEs_rec(local_rulers.size());
	
		
		std::vector<std::uint64_t> requests(local_rulers.size());
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
		for (std::uint64_t i = 0; i < local_rulers.size(); i++)
		{
			map_ruler_to_its_index[local_rulers[i]] = i;
			std::int32_t targetPE = calculate_targetPE(mst[local_rulers[i]]);
			targetPEs_rec[i] = targetPE;
			num_packets_per_PE[targetPE]++;
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		for (std::uint64_t i = 0; i < local_rulers.size(); i++)
		{
			std::int32_t targetPE = calculate_targetPE(mst[local_rulers[i]]);
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			requests[packet_index] = mst[local_rulers[i]];
		}
		
		timer.switch_category("communication");
		auto recv = comm.alltoallv(kamping::send_buf(requests), kamping::send_counts(num_packets_per_PE));
		timer.switch_category("local_work");

		std::vector<std::uint64_t> recv_requests = recv.extract_recv_buffer();
		
		
		//answers können inplace in requests eingetragen werden
		for (std::uint64_t i = 0; i < recv_requests.size(); i++)
		{
			recv_requests[i] = map_ruler_to_its_index[recv_requests[i]-node_offset] + prefix_sum_num_vertices_per_PE[rank];
		}
		timer.switch_category("communication");
		std::vector<std::uint64_t> recv_answers = comm.alltoallv(kamping::send_buf(recv_requests), kamping::send_counts(recv.extract_recv_counts())).extract_recv_buffer();
		timer.switch_category("local_work");
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
		for (std::uint64_t i = 0; i < local_rulers.size(); i++)
		{
			std::int32_t targetPE = calculate_targetPE(mst[local_rulers[i]]);
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			s_rec[i] = recv_answers[packet_index];
			r_rec[i] = del[local_rulers[i]];
		}
		
		timer.add_checkpoint("rekursion");

		irregular_pointer_doubling algorithm(s_rec, r_rec, targetPEs_rec, prefix_sum_num_vertices_per_PE);
		std::vector<std::int64_t> ranks = algorithm.start(comm);
		timer.add_checkpoint("finalen_ranks_berechnen");

		
		//rank[i + node_offset] = rank[mst[i]] + del[i], and requests[i] = rank[mst[i]] is goal
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
		requests.resize(num_local_vertices);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{
			std::int32_t targetPE = calculate_targetPE(mst[i]);
			num_packets_per_PE[targetPE]++;
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{
			std::int32_t targetPE = calculate_targetPE(mst[i]);
			std::uint64_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			requests[packet_index] = mst[i];
		}
		
		timer.switch_category("communication");
		auto recv2 = comm.alltoallv(kamping::send_buf(requests), kamping::send_counts(num_packets_per_PE));
		timer.switch_category("local_work");
		recv_requests = recv2.extract_recv_buffer();
		num_packets_per_PE = recv2.extract_recv_counts();
		for (std::uint64_t i = 0; i < recv_requests.size(); i++)
		{
			recv_requests[i] = ranks[map_ruler_to_its_index[recv_requests[i] - node_offset]];
		}
		
		timer.switch_category("communication");
		recv_answers = comm.alltoallv(kamping::send_buf(recv_requests), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
		timer.switch_category("local_work");
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
		
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{
			std::int32_t targetPE = calculate_targetPE(mst[i]);
			std::uint64_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			del[i] = size * num_local_vertices - 1 - (del[i] + recv_answers[packet_index]);
		}
		
		timer.finalize(comm, num_local_vertices, dist_rulers);

	
		return del;
	
	}
	
	

	
	
	std::string packet_to_string(packet packet)
	{
		return "(" + std::to_string(packet.ruler_source) + "," + std::to_string(packet.destination) + "," + std::to_string(packet.distance) + ")";
	}
	
	bool packet_will_be_forwarded(packet& packet)
	{
		std::uint64_t local_index = packet.destination - node_offset;
		return !is_final(local_index) && !is_ruler(local_index);
	}
	
	void mark_as_ruler(std::uint64_t local_index)
	{
		s[local_index] =  mark(s[local_index],0);
	}
	
	bool is_ruler(std::uint64_t local_index)
	{
		return is_marked(s[local_index],0);
	}
	
	void mark_as_reached(std::uint64_t local_index)
	{
		s[local_index] =  mark(s[local_index],1);
	}
	
	bool is_reached(std::uint64_t local_index)
	{
		return is_marked(s[local_index],1);
	}
	
	bool is_final(std::uint64_t local_index)
	{
		return local_index + node_offset == unmask(s[local_index]);
	}
	
	std::int32_t calculate_targetPE(std::uint64_t global_index)
	{
		return unmask(global_index) / num_local_vertices;
	}
	
	//sagen wir mal, die obersten 4 bit stehen zum markieren frei
	
	std::uint64_t unmask(std::uint64_t value)
	{
		return value & 0xfffffffffffffff;
	}
	
	// the nth most significant bit will be marked, n>= 0
	std::uint64_t mark(std::uint64_t index, int n)
	{
		return index | (((std::uint64_t) 0x8000000000000000) >> n);
	}
	
	std::uint64_t unmark(std::uint64_t index, int n)
	{
		return index & (0xffffffffffffffff & (~(((std::uint64_t) 0x8000000000000000) >> n)));
	}
	
	bool is_marked(std::uint64_t index, int n)
	{
		return (index & (((std::uint64_t) 0x8000000000000000) >> n)) != 0;
	}
	
 
	void calculate_send_displacements_and_reset_num_packets_per_PE(std::vector<std::int32_t>& send_displacements, std::vector<std::int32_t	>& num_packets_per_PE)
	{
		send_displacements[0]=0;
		for (std::int32_t i = 1; i < size + 1; i++)
			send_displacements[i] = send_displacements[i-1] + num_packets_per_PE[i-1];
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
	}
	
	
	private:
	std::uint64_t node_offset;
	std::uint64_t num_local_vertices;
	std::uint64_t rank, size;
	std::vector<std::uint64_t> s; //s einfach immer übergeben, sonst wird da viel zu viel rumkpiert
	
	std::uint64_t dist_rulers;
};


	