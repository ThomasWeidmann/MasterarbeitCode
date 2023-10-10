#include "wood_irregular_pointer_doubling.cpp"

class wood_regular_ruling_set2 //this is for trees
{
	struct packet{
		std::uint64_t ruler_source;
		std::uint64_t destination;
		std::uint32_t distance;
	};

	struct adj_arr{
		std::vector<std::uint64_t> edges;
		std::vector<std::uint64_t> bounds;
	};
	
	struct result{
			std::uint64_t root;
			std::int64_t distance;
	};
	
	public:
	
	wood_regular_ruling_set2(std::vector<std::uint64_t>& s, std::uint64_t comm_rounds, kamping::Communicator<>& comm)
	{
		std::vector<std::string> categories = {"local_work", "communication"};
		timer timer("graph_umdrehen", categories, "local_work");
		
		size = comm.size();
		rank = comm.rank();
		num_local_vertices = s.size();
		this->comm_rounds = comm_rounds;
		node_offset = rank * num_local_vertices;
		//now turn around s array
		
		adj_arr adj_arr = calculate_adj_arr(s, comm, node_offset);
		edges = adj_arr.edges;
		bounds = adj_arr.bounds;
		
		start(comm, timer);
	}
	
	//all edges will be turn around, therefore we have indegree 1 and outdegree can be any integer. Additionally adj_arr will have no self edges
	adj_arr calculate_adj_arr(std::vector<std::uint64_t>& s, kamping::Communicator<>& comm, std::uint64_t node_offset)
	{
		struct edge{
			std::uint64_t source;
			std::uint64_t destination;
		};
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		for (std::uint64_t i = 0; i < s.size(); i++)
		{
			if (i + node_offset == s[i])
				continue;
			
			std::int32_t targetPE = s[i] / s.size();
			num_packets_per_PE[targetPE]++;		
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		std::vector<edge> edges(s.size());
		
		for (std::uint64_t i = 0; i < s.size(); i++)
		{
			if (i + node_offset == s[i])
				continue;
			
			std::int32_t targetPE = s[i] / s.size();
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			
			edges[packet_index].source = i + node_offset;
			edges[packet_index].destination = s[i];
		}
		
	

		
		auto recv = comm.alltoallv(kamping::send_buf(edges), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
	
		std::vector<std::uint64_t> turned_edges(recv.size());
		std::vector<std::uint64_t> edges_per_node(s.size(),0);
		
		for (std::uint64_t i = 0; i < recv.size(); i++)
		{
			edges_per_node[recv[i].destination - node_offset]++;
		}
		std::vector<std::uint64_t> bounds(s.size() + 1, 0);
		for (std::uint64_t i = 1; i <= s.size(); i++)
			bounds[i] = bounds[i-1] + edges_per_node[i-1];
		std::fill(edges_per_node.begin(), edges_per_node.end(), 0);
		for (std::uint64_t i = 0; i < recv.size(); i++)
		{
			std::uint64_t target_node = recv[i].destination - node_offset;
			std::uint64_t packet_index = bounds[target_node] + edges_per_node[target_node]++;
			turned_edges[packet_index] = recv[i].source;
	
			
		}
		/*
		//print for testing
		std::cout << "s: ";
		for (int i = 0; i < s.size(); i++)
			std::cout << s[i] << ",";
		std::cout << "\nbounds: ";
		for (int i = 0; i < bounds.size(); i++)
			std::cout << bounds[i] << ",";
		std::cout << "\nturned edges: ";
		for (int i = 0; i < turned_edges.size(); i++)
			std::cout << turned_edges[i] << ",";
		std::cout << std::endl;*/
		adj_arr adj_arr;
		adj_arr.edges = turned_edges;
		adj_arr.bounds = bounds;
		return adj_arr;

	}
	
	
	void start(kamping::Communicator<>& comm, timer timer)
	{
		
		
		
		timer.add_checkpoint("ruler_pakete_senden");

		std::uint64_t expected_num_packets = num_local_vertices/comm_rounds;
		std::vector<packet> out_buffer(0);
		
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		std::vector<std::uint64_t> local_rulers(0);
		std::uint64_t ruler_index = 0; //this means that first free rulers has an index >= ruler_index
	
		std::uint64_t num_packets = 0;
		while (num_packets < expected_num_packets && ruler_index < num_local_vertices)
		{
			
			
			
			local_rulers.push_back(ruler_index);
		
			
			for (std::uint64_t i = unmask(bounds[ruler_index]); i < unmask(bounds[ruler_index+1]); i++)
			{
				std::int32_t targetPE = calculate_targetPE(edges[i]);
				num_packets_per_PE[targetPE]++;
				num_packets++;
			}
			ruler_index++;
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		out_buffer.resize(send_displacements[size]);
		for (std::uint64_t i = 0; i < local_rulers.size(); i++)
		{
			std::uint64_t local_index = local_rulers[i];
			
			for (std::uint64_t i = unmask(bounds[local_index]); i < unmask(bounds[local_index+1]); i++)
			{
				std::int32_t targetPE = calculate_targetPE(edges[i]);
				std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
				
				out_buffer[packet_index].ruler_source = local_index + node_offset;
				out_buffer[packet_index].destination = edges[i];
				out_buffer[packet_index].distance = 1;
				
			}
			
		
			
			
			
			
			mark_as_ruler(local_index);
		}

	
		std::vector<std::uint64_t> mst(num_local_vertices);
		std::iota(mst.begin(),mst.end(),node_offset); 
		std::vector<std::uint64_t> del(num_local_vertices,0);
	

		//for (std::uint64_t iteration = 0; iteration <= comm_rounds; iteration++)
		timer.add_checkpoint("pakete_verfolgen");

		bool work_left = true;
		while (any_PE_has_work(comm, work_left))
		{
			/*
			std::cout << rank << " in iteration " << iteration << " with following packages:\n";
			for (packet& packet: out_buffer)
				std::cout << "(" << packet.ruler_source << "," << packet.destination << "," << packet.distance << "),";
			std::cout << std::endl;*/
			
			
			std::vector<packet> recv_buffer = comm.alltoallv(kamping::send_buf(out_buffer), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
			std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
			
			work_left = recv_buffer.size() > 0;
			
			std::uint64_t num_forwarded_packages = 0;
			for (packet& packet: recv_buffer)
			{
				
				std::uint64_t local_index = packet.destination - node_offset;
				
				
				mark_as_reached(local_index);
				if (!is_ruler(local_index))
				{
					for (std::uint64_t i = unmask(bounds[local_index]); i < unmask(bounds[local_index+1]); i++)
					{
						std::int32_t targetPE = calculate_targetPE(edges[i]);
						num_packets_per_PE[targetPE]++;
						num_forwarded_packages++;
						
					}
				}

			}
			
			std::vector<std::uint64_t> rulers_to_send_packages(0);
			while (num_forwarded_packages < expected_num_packets)
			{
				while (ruler_index < num_local_vertices && (is_reached(ruler_index) || is_ruler(ruler_index))) ruler_index++;
				
				if (ruler_index == num_local_vertices)
				{
					break;
				}
				local_rulers.push_back(ruler_index);
				mark_as_ruler(ruler_index);
				rulers_to_send_packages.push_back(ruler_index);
				for (std::uint64_t i = unmask(bounds[ruler_index]); i < unmask(bounds[ruler_index+1]); i++)
				{
					std::int32_t targetPE = calculate_targetPE(edges[i]);
					num_packets_per_PE[targetPE]++;
					num_forwarded_packages++;
		
				}
		
			}
			
			calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
			out_buffer.resize(send_displacements[size]); 
			
			for (packet& packet: recv_buffer) 
			{
				std::uint64_t local_index = packet.destination - node_offset;
				mst[local_index] = packet.ruler_source;
				del[local_index] = packet.distance;
				
				
				
				if (!is_ruler(local_index))
				{
					for (std::uint64_t i = unmask(bounds[local_index]); i < unmask(bounds[local_index+1]); i++)
					{
						std::int32_t targetPE = calculate_targetPE(edges[i]);
						std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;

						out_buffer[packet_index].ruler_source = packet.ruler_source; 
						out_buffer[packet_index].destination = unmask(edges[i]);
						out_buffer[packet_index].distance = packet.distance + 1;
						
						//std::cout << packet_to_string(out_buffer[packet_index]) << " forwarded in iteration " << iteration << "\n";

					}					
				}	

			}
			
			for (std::uint64_t i = 0; i < rulers_to_send_packages.size(); i++)
			{
				std::uint64_t local_index = rulers_to_send_packages[i];
				for (std::uint64_t i = unmask(bounds[local_index]); i < unmask(bounds[local_index+1]); i++)
				{
					
					
					std::int32_t targetPE = calculate_targetPE(edges[i]);
					std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;

					out_buffer[packet_index].ruler_source = local_index + node_offset;
					out_buffer[packet_index].destination = unmask(edges[i]);
					out_buffer[packet_index].distance = 1;
					
					//std::cout << packet_to_string(out_buffer[packet_index]) << " started in iteration " << iteration << "\n";

				}
				
				
			}
			
		}	
		timer.add_checkpoint("rekursion_vorbereiten");

		/*
		for (int i = 0; i<num_local_vertices; i++)
			std::cout << "final: mst[" << i + node_offset << "]=" << mst[i] << ", del[" << i+node_offset << "]=" << del[i] << (is_ruler(i)?" sruler\n":"\n");
		*/
		std::vector<std::uint64_t> num_local_vertices_per_PE;
		comm.allgather(kamping::send_buf(local_rulers.size()), kamping::recv_buf(num_local_vertices_per_PE));
		std::vector<std::uint64_t> prefix_sum_num_vertices_per_PE(size + 1,0);
		for (std::uint32_t i = 1; i < size + 1; i++)
		{
			prefix_sum_num_vertices_per_PE[i] = prefix_sum_num_vertices_per_PE[i-1] + num_local_vertices_per_PE[i-1];
		}
		
		
		std::vector<std::uint64_t> map_ruler_to_its_index(num_local_vertices);
		std::vector<std::uint64_t> s_rec(local_rulers.size());
		std::vector<std::uint64_t> r_rec(local_rulers.size());
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
		
		auto recv = comm.alltoallv(kamping::send_buf(requests), kamping::send_counts(num_packets_per_PE));
		
		std::vector<std::uint64_t> recv_requests = recv.extract_recv_buffer();
		
		
		//answers können inplace in requests eingetragen werden
		for (std::uint64_t i = 0; i < recv_requests.size(); i++)
		{
			recv_requests[i] = map_ruler_to_its_index[recv_requests[i]-node_offset] + prefix_sum_num_vertices_per_PE[rank];
		}
		std::vector<std::uint64_t> recv_answers = comm.alltoallv(kamping::send_buf(recv_requests), kamping::send_counts(recv.extract_recv_counts())).extract_recv_buffer();
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
		for (std::uint64_t i = 0; i < local_rulers.size(); i++)
		{
			std::int32_t targetPE = calculate_targetPE(mst[local_rulers[i]]);
			std::int32_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			s_rec[i] = recv_answers[packet_index];
			r_rec[i] = del[local_rulers[i]];
			
			if (r_rec[i] == 0)
				s_rec[i] = i + prefix_sum_num_vertices_per_PE[rank];
		}
		
		/*
		for (int i = 0; i < local_rulers.size(); i++)
			std::cout << i + prefix_sum_num_vertices_per_PE[rank] << " s_rec:" << s_rec[i] << ", r_rec:" << r_rec[i] << std::endl;
		*/
		timer.add_checkpoint("rekursion");

		wood_irregular_pointer_doubling recursion(s_rec, r_rec, targetPEs_rec, prefix_sum_num_vertices_per_PE, comm);
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
		timer.add_checkpoint("finalen_ranks_berechnen");







		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{
			std::int32_t targetPE = calculate_targetPE(mst[i]);
			num_packets_per_PE[targetPE]++;
		}
		
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		//every node i now has a ruler as checkpoint in mst array aka mst[i]. We now want from mst[i] the distance to root and root
		//there is one inderection needed to translate the recursive index of the root to the global index non recursive of the root
		struct request{
			std::uint64_t node;
			std::uint64_t ruler;
		};
		std::vector<request> request_ruler(num_local_vertices);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{
			std::int32_t targetPE = calculate_targetPE(mst[i]);
			std::uint64_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			request_ruler[packet_index].ruler = mst[i];
			request_ruler[packet_index].node = i + node_offset;
		}
		//first we send requests to mst[i], because only this PE knows the distance to root and root node from mst[i]
		auto recv_request_ruler = comm.alltoallv(kamping::send_buf(request_ruler), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
		for (std::uint64_t i = 0; i < recv_request_ruler.size(); i++)
		{
			std::int32_t targetPE = recursion.targetPEs[map_ruler_to_its_index[recv_request_ruler[i].ruler - node_offset]];
			num_packets_per_PE[targetPE]++;
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		struct forwarded_request{
			std::uint64_t node;
			std::uint64_t root; //this will be global but recursive index
			std::uint64_t distance;
		};
		
		std::vector<forwarded_request> request_ruler_forwarded(recv_request_ruler.size());
		
		for (std::uint64_t i = 0; i < recv_request_ruler.size(); i++)
		{
			std::int32_t targetPE = recursion.targetPEs[map_ruler_to_its_index[recv_request_ruler[i].ruler - node_offset]];
			std::uint64_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			request_ruler_forwarded[packet_index].node = recv_request_ruler[i].node;
			request_ruler_forwarded[packet_index].root = recursion.q[map_ruler_to_its_index[recv_request_ruler[i].ruler - node_offset]];
			request_ruler_forwarded[packet_index].distance = recursion.r[map_ruler_to_its_index[recv_request_ruler[i].ruler - node_offset]];
			
			
		}
		auto recv_request_ruler_forwarded = comm.alltoallv(kamping::send_buf(request_ruler_forwarded), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer();
		std::fill(num_packets_per_PE.begin(), num_packets_per_PE.end(), 0);
		for (std::uint64_t i = 0; i < recv_request_ruler_forwarded.size(); i++)
		{
			std::int32_t targetPE = calculate_targetPE(recv_request_ruler_forwarded[i].node);
			num_packets_per_PE[targetPE]++;
		}
		calculate_send_displacements_and_reset_num_packets_per_PE(send_displacements, num_packets_per_PE);
		struct answer{
			std::uint64_t node;
			std::uint64_t root; //this will be global index
			std::uint64_t distance;
		};
		std::vector<answer> answer_ruler(recv_request_ruler_forwarded.size());
		for (std::uint64_t i = 0; i < recv_request_ruler_forwarded.size(); i++)
		{
			std::int32_t targetPE = calculate_targetPE(recv_request_ruler_forwarded[i].node);
			std::uint64_t packet_index = send_displacements[targetPE] + num_packets_per_PE[targetPE]++;
			answer_ruler[packet_index].node = recv_request_ruler_forwarded[i].node;
			answer_ruler[packet_index].distance = recv_request_ruler_forwarded[i].distance;
			answer_ruler[packet_index].root = local_rulers[recv_request_ruler_forwarded[i].root - prefix_sum_num_vertices_per_PE[rank]] + node_offset;
			/*
			std::cout << rank << " bekommt request (node,root): (" << recv_request_ruler_forwarded[i].node << "," << recv_request_ruler_forwarded[i].root 
			<< ") und wird weitergeleitet an " << targetPE << " mit : (" << answer_ruler[packet_index].node << "," 
			<< answer_ruler[packet_index].root << ")" << std::endl;*/
		}
		auto recv_ruler_answers = comm.alltoallv(kamping::send_buf(answer_ruler), kamping::send_counts(num_packets_per_PE)).extract_recv_buffer(); //size = num_local_vertices
		
		
		
	
		result_dist = std::vector<std::int64_t>(num_local_vertices);
		result_root = std::vector<std::uint64_t>(num_local_vertices);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{
			std::uint64_t local_index = recv_ruler_answers[i].node - node_offset;
			result_root[local_index] = recv_ruler_answers[i].root;
			result_dist[local_index] = recv_ruler_answers[i].distance + del[local_index];
		}
		
		timer.finalize(comm, num_local_vertices, comm_rounds);

		/*
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
			std::cout << i + node_offset << " has distance " << results[i].distance << " of its root " << results[i].root << std::endl;
		*/
		
	}
	
	

	std::int32_t calculate_targetPE(std::uint64_t global_index)
	{
		return unmask(global_index) / num_local_vertices;
	}
	
	std::string packet_to_string(packet packet)
	{
		return "(" + std::to_string(packet.ruler_source) + "," + std::to_string(packet.destination) + "," + std::to_string(packet.distance) + ")";
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
	
	void mark_as_ruler(std::uint64_t local_index)
	{
		bounds[local_index] =  mark(bounds[local_index],0);
	}
	
	bool is_ruler(std::uint64_t local_index)
	{
		return is_marked(bounds[local_index],0);
	}
	
	void mark_as_reached(std::uint64_t local_index)
	{
		bounds[local_index] =  mark(bounds[local_index],1);
	}
	
	bool is_reached(std::uint64_t local_index)
	{
		return is_marked(bounds[local_index],1);
	}
	
	bool is_leaf(std::uint64_t local_index)
	{
		return unmask(bounds[local_index]) == unmask(bounds[local_index + 1]);
	}
	
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
	
	
	public:
	std::uint64_t node_offset;
	std::uint64_t num_local_vertices;
	std::uint64_t rank, size;
	std::vector<std::uint64_t> edges;
	std::vector<std::uint64_t> bounds;
	
	std::vector<std::int64_t> result_dist;
	std::vector<std::uint64_t> result_root;
	
	std::uint64_t comm_rounds;
};


	