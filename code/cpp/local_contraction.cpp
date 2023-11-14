class local_contraction
{

	public:
	
	local_contraction()
	{

	}
	
	
	void start(kamping::Communicator<>& comm, std::vector<std::uint64_t>& s)
	{
		size = comm.size();
		rank = comm.rank();
		num_local_vertices = s.size();
		node_offset = rank * num_local_vertices;
		
		std::cout << rank << " with s arr: ";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << s[i] << " ";
		std::cout << std::endl;
		
		
		std::vector<bool> has_local_successor(num_local_vertices, false);
		std::vector<bool> has_local_predecessor(num_local_vertices, false);

		
		std::vector<std::int32_t> num_packets_per_PE(size,0);
		std::vector<std::int32_t> send_displacements(size + 1,0);
		
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{
			std::int32_t targetPE = s[i] / num_local_vertices;
			
			if (targetPE == rank)
			{
				has_local_successor[i] = true;
				has_local_predecessor[s[i] - node_offset] = true;
				//std::cout << i << " has local predecessor" << std::endl;
			}
		}
		
		struct reduced_nodes {
			std::uint64_t i;
			std::uint64_t s;
			std::int64_t r;
		};
		std::vector<reduced_nodes> nodes(0);
		for (std::uint64_t i = 0; i < num_local_vertices; i++)
		{
			
			if (!has_local_predecessor[i])
			{
				std::uint64_t succ = s[i];
				std::int64_t r = 1;
				while (((succ / num_local_vertices) == rank) && (succ != s[succ - node_offset])) 
				{
					succ = s[succ - node_offset];
					r++;
				}
				std::cout << "final node " << i + node_offset << " zeigt auf " << succ << std::endl;
			}
				//std::cout << i + node_offset << " hat keinen local predecessor" << std::endl;
			
		}

		
	}
	
	

		
	private:
	
	std::uint64_t node_offset;
	std::uint64_t num_local_vertices;
	std::int32_t rank, size;
};