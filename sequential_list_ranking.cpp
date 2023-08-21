class sequential_list_ranking
{
	public:
	
	//if this PE has final node, then final node is set to a valid value, otherweise it is -1
	sequential_list_ranking(std::vector<std::int32_t>& successors)
	{
		s = successors;
	}
	
	
	std::vector<std::int32_t> start(kamping::Communicator<>& comm)
	{

		
		std::int32_t num_global_vertices = s.size();
	
		std::vector<bool> has_pred(num_global_vertices, 0);
		for (std::int32_t i=0; i < num_global_vertices; i++)
			has_pred[s[i]] = true;
		
		std::int32_t start_node=0;
		for (std::int32_t i=0; i < num_global_vertices; i++)
		{
			if (!has_pred[i])
			{
				start_node = i;
			}
		}


		
		std::vector<std::int32_t> r(num_global_vertices,0);
		std::int32_t rank = num_global_vertices - 1;
		std::int32_t node = start_node;
		for (int i = 0; i < num_global_vertices-1; i++)
		{
			r[node] = rank--;
			node = s[node];
			
		}
		
	
		
		return r;
	}
	
	
	private:

	std::vector<std::int32_t> s;
};
		