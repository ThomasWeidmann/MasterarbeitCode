
class irregular_ruling_set2_rec
{
	struct packet{
		std::uint64_t ruler_source;
		std::uint64_t destination;
		std::uint64_t distance;
	};


	
	public:
	
	irregular_ruling_set2_rec(std::vector<std::uint64_t>& s, std::vector<std::uint64_t> r, std::vector<std::uint32_t> targetPEs, std::uint64_t dist_rulers, std::uint64_t node_offset)
	{
		this->s = s;
		this->dist_rulers = dist_rulers;
		this->r = r;
		this->targetPEs = targetPEs;
		this->node_offset = node_offset;
	}
	
	
	std::vector<std::uint64_t> start(kamping::Communicator<>& comm)
	{
		rank = comm.rank();
		size = comm.rank();
		num_local_vertices = s.size();
		
		//todo
		
	}
	
	
	private:
	std::uint64_t node_offset;
	std::uint64_t num_local_vertices;
	std::uint64_t rank, size;
	std::vector<std::uint64_t> s;
	std::vector<std::uint64_t> r;
	std::vector<std::uint32_t> targetPEs;
	std::uint64_t dist_rulers;
};