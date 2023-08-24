//this class is implemented like in the paper "ultimate parallel list ranking" to minimize number of communication steps in trade of local work


class regular_ruling_set2
{
	public:
	
	regular_ruling_set2(std::vector<std::int64_t>& s, std::int32_t dist_rulers)
	{
		this->s = s;
		this->dist_rulers = dist_rulers;
	}
	
	
	std::vector<std::int64_t> start(kamping::Communicator<>& comm)
	{
		
		for (std::int32_t i = 0; i < s.size() && i < 10; i++)
			std::cout << s[i] << " ";
		std::cout << ", auf rank " << comm.rank() << std::endl;
		return s;
	
	}
	
	
	private:
	std::vector<std::int64_t> s;
	std::int32_t dist_rulers;
};
	