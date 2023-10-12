#pragma once

class list_ranking
{
	public:
	
	virtual std::vector<std::int64_t> start(kamping::Communicator<>& comm, std::vector<std::uint64_t>& successors) = 0;
	
	
};