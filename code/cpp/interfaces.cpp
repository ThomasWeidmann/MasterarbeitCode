#pragma once

template <typename T> 
class communicator
{
	public:
	
	virtual std::vector<std::int32_t> extract_recv_counts() = 0;
	virtual std::vector<T> extract_recv_buffer() = 0;
	
	
};
