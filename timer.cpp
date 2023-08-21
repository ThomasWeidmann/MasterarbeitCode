#pragma once


#include <chrono>

#include <iostream>
#include <fstream>

#include "kamping/checking_casts.hpp"
#include "kamping/collectives/alltoall.hpp"
#include "kamping/communicator.hpp"
#include "kamping/environment.hpp"
#include "kamping/collectives/allgather.hpp"
#include "kamping/collectives/bcast.hpp"
#include "kamping/collectives/gather.hpp"

#include <algorithm>
#include <cmath>
#include <vector>



class timer
{
	public:
	
	timer(std::string first_checkpoint)
	{

		_times = std::vector<uint64_t>(1);
		_times[0] =  get_time();

		_names = std::vector<std::string>(1);
		_names[0] = first_checkpoint;
		
	}
	
	void add_checkpoint(std::string checkpoint)
	{
		//better pushback

		_times.push_back(get_time());
		
		_names.push_back(checkpoint);
	}
	
	void finalize(kamping::Communicator<>& comm, std::string parameters)
	{

		_times.push_back(get_time());
		
		
		
		print(comm, parameters);
	}
	
	
	void print(kamping::Communicator<>& comm, std::string parameters)
	{
		using namespace kamping;
		
		std::vector<uint64_t> relative_times(_names.size());
		for (std::int32_t i = 0; i < _names.size(); i++)
		{
			relative_times[i] = _times[i+1] - _times[i];
		}
		
		std::vector<uint64_t> all_relative_times;
	
		comm.gather(send_buf(relative_times), recv_buf(all_relative_times), root(0));

	
		if (comm.rank() == 0)
		{
			std::string output = "{\n\"parameters\":" + parameters + ",\n" ; 
	
			std::int32_t size = comm.size();
			
			std::int32_t total_time = get_time() - _times[0];
			
			std::vector<uint64_t> all_relative_times_from_one_checkpoint(size);
			for (int i = 0; i < _names.size(); i++)
			{
				
				for (int j = 0; j < size; j++)
				{
					all_relative_times_from_one_checkpoint[j] = all_relative_times[i + j*_names.size()];
				}
				
				output += "\"" + _names[i] + "\"" + ":[" + std::to_string(all_relative_times_from_one_checkpoint[0]);
				for (int i = 1; i < size; i++)
					output += "," + std::to_string(all_relative_times_from_one_checkpoint[i]);
				output += "],\n";
					
			}

			output += "\"total time\":" + std::to_string(total_time) + "\n}\n";
			
			std::ofstream myfile;
			myfile.open ("results.txt",  std::ios::app);
			myfile << output;
			myfile.close();
		}
		
		
	}

	uint64_t get_time()
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}

	private:

	std::vector<uint64_t> _times;
	std::vector<std::string> _names;
		
	
};