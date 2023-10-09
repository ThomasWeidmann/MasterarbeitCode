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
	
	timer(std::string first_checkpoint, std::vector<std::string> categories, std::string start_category)
	{

		_times = std::vector<uint64_t>(1);
		_times[0] =  get_time();

		_names = std::vector<std::string>(1);
		_names[0] = first_checkpoint;
		
		current_category = start_category;
		category_timestamp = get_time();
		
		map = std::unordered_map<std::string, uint64_t>();
		for (int i = 0; i < categories.size(); i++)
			map[categories[i]] = 0;
		
	}
	
	void switch_category(std::string category)
	{
		std::uint64_t current_time = get_time();
		map[current_category] = map[current_category] +  current_time - category_timestamp;
		category_timestamp = current_time;
		current_category = category;
		
	}
	
	void add_checkpoint(std::string checkpoint)
	{
		//better pushback

		_times.push_back(get_time());
		
		_names.push_back(checkpoint);
	}
	
	void finalize(kamping::Communicator<>& comm, std::uint64_t num_local_vertices, std::uint64_t dist_rulers)
	{

		_times.push_back(get_time());

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
			std::string output = "{\n\"num_local_vertices\":" + std::to_string(num_local_vertices) + ",\n" ; 
			output += "\"dist_rulers\":" + std::to_string(dist_rulers) + ",\n";
			output += "\"p\":" + std::to_string(comm.size()) + ",\n";
			
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

			output += "\"total time\":" + std::to_string(total_time) + ",\n";
			
			std::string all_categories = "[";
			for (const auto& [key, value] : map)
			{
				all_categories += "\"" + key + "\",";
				
				output += "\"" + key + "\":" + std::to_string(value) + ",\n";
			}
			all_categories.pop_back();
			all_categories += "]";
			
			output += "\"categories\":" + all_categories + "\n},\n";
			
			
			
			std::cout << output;
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
	std::string current_category;
	std::uint64_t category_timestamp;
	std::unordered_map<std::string, uint64_t> map;
	
	
	
	std::vector<uint64_t> _times;
	std::vector<std::string> _names;
		
	
};