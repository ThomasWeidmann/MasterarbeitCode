#include <fmt/format.h>
#include <fmt/ranges.h>
#include "CLI11.hpp"
//#include <CLI/CLI.hpp>
#include <message-queue/buffered_queue.hpp>
#include <random>

class example
{
	public:
	
	example()
	{
		
	}
	
	
	
	void test(std::vector<std::uint64_t>& s, std::uint64_t dist_rulers, kamping::Communicator<>& comm)
	{
		this->s = s;
		
		struct packet{
			std::uint64_t ruler_source;
			std::uint64_t destination;
			std::uint32_t distance;
		};
		
		
		rank = comm.rank();
		size = comm.size();
		num_local_vertices = s.size();
		node_offset = rank * num_local_vertices;
		
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
		
		//queue initialized
		auto printing_cleaner = [](auto& buf, message_queue::PEID receiver) {
			message_queue::atomic_debug(fmt::format("Preparing buffer {} to {}.", buf, receiver));
		};

		//packet = {ruler source, destination} because compile problems
		auto queue = message_queue::make_buffered_queue<std::pair<std::uint64_t,std::uint64_t>>(MPI_COMM_WORLD, printing_cleaner);		
		
		auto on_message = [&](message_queue::Envelope<std::pair<std::uint64_t,std::uint64_t>> auto envelope) {
            if (true) {
                auto begin = envelope.message.begin();
                std::stringstream ss;
                ss << "Message " << *(begin + 2) << " from " << *begin << " arrived after " << *(begin + 1) << " hops.";
                message_queue::atomic_debug(ss.str());
            } else {
                envelope.message[1]++;
                queue.post_message(std::move(envelope.message), 0);
            }
        };
		
		//
		
		for (std::uint64_t i = 0; i < out_buffer_size; i++)
		{
			while (is_final(ruler_index)) ruler_index++;
			
			local_rulers.push_back(ruler_index);
			rulers_to_send_packages[i] = ruler_index;
			std::int32_t targetPE = calculate_targetPE(s[ruler_index]);
			num_packets_per_PE[targetPE]++;
			ruler_index++;
			std::cout << rank << " posts message (" << ruler_index + node_offset << "," << s[ruler_index] << ") to PE " << targetPE << std::endl;
			queue.post_message(std::pair{ruler_index + node_offset, s[ruler_index]}, targetPE);
			
		}
	}
	
	
	
	
	bool is_final(std::uint64_t local_index)
	{
		return local_index + node_offset == unmask(s[local_index]);
	}
	
	std::int32_t calculate_targetPE(std::uint64_t global_index)
	{
		return unmask(global_index) / num_local_vertices;
	}
	
	std::uint64_t unmask(std::uint64_t value)
	{
		return value & 0xfffffffffffffff;
	}
	
	private:
	
	std::uint64_t node_offset;
	std::uint64_t num_local_vertices;
	std::uint64_t rank, size;
	std::vector<std::uint64_t> s;
	std::uint64_t dist_rulers;
	std::uint32_t num_iterations;


};