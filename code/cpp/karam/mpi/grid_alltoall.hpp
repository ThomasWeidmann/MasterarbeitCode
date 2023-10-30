#include <algorithm>
#include <numeric>

#include <kamping/collectives/alltoall.hpp>
#include <kamping/collectives/barrier.hpp>
#include <kamping/communicator.hpp>

#include "indirect_all_to_all/grid_all_to_all_impl.hpp"

namespace karam::mpi {

template <typename SendBuffer, typename DestinationHandler>
auto grid_mpi_all_to_all(
  SendBuffer const&       send_buf,
  DestinationHandler&&    get_final_destination,
  karam::mpi::GridCommunicator const& grid_comm
) {
  auto mpi_result_rowwise  = rowwise_exchange<true>(send_buf, get_final_destination, grid_comm);
  auto rowwise_recv_buf    = mpi_result_rowwise.extract_recv_buffer();
  auto rowwise_recv_counts = mpi_result_rowwise.extract_recv_counts();
  auto rowwise_recv_displs = mpi_result_rowwise.extract_recv_displs();

  return columnwise_exchange(rowwise_recv_buf, grid_comm);
}
//from here on my code

template <typename SendBuffer>
auto grid_mpi_all_to_all2(
  SendBuffer const&       send_buf,
  std::vector<std::int32_t>    send_counts,
  karam::mpi::GridCommunicator const& grid_comm,
  kamping::Communicator<>& comm
) {
  //auto mpi_result_rowwise  = rowwise_exchange<true>(send_buf, get_final_destination, grid_comm);
  //auto rowwise_recv_buf    = mpi_result_rowwise.extract_recv_buffer();
  //auto rowwise_recv_counts = mpi_result_rowwise.extract_recv_counts();
  //auto rowwise_recv_displs = mpi_result_rowwise.extract_recv_displs();

  //return columnwise_exchange(rowwise_recv_buf, grid_comm);
}

template <typename SendBuffer>
auto my_grid_all_to_all(
  SendBuffer const&       send_buf,
  std::vector<std::int32_t>    send_counts,
  karam::mpi::GridCommunicator const& grid_comm,
  kamping::Communicator<>& comm
) {
	using T                                               = typename SendBuffer::value_type;
	 
	std::vector<std::int32_t> send_counts_row = std::vector<std::int32_t>(grid_comm.row_comm().size(),0);
	for (std::int32_t p = 0; p < comm.size(); p++)
	{
		std::int32_t targetPE = grid_comm.proxy_col_index(static_cast<std::size_t>(p));
		send_counts_row[targetPE] += send_counts[p];
	}
	
	auto send_displacements = send_counts_row;
	std::exclusive_scan(send_counts_row.begin(), send_counts_row.end(), send_displacements.begin(), 0ull);
	auto       index_displacements = send_displacements;
	
	utils::default_init_vector<IndirectMessage<T>> contiguous_send_buf(send_buf.size()); 
	
	std::uint64_t index = 0;
	for (std::int32_t p = 0; p < comm.size(); p++)
	{

		for (std::uint64_t i = 0; i < send_counts[p]; i++)
		{
			auto const final_destination = p;
			auto const destination_in_row = grid_comm.proxy_col_index(static_cast<std::size_t>(final_destination));
			auto const idx = index_displacements[destination_in_row]++;

			contiguous_send_buf[static_cast<std::size_t>(idx)] = IndirectMessage<T>(
				static_cast<std::uint32_t>(comm.rank()),
				static_cast<std::uint32_t>(final_destination),
				send_buf[index++]
			  );
			
			
		}
		
	}

	auto mpi_result_rowwise = grid_comm.row_comm().alltoallv(
		kamping::send_buf(contiguous_send_buf),
		kamping::send_counts(send_counts_row)
	);

  //auto mpi_result_rowwise  = rowwise_exchange<true>(send_buf, get_final_destination, grid_comm);
	auto rowwise_recv_buf    = mpi_result_rowwise.extract_recv_buffer();


	return columnwise_exchange(rowwise_recv_buf, grid_comm);
}


template <bool use_indirect_wrapper = true, typename SendBuffer, typename DestinationHandler>
auto my_rowwise_exchange(
  SendBuffer const&       send_buf,
  DestinationHandler&&    get_final_destination,
  GridCommunicator const& grid_comm
) {
  using T                                               = typename SendBuffer::value_type;
  kamping::Communicator<> const& comm                   = kamping::comm_world();
  auto                           get_destination_in_row = [&](auto const& elem) {
    auto final_destination = get_final_destination(elem);
    return grid_comm.proxy_col_index(static_cast<std::size_t>(final_destination));
  };
  auto const send_counts = compute_send_counts_from_simple_buffer(
    send_buf,
    get_destination_in_row,
    grid_comm.row_comm().size()
  );
  auto send_displacements = send_counts;
  std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displacements.begin(), 0ull);
  auto       index_displacements = send_displacements;
  auto const total_send_count =
    static_cast<std::size_t>(send_displacements.back() + send_counts.back());

  using MsgType = std::conditional_t<use_indirect_wrapper, IndirectMessage<T>, T>;
  std::cout << comm.rank() << " with " << total_send_count << std::endl;
  utils::default_init_vector<MsgType> contiguous_send_buf(total_send_count);
  for (auto const& elem: send_buf) {
    auto const final_destination = get_final_destination(elem);
    auto const destination_in_row =
      grid_comm.proxy_col_index(static_cast<std::size_t>(final_destination));
    auto const idx = index_displacements[destination_in_row]++;
    if constexpr (use_indirect_wrapper) {
      contiguous_send_buf[static_cast<std::size_t>(idx)] = MsgType(
        static_cast<std::uint32_t>(comm.rank()),
        static_cast<std::uint32_t>(final_destination),
        elem
      );
    } else {
      contiguous_send_buf[static_cast<std::size_t>(idx)] = elem;
    }
  }
	
std::cout << comm.rank() << " with send_buf:";
for (int i = 0; i < contiguous_send_buf.size(); i++)
	std::cout << contiguous_send_buf[i] << " ";
std::cout << std::endl;

std::cout << comm.rank() << " with send_counts:";
for (int i = 0; i < send_counts.size(); i++)
	std::cout << send_counts[i] << " ";
std::cout << std::endl;
	
  return grid_comm.row_comm().alltoallv(
    kamping::send_buf(contiguous_send_buf),
    kamping::send_counts(send_counts)
  );
}


}



