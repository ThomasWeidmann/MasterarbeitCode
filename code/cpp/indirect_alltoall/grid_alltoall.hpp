#include "grid_all_to_all_impl.hpp"

template <typename SendBuffer, typename DestinationHandler>
auto grid_mpi_all_to_all(
  SendBuffer const&       send_buf,
  DestinationHandler&&    get_final_destination,
  GridCommunicator const& grid_comm
) {
  auto mpi_result_rowwise  = rowwise_exchange<true>(send_buf, get_final_destination, grid_comm);
  auto rowwise_recv_buf    = mpi_result_rowwise.extract_recv_buffer();
  auto rowwise_recv_counts = mpi_result_rowwise.extract_recv_counts();
  auto rowwise_recv_displs = mpi_result_rowwise.extract_recv_displs();

  return columnwise_exchange(rowwise_recv_buf, grid_comm);
}
