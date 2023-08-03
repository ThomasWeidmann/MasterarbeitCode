
#include <math.h> 
#include <kagen.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

int rank, size;

struct two_int {
	std::int32_t i1;
	std::int32_t i2;
};

struct three_int {
	std::int32_t i1;
	std::int32_t i2;
	std::int32_t i3;
};

struct unidirectional_path {
  std::vector<std::int32_t> s;
  std::int32_t num_global_vertices;
  std::int32_t num_local_vertices;
  std::vector<std::int32_t> num_vertices_per_pe;
  std::vector<std::int32_t> prefix_sum_num_vertices_per_pe;
};

struct bidirectional_path {
  std::vector<std::int32_t> s;
  std::vector<std::int32_t> p;
  std::int32_t num_global_vertices;
  std::int32_t num_local_vertices;
  std::vector<std::int32_t> num_vertices_per_pe;
  std::vector<std::int32_t> prefix_sum_num_vertices_per_pe;
};

std::vector<std::int32_t> independent_set_removal(bidirectional_path bidirectional_path);
void independent_set_removal_rec(bidirectional_path bidirectional_path, std::vector<std::int32_t> r);
void pointer_doubling(unidirectional_path unidirectional_path);
unidirectional_path generate_unidirectional_path(std::int32_t num_global_vertices);
bidirectional_path generate_bidirectional_path(std::int32_t num_global_vertices);
void sparse_ruling_set(unidirectional_path unidirectional_path, std::int32_t num_rulers);

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	srand((unsigned) time(NULL) + rank);
	if (argc != 2){
		std::cout << "Es gibt einen notwendigen Parameter der Anzahl der gesamten Knoten angibt!" << std::endl;
		return 0;
	}
	std::int32_t num_global_vertices = atoi(argv[1]);

	unidirectional_path unidirectional_path = generate_unidirectional_path(num_global_vertices);
	sparse_ruling_set(unidirectional_path, 2);
    MPI_Finalize();
}

std::vector<std::int32_t> independent_set_removal(bidirectional_path bidirectional_path)
{
	std::vector<std::int32_t> s = bidirectional_path.s;
	std::int32_t num_local_vertices = bidirectional_path.num_local_vertices;
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe = bidirectional_path.prefix_sum_num_vertices_per_pe;
	
	std::vector<std::int32_t> r(num_local_vertices);
	for (std::int32_t i = 0; i < num_local_vertices; i++)
	{
		r[i] = 1;
		if (s[i] == i + prefix_sum_num_vertices_per_pe[rank])
			r[i] = 0;
	}
	
	independent_set_removal_rec(bidirectional_path, r);
	
	return r;
}


std::int32_t global_index(std::int32_t local_index, std::vector<std::int32_t> prefix_sum_num_vertices_per_pe)
{
	return local_index + prefix_sum_num_vertices_per_pe[rank];
}

std::int32_t calculate_targetPE(std::int32_t global_index, std::int32_t num_global_vertices)
{
	std::int32_t lower_bound_num_local_vertices = num_global_vertices / size;
	std::int32_t upper_bound_num_local_vertices = lower_bound_num_local_vertices + 1;
	std::int32_t rest = num_global_vertices - size * lower_bound_num_local_vertices;
	
	if (global_index < rest * upper_bound_num_local_vertices)
		return global_index / upper_bound_num_local_vertices;
	
	return rest + (global_index - (rest * upper_bound_num_local_vertices)) / lower_bound_num_local_vertices;
}

void sparse_ruling_set(unidirectional_path unidirectional_path, std::int32_t dist_rulers)
{
	
	
	

	//an element is a ruler, iff its local_index % dist_rulers == 0, this way we assure every PE has same number of rulers
	std::vector<std::int32_t> s = unidirectional_path.s;
	std::int32_t num_global_vertices = unidirectional_path.num_global_vertices;
	std::int32_t num_local_vertices = unidirectional_path.num_local_vertices;
	std::vector<std::int32_t> num_vertices_per_pe = unidirectional_path.num_vertices_per_pe;
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe = unidirectional_path.prefix_sum_num_vertices_per_pe;
	
	std::int32_t num_local_rulers = (num_local_vertices - 1) / dist_rulers + 1;
	//length of request between any two PES, expect value is num_local_rulers/ size
	
	
	std::vector<std::int32_t> targetPEs(num_local_vertices);
	for (std::int32_t i = 0; i < num_local_vertices; i++)
		targetPEs[i] = calculate_targetPE(s[i], num_global_vertices);
	
	//first all rulers write their packets in out_buffer
	uint64_t start_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

	std::vector<std::int32_t> sizes_out_buffer(size);
	for (std::int32_t local_ruler_index = 0; local_ruler_index < num_local_vertices; local_ruler_index += dist_rulers)
	{
		std::int32_t targetPE = targetPEs[local_ruler_index];
		sizes_out_buffer[targetPE]++;	
	}
	
	std::vector<std::int32_t> sizes_recv_buffer(size);
	MPI_Alltoall(&sizes_out_buffer[0], 1, MPI_INT, &sizes_recv_buffer[0], 1, MPI_INT, MPI_COMM_WORLD);
	std::vector<two_int> linear_recv_buffer(num_local_vertices);
	std::vector<two_int*> recv_buffer(size);
	std::vector<std::int32_t> displacements_recv_buffer(size);
	displacements_recv_buffer[0]=0;
	recv_buffer[0] = &linear_recv_buffer[0];
	for (std::int32_t i = 1; i < size; i++)
	{
		recv_buffer[i] = recv_buffer[i-1] + sizes_recv_buffer[i-1];
		displacements_recv_buffer[i] = displacements_recv_buffer[i-1] + sizes_recv_buffer[i-1];
	}
	std::vector<two_int> linear_out_buffer(num_local_vertices);
	std::vector<two_int*> out_buffer(size);
	std::vector<std::int32_t> displacements_out_buffer(size);
	
	displacements_out_buffer[0] = 0;
	out_buffer[0] = &linear_out_buffer[0];
	for (std::int32_t i = 1; i < size; i++)
	{
		out_buffer[i] = out_buffer[i-1] + sizes_out_buffer[i-1];
		displacements_out_buffer[i] = displacements_out_buffer[i-1] + sizes_out_buffer[i-1];
	}
	std::fill(sizes_out_buffer.begin(), sizes_out_buffer.end(), 0);
	for (std::int32_t local_ruler_index = 0; local_ruler_index < num_local_vertices; local_ruler_index += dist_rulers)
	{
		std::int32_t targetPE = targetPEs[local_ruler_index];
		out_buffer[targetPE][sizes_out_buffer[targetPE]].i1 = s[local_ruler_index];
		out_buffer[targetPE][sizes_out_buffer[targetPE]].i2 = global_index(local_ruler_index, prefix_sum_num_vertices_per_pe);
		sizes_out_buffer[targetPE]++;	
	}
	
	MPI_Alltoallv(&linear_out_buffer[0], &sizes_out_buffer[0],
                  &displacements_out_buffer[0], MPI_LONG, &linear_recv_buffer[0],
                  &sizes_recv_buffer[0], &displacements_recv_buffer[0], MPI_LONG, MPI_COMM_WORLD);
	
	/*
	if (rank  == 0)
	{
		std::cout << "sizes recv buffer:\n";
		for (int i = 0; i < size; i++)
			std::cout << i << " " << sizes_recv_buffer << "\n";
		std::cout << std::endl;
	}
	
	
	
	std::vector<two_int> linear_recv_buffer(num_local_vertices);

	
	
	
	if (rank  == 0)
	{
		std::cout << "adress recv buffer:\n";
		for (int i = 0; i < size; i++)
			std::cout << i << " " << sizes_recv_buffer << "\n";
		std::cout << std::endl;
	}
	
	
	if (rank == 0)
		std::cout << out_buffer[0] << std::endl;*/
		
	if (rank == 0)  std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - start_time << std::endl;

}

void independent_set_removal_rec(bidirectional_path bidirectional_path, std::vector<std::int32_t> r)
{
	std::vector<std::int32_t> s = bidirectional_path.s;
	std::vector<std::int32_t> p = bidirectional_path.p;
	std::int32_t num_global_vertices = bidirectional_path.num_global_vertices;
	std::int32_t num_local_vertices = bidirectional_path.num_local_vertices;
	std::vector<std::int32_t> num_vertices_per_pe = bidirectional_path.num_vertices_per_pe;
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe = bidirectional_path.prefix_sum_num_vertices_per_pe;
	
	
	std::vector<std::int32_t> c(num_local_vertices); //last bit of c[i] is coin flip, second last bit of c[i] is 1 iff i \in i
	
	
	for (uint32_t i = 0; i < num_local_vertices; i++)
		c[i] = rand() % 2;
	
	
	std::vector<std::int32_t> c_global(num_global_vertices);
	std::vector<std::int32_t> s_global(num_global_vertices);
	std::vector<std::int32_t> p_global(num_global_vertices);
	std::vector<std::int32_t> r_global(num_global_vertices);

	
	MPI_Allgatherv(&*c.begin(), num_local_vertices, MPI_INT, &*c_global.begin(), &*num_vertices_per_pe.begin(), &*prefix_sum_num_vertices_per_pe.begin(), MPI_INT, MPI_COMM_WORLD);
	MPI_Allgatherv(&*s.begin(), num_local_vertices, MPI_INT, &*s_global.begin(), &*num_vertices_per_pe.begin(), &*prefix_sum_num_vertices_per_pe.begin(), MPI_INT, MPI_COMM_WORLD);
	MPI_Allgatherv(&*p.begin(), num_local_vertices, MPI_INT, &*p_global.begin(), &*num_vertices_per_pe.begin(), &*prefix_sum_num_vertices_per_pe.begin(), MPI_INT, MPI_COMM_WORLD);
	MPI_Allgatherv(&*r.begin(), num_local_vertices, MPI_INT, &*r_global.begin(), &*num_vertices_per_pe.begin(), &*prefix_sum_num_vertices_per_pe.begin(), MPI_INT, MPI_COMM_WORLD);


	
	//noch wird f von jedem PE ganz berechnet und nicht jeder berechnet seinen Teil
	std::vector<std::int32_t> f(num_global_vertices);
	std::vector<std::int32_t> i_arr(num_global_vertices); //set of nodes in I
	std::vector<std::int32_t> not_i_arr(num_global_vertices); //set of nodes not in I
	
	std::int32_t i_arr_index = 0;
	std::int32_t not_i_arr_index = 0; //i_arr_index + not_i_arr_index == num_global_vertices after following for loop
	
	uint32_t f_val = 0;
	for (uint32_t i = 0; i < num_global_vertices; i++)
	{
		if (p_global[i] != i && s_global[i] != i && (c_global[i] & 1) == 1 && c_global[s_global[i]] == 0) //iff i \in I
		{
			c_global[i] |= 2;
			i_arr[i_arr_index++] = i;
		}
		else
		{
			f[i] = f_val++;
			not_i_arr[not_i_arr_index++] = i;
		}
	}
	
	
	std::int32_t num_local_vertices_rec = 0;
	for (uint32_t i = prefix_sum_num_vertices_per_pe[rank]; i < prefix_sum_num_vertices_per_pe[rank + 1]; i++)
		if ((c_global[i] & 2) == 0)
			num_local_vertices_rec++;
	
	std::vector<std::int32_t> num_vertices_per_pe_rec(size);
	MPI_Allgather(&num_local_vertices_rec, 1, MPI_INT, &*num_vertices_per_pe_rec.begin(), 1, MPI_INT, MPI_COMM_WORLD); 
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe_rec(size + 1);
	prefix_sum_num_vertices_per_pe_rec[0]=0;
	for (std::int32_t i = 0; i < size; i++)
		prefix_sum_num_vertices_per_pe_rec[i+1] = prefix_sum_num_vertices_per_pe_rec[i] + num_vertices_per_pe_rec[i];
	std::int32_t num_global_vertices_rec = prefix_sum_num_vertices_per_pe_rec[size];
	
	std::vector<std::int32_t> s_rec(num_local_vertices_rec);
	std::vector<std::int32_t> p_rec(num_local_vertices_rec);
	std::vector<std::int32_t> r_rec(num_local_vertices_rec);
	
	for (uint32_t i = prefix_sum_num_vertices_per_pe_rec[rank]; i < prefix_sum_num_vertices_per_pe_rec[rank + 1]; i++)
	{
		if ((c_global[s_global[not_i_arr[i]]] & 2) != 0) //iff s[i] in I
		{
			s_rec[f[not_i_arr[i]]] = f[s_global[s_global[not_i_arr[i]]]];
			r_rec[f[not_i_arr[i]]] = r_global[not_i_arr[i]] +  r_global[s_global[not_i_arr[i]]];
		}
		else
		{
			s_rec[f[not_i_arr[i]]] = f[s_global[not_i_arr[i]]];
			r_rec[f[not_i_arr[i]]] = r_global[not_i_arr[i]];
		}
		
		if ((c_global[p_global[not_i_arr[i]]] & 2) != 0) //iff p[i] in I
			p_rec[f[not_i_arr[i]]] = f[p_global[p_global[not_i_arr[i]]]];
		else
			p_rec[f[not_i_arr[i]]] = f[p_global[not_i_arr[i]]];
	}
	
	/*
	bidirectional_path bidirectional_path_rec;
	bidirectional_path_rec.s = s_rec;
	bidirectional_path_rec.p = p_rec;
	bidirectional_path_rec.num_global_vertices = num_global_vertices_rec;
	bidirectional_path_rec.num_local_vertices = num_local_vertices_rec;
	bidirectional_path_rec.num_vertices_per_pe = num_vertices_per_pe_rec;
	bidirectional_path_rec.prefix_sum_num_vertices_per_pe = prefix_sum_num_vertices_per_pe_rec;
	
	independent_set_removal_rec(bidirectional_path_rec, r_rec)
	
	//hier muss jetzt nach pseudocode von sanders das r array aktualisiert werden
	*/	
	
	if (rank == 0) //prints everything
	{
		std::cout << "s array: ";
		for (int i = 0; i < num_global_vertices; i++)
			std::cout << s_global[i] << " ";
		std::cout << std::endl;
		
		std::cout << "p array: ";
		for (int i = 0; i < num_global_vertices; i++)
			std::cout << p_global[i] << " ";
		std::cout << std::endl;
		
		std::cout << "i array: ";
		for (int i = 0; i < i_arr_index; i++)
			std::cout << i_arr[i] << " ";
		std::cout << std::endl;
		
		std::cout << "not_i array: ";
		for (int i = 0; i < not_i_arr_index; i++)
			std::cout << not_i_arr[i] << " ";
		std::cout << std::endl;
		
		std::cout << "f array: ";
		for (int i = 0; i < num_global_vertices; i++)
			std::cout << f[i] << " ";
		std::cout << std::endl;
		
		std::cout << "boolean iff i in I array: ";
		for (int i = 0; i < num_global_vertices; i++)
			std::cout << ((c_global[i] & 2) != 0) << " ";
		std::cout << std::endl;
	}
}


bidirectional_path generate_bidirectional_path(std::int32_t num_global_vertices)
{
	kagen::KaGen gen(MPI_COMM_WORLD);
	
	auto path = gen.GenerateDirectedPath(num_global_vertices, true);
    std::int32_t num_local_vertices = path.vertex_range.second - path.vertex_range.first;
	

	std::vector<std::int32_t> num_vertices_per_pe(size);
	MPI_Allgather(&num_local_vertices, 1, MPI_INT, &*num_vertices_per_pe.begin(), 1, MPI_INT, MPI_COMM_WORLD); 
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe(size + 1);
	prefix_sum_num_vertices_per_pe[0]=0;
	for (std::int32_t i = 0; i < size; i++)
		prefix_sum_num_vertices_per_pe[i+1] = prefix_sum_num_vertices_per_pe[i] + num_vertices_per_pe[i];
	
	std::vector<std::int32_t> s(num_local_vertices);
	for (std::int32_t i = 0; i < num_local_vertices; i++)
		s[i] = i + prefix_sum_num_vertices_per_pe[rank]; //last edge has pointer to itself and not no pointer like kagen does
	
	for (auto const& [src, dst]: path.edges){
		s[src - prefix_sum_num_vertices_per_pe[rank]] = dst; 
	}
	
	std::vector<std::int32_t> s_global(num_global_vertices); 
	MPI_Allgatherv(&*s.begin(), num_local_vertices, MPI_INT, &*s_global.begin(), &*num_vertices_per_pe.begin(), &*prefix_sum_num_vertices_per_pe.begin(), MPI_INT, MPI_COMM_WORLD);
	std::vector<std::int32_t> p(num_local_vertices);
	
	for (std::int32_t i = 0; i < num_local_vertices; i++)
		p[i] = i + prefix_sum_num_vertices_per_pe[rank];
		
	for (std::int32_t i = 0; i < num_global_vertices; i++)
		if (i != s_global[i])
			if (path.vertex_range.first <= s_global[i] && s_global[i] < path.vertex_range.second)
				p[s_global[i] - prefix_sum_num_vertices_per_pe[rank]] = i;
	
	bidirectional_path bidirectional_path;
	bidirectional_path.s = s;
	bidirectional_path.p = p;
	bidirectional_path.num_global_vertices = num_global_vertices;
	bidirectional_path.num_local_vertices = num_local_vertices;
	bidirectional_path.num_vertices_per_pe = num_vertices_per_pe;
	bidirectional_path.prefix_sum_num_vertices_per_pe = prefix_sum_num_vertices_per_pe;
	
	return bidirectional_path;
}


unidirectional_path generate_unidirectional_path(std::int32_t num_global_vertices)
{
	kagen::KaGen gen(MPI_COMM_WORLD);
	
	auto path = gen.GenerateDirectedPath(num_global_vertices, true);
    std::int32_t num_local_vertices = path.vertex_range.second - path.vertex_range.first;
	

	std::vector<std::int32_t> num_vertices_per_pe(size);
	MPI_Allgather(&num_local_vertices, 1, MPI_INT, &*num_vertices_per_pe.begin(), 1, MPI_INT, MPI_COMM_WORLD); 
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe(size + 1);
	prefix_sum_num_vertices_per_pe[0]=0;
	for (std::int32_t i = 0; i < size; i++)
		prefix_sum_num_vertices_per_pe[i+1] = prefix_sum_num_vertices_per_pe[i] + num_vertices_per_pe[i];
	
	std::vector<std::int32_t> s(num_local_vertices);
	for (std::int32_t i = 0; i < num_local_vertices; i++)
		s[i] = i + prefix_sum_num_vertices_per_pe[rank]; //last edge has pointer to itself and not no pointer like kagen does

	for (auto const& [src, dst]: path.edges){
		s[src - prefix_sum_num_vertices_per_pe[rank]] = dst; 
	}

	unidirectional_path unidirectional_path;
	unidirectional_path.s = s;
	unidirectional_path.num_global_vertices = num_global_vertices;
	unidirectional_path.num_local_vertices = num_local_vertices;
	unidirectional_path.num_vertices_per_pe = num_vertices_per_pe;
	unidirectional_path.prefix_sum_num_vertices_per_pe = prefix_sum_num_vertices_per_pe;
	
	return unidirectional_path;
}

void pointer_doubling(unidirectional_path unidirectional_path)
{
	std::vector<std::int32_t> s = unidirectional_path.s;
	std::int32_t num_global_vertices = unidirectional_path.num_global_vertices;
	std::int32_t num_local_vertices = unidirectional_path.num_local_vertices;
	std::vector<std::int32_t> num_vertices_per_pe = unidirectional_path.num_vertices_per_pe;
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe = unidirectional_path.prefix_sum_num_vertices_per_pe;
	
	if (size % 2 != 0)
	{
		std::cout << "pointer doubling is just defined for even numer of PEs!" << std::endl;
		return;
	}
	std::vector<std::int32_t> value_to_PE(num_global_vertices);
	std::int32_t position = 0;
	for (std::int32_t i = 0; i < size; i++)
		for (std::int32_t j = 0; j < num_vertices_per_pe[i]; j++)
			value_to_PE[position++] = i;
		
		
	two_int r_q[num_local_vertices]; //i1 is r, i2 is q
	for (std::int32_t i = 0; i < num_local_vertices; i++)
	{
		r_q[i].i1 = 1;
		if (s[i] == i + prefix_sum_num_vertices_per_pe[rank])
			r_q[i].i1 = 0;
		r_q[i].i2 =s[i];
	}
	
	std::int32_t max_number_of_request_per_PE = num_global_vertices; //is the max length of request between any two PES, expect value is num_global_vertices/ (size^2)
	std::vector<std::int32_t> num_elements_round1(size, 0);
	std::vector<std::vector<three_int>> round1(size, std::vector<three_int>(max_number_of_request_per_PE));
	std::vector<std::int32_t> num_elements_round2(size, 0);
	std::vector<std::vector<three_int>> round2(size, std::vector<three_int>(max_number_of_request_per_PE));
	std::vector<std::int32_t> num_elements_round3(size, 0);
	std::vector<std::vector<three_int>> round3(size, std::vector<three_int>(max_number_of_request_per_PE));
	std::vector<std::int8_t> active_vertices(num_local_vertices, 0); //active_vertices[i] == 0 iff local node i is active
	three_int buffer[max_number_of_request_per_PE];
	
	/*testing
	std::cout << "PE " << rank << " has following edges values\n";
		for (std::int32_t i = 0; i < num_local_vertices; i++)
			std::cout << "(" << i + prefix_sum_num_vertices_per_pe[rank] << "," << r_q[i].i2 << "),";
	std::cout << std::endl;
	*/
	
	
	/*testing
	std::cout << "PE " << rank << " requests the following packets\n";
	for (std::int32_t i = 0; i < size; i++)
	{
		std::cout << "from PE " << i << ":";
		for (std::int32_t j = 0; j < num_elements_round1[i]; j++)
			std::cout << "(" << round1[i][j].i1 << "," << round1[i][j].i2 << "," << round1[i][j].i3 << "),";
		std::cout << std::endl;
	}*/
	uint64_t start_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	uint64_t communication_time = 0;
	
	std::int32_t num_iterations =( log(num_global_vertices) / log(3)) + 2;
	
	for (int iteration = 0; iteration < num_iterations; iteration++)
	{
		
		
		std::fill(num_elements_round1.begin(), num_elements_round1.end(), 0);
		std::fill(num_elements_round2.begin(), num_elements_round2.end(), 0);
		std::fill(num_elements_round3.begin(), num_elements_round3.end(), 0);
		
		//packages for round1
		for (std::int32_t i = 0; i < num_local_vertices; i++)
		{
			if (active_vertices[i] != 0)
				continue;
			
			std::int32_t targetPE = value_to_PE[r_q[i].i2];
		
			round1[targetPE][num_elements_round1[targetPE]] = {i + prefix_sum_num_vertices_per_pe[rank], r_q[i].i2, r_q[i].i1};
		/*
			round1[targetPE][num_elements_round1[targetPE]].i1 = i + prefix_sum_num_vertices_per_pe[rank];
			round1[targetPE][num_elements_round1[targetPE]].i2 = r_q[i].i2;
			round1[targetPE][num_elements_round1[targetPE]].i3 = r_q[i].i1;*/
			num_elements_round1[targetPE]++;
			
		}
		for (std::int32_t i = 0; i < size; i++)
		{
			round1[i][num_elements_round1[i]].i1 = -1; //this is end of list
			num_elements_round1[i]++;
		}
		
			
		//1 factor algorithm for round1
		for (std::int32_t i = 0; i < size - 1; i++)
		{
			std::int32_t partner;
			std::int32_t idle = (size * ((size - 2)/2 +1) * i) % (size - 1);
			if (rank == size - 1)
				partner = idle;
			else if (rank == idle)
				partner = size - 1;
			else 
				partner = (i - rank + size - 1) % (size - 1);
			
	
			uint64_t local_communication_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
			if (rank < partner)
			{
				MPI_Send(&round1[partner][0], 3 * num_elements_round1[partner], MPI_INT, partner, 0, MPI_COMM_WORLD);
				MPI_Recv(&buffer[0], 3 * max_number_of_request_per_PE, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			else
			{
				MPI_Recv(&buffer[0], 3 * max_number_of_request_per_PE, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(&round1[partner][0], 3 * num_elements_round1[partner], MPI_INT, partner, 0, MPI_COMM_WORLD);
			}
			communication_time += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - local_communication_time;

			/*testing
			if (iteration == 1){
			for (int i = 0; i < num_elements_round1[partner] - 1; i++)
				std::cout << "ROUND 1: PE " << rank << " sent following packet to PE " << partner << " with following triple ("<< round1[partner][i].i1 <<","<< round1[partner][i].i2<< ","<< round1[partner][i].i3 <<  ")" << std::endl;
		
		
			for (int i = 0; buffer[i].i1 >= 0; i++)
				std::cout << "ROUND 1: PE " << rank << " received following packet from PE " << partner << " with following triple ("<< buffer[i].i1 <<","<< buffer[i].i2<< ","<< buffer[i].i3 <<  ")" << std::endl;
			}*/
			
			
			//prepare round2
			for (std::int32_t i = 0; buffer[i].i1 >= 0; i++)
			{
				std::int32_t local_index = buffer[i].i2 - prefix_sum_num_vertices_per_pe[rank];
				std::int32_t target_node = r_q[local_index].i2; //aka q(q(i))
				std::int32_t targetPE = value_to_PE[target_node];
				
				round2[targetPE][num_elements_round2[targetPE]] = {buffer[i].i1, target_node, buffer[i].i3 + r_q[local_index].i1};
				/*
				round2[targetPE][num_elements_round2[targetPE]].i1 = buffer[i].i1;
				round2[targetPE][num_elements_round2[targetPE]].i2 = target_node;
				round2[targetPE][num_elements_round2[targetPE]].i3 = buffer[i].i3 + r_q[local_index].i1;
				*/
				num_elements_round2[targetPE]++;

			}
			
		}
		
		//requests for itself
		for (std::int32_t i = 0; i < num_elements_round1[rank] - 1; i++)
		{
			std::int32_t local_index = round1[rank][i].i2 - prefix_sum_num_vertices_per_pe[rank];
			std::int32_t target_node = r_q[local_index].i2; //aka q(q(i))
			std::int32_t targetPE = value_to_PE[target_node];
			
			round2[targetPE][num_elements_round2[targetPE]].i1 = round1[rank][i].i1;
			round2[targetPE][num_elements_round2[targetPE]].i2 = target_node;
			round2[targetPE][num_elements_round2[targetPE]].i3 = round1[rank][i].i3 + r_q[local_index].i1;
			
			num_elements_round2[targetPE]++;
			
		}
		for (std::int32_t i = 0; i < size; i++)
			round2[i][num_elements_round2[i]++].i1 = -1; //this is end of list
		/*testing
		std::cout << "test ROUND 2: PE " << rank << " requests the following packets\n";
		for (std::int32_t i = 0; i < size; i++)
		{
			std::cout << "test ROUND 2: from PE " << i << ":";
			for (std::int32_t j = 0; j < num_elements_round2[i]; j++)
				std::cout << "(" << round2[i][j].i1 << "," << round2[i][j].i2 << "," << round2[i][j].i3 << "),";
			std::cout << std::endl;
		}*/
		
		
		
		//1 factor algorithm for round2
		for (std::int32_t i = 0; i < size - 1; i++)
		{
			std::int32_t partner;
			std::int32_t idle = (size * ((size - 2)/2 +1) * i) % (size - 1);
			if (rank == size - 1)
				partner = idle;
			else if (rank == idle)
				partner = size - 1;
			else 
				partner = (i - rank + size - 1) % (size - 1);
			
			uint64_t local_communication_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
			if (rank < partner)
			{
				MPI_Send(&round2[partner][0], 3 * num_elements_round2[partner], MPI_INT, partner, 0, MPI_COMM_WORLD);
				MPI_Recv(&buffer[0], 3 * max_number_of_request_per_PE, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			else
			{
				MPI_Recv(&buffer[0], 3 * max_number_of_request_per_PE, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(&round2[partner][0], 3 * num_elements_round2[partner], MPI_INT, partner, 0, MPI_COMM_WORLD);
			}
			communication_time += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - local_communication_time;

			/*testing
			if (iteration == 1){
			for (int i = 0; i < num_elements_round2[partner]-1; i++)
				std::cout << "ROUND 2: PE " << rank << " sent following packet to PE " << partner << " with following triple ("<< round2[partner][i].i1 <<","<< round2[partner][i].i2<< ","<< round2[partner][i].i3 <<  ")" << std::endl;
			
		
			for (int i = 0; buffer[i].i1 >= 0; i++)
				std::cout << "ROUND 2: PE " << rank << " received following packet from PE " << partner << " with following triple ("<< buffer[i].i1 <<","<< buffer[i].i2<< ","<< buffer[i].i3 <<  ")" << std::endl;
			}*/
			
			//prepare round3
			for (std::int32_t i = 0; buffer[i].i1 >= 0; i++)
			{
				std::int32_t local_index = buffer[i].i2 - prefix_sum_num_vertices_per_pe[rank];
				std::int32_t target_node = r_q[local_index].i2; //aka q(q(q(i)))
				std::int32_t targetPE = value_to_PE[buffer[i].i1];
				
				round3[targetPE][num_elements_round3[targetPE]] = {buffer[i].i1, target_node, buffer[i].i3 + r_q[local_index].i1};
				/*
				round3[targetPE][num_elements_round3[targetPE]].i1 = buffer[i].i1;
				round3[targetPE][num_elements_round3[targetPE]].i2 = target_node;
				round3[targetPE][num_elements_round3[targetPE]].i3 = buffer[i].i3 + r_q[local_index].i1;*/
				num_elements_round3[targetPE]++;
			}
			
		}
		
		//requests for itself
		for (std::int32_t i = 0; i < num_elements_round2[rank] - 1; i++)
		{
			std::int32_t local_index = round2[rank][i].i2 - prefix_sum_num_vertices_per_pe[rank];
			std::int32_t target_node = r_q[local_index].i2; //aka q(q(q(i)))
			std::int32_t targetPE = value_to_PE[round2[rank][i].i1];
			
			round3[targetPE][num_elements_round3[targetPE]].i1 = round2[rank][i].i1;
			round3[targetPE][num_elements_round3[targetPE]].i2 = target_node;
			round3[targetPE][num_elements_round3[targetPE]].i3 = round2[rank][i].i3 + r_q[local_index].i1;
			
			
			num_elements_round3[targetPE]++;
		}
		
		/*testing
		if (iteration == 1){
		std::cout << "ROUND 3: PE " << rank << " requests the following packets\n";
		for (std::int32_t i = 0; i < size; i++)
		{
			std::cout << "ROUND 3: from PE " << i << ":";
			for (std::int32_t j = 0; j < num_elements_round3[i]; j++)
				std::cout << "(" << round3[i][j].i1 << "," << round3[i][j].i2 << "," << round3[i][j].i3 << "),";
			std::cout << std::endl;
		}
		}*/
		
		for (std::int32_t i = 0; i < size; i++)
			round3[i][num_elements_round3[i]++].i1 = -1; //this is end of list

		
		//1 factor algorithm for round3
		for (std::int32_t i = 0; i < size - 1; i++)
		{
			std::int32_t partner;
			std::int32_t idle = (size * ((size - 2)/2 +1) * i) % (size - 1);
			if (rank == size - 1)
				partner = idle;
			else if (rank == idle)
				partner = size - 1;
			else 
				partner = (i - rank + size - 1) % (size - 1);
			
			uint64_t local_communication_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
			if (rank < partner)
			{
				MPI_Send(&round3[partner][0], 3 * num_elements_round3[partner], MPI_INT, partner, 0, MPI_COMM_WORLD);
				MPI_Recv(&buffer[0], 3 * max_number_of_request_per_PE, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			else
			{
				MPI_Recv(&buffer[0], 3 * max_number_of_request_per_PE, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(&round3[partner][0], 3 * num_elements_round3[partner], MPI_INT, partner, 0, MPI_COMM_WORLD);
			}
			communication_time += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - local_communication_time;

			/*testing
			if (iteration == 1){
			for (int i = 0; i < num_elements_round3[partner]-1; i++)
				std::cout << "ROUND 3: PE " << rank << " sent following packet to PE " << partner << " with following triple ("<< round3[partner][i].i1 <<","<< round3[partner][i].i2<< ","<< round3[partner][i].i3 <<  ")" << std::endl;
			
		
			for (int i = 0; buffer[i].i1 >= 0; i++)
				std::cout << "ROUND 3: PE " << rank << " received following packet from PE " << partner << " with following triple ("<< buffer[i].i1 <<","<< buffer[i].i2<< ","<< buffer[i].i3 <<  ")" << std::endl;
			}*/
			
			//register values
			for (std::int32_t i = 0; buffer[i].i1 >= 0; i++)
			{
				std::int32_t local_index = buffer[i].i1 - prefix_sum_num_vertices_per_pe[rank];
				
				if (r_q[local_index].i1 == buffer[i].i3)
					active_vertices[local_index] = 1;
				
				r_q[local_index].i1 = buffer[i].i3;
				r_q[local_index].i2 = buffer[i].i2;	
				
				
			}
			
		}
		
		//requests for itself
		for (std::int32_t i = 0; i < num_elements_round3[rank] - 1; i++)
		{
			std::int32_t local_index = round3[rank][i].i1 - prefix_sum_num_vertices_per_pe[rank];
			
			if (r_q[local_index].i1 == round3[rank][i].i3)
				active_vertices[local_index] = 1;
		
			r_q[local_index].i1 = round3[rank][i].i3;
			r_q[local_index].i2 = round3[rank][i].i2;	
			
		}
		
		/*
		std::cout << "final values (node,q,r) after iteration " << iteration << " for PE " << rank << ":";
		for (int i = 0; i < num_local_vertices; i++)
			std::cout << "(" << i + prefix_sum_num_vertices_per_pe[rank] << "," << r_q[i].i2 << "," << r_q[i].i1 << "),";
		std::cout << std::endl;
		*/
	}
	if (rank == 0)  std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - start_time << " of which is " << communication_time << " communication_time" << std::endl;

}

void pointer_doubling2(unidirectional_path unidirectional_path)
{
	std::vector<std::int32_t> s = unidirectional_path.s;
	std::int32_t num_global_vertices = unidirectional_path.num_global_vertices;
	std::int32_t num_local_vertices = unidirectional_path.num_local_vertices;
	std::vector<std::int32_t> num_vertices_per_pe = unidirectional_path.num_vertices_per_pe;
	std::vector<std::int32_t> prefix_sum_num_vertices_per_pe = unidirectional_path.prefix_sum_num_vertices_per_pe;
	
	uint64_t start_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	uint64_t communication_time = 0;
	
	std::vector<std::int32_t> r_buffer(num_global_vertices);
	std::vector<std::int32_t> r(num_local_vertices);
	std::vector<std::int32_t> q_buffer(num_global_vertices);
	std::vector<std::int32_t> q(num_local_vertices);
	
	for (std::int32_t i = 0; i < num_local_vertices; i++)
	{
		r[i]=1;
		if (s[i] == i + prefix_sum_num_vertices_per_pe[rank])
			r[i]= 0;
		q[i]=s[i];
	}
	
	std::uint32_t num_iterations = std::bit_width((std::uint32_t) num_global_vertices) + 1;
	for (int i = 0; i < num_iterations; i++)
	{
		uint64_t local_communication_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		MPI_Allgatherv(&*q.begin(), num_local_vertices, MPI_INT, &*q_buffer.begin(), &*num_vertices_per_pe.begin(), &*prefix_sum_num_vertices_per_pe.begin(), MPI_INT, MPI_COMM_WORLD);
		MPI_Allgatherv(&*r.begin(), num_local_vertices, MPI_INT, &*r_buffer.begin(), &*num_vertices_per_pe.begin(), &*prefix_sum_num_vertices_per_pe.begin(), MPI_INT, MPI_COMM_WORLD);
		communication_time += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - local_communication_time;
		
		/* Wenn das hier auskommentiert wird werden das q und r array geprintet
		if (rank == 0){
			std::cout << "this is iteration " << i << std::endl;
			std::cout << "q array: ";
			for (int j = 0; j < num_global_vertices; j++)
				std::cout << q_buffer[j] << " ";
			std::cout << std::endl;
			std::cout << "r array: ";
			for (int j = 0; j < num_global_vertices; j++)
				std::cout << r_buffer[j] << " ";
			std::cout << std::endl;
		}*/

		for (std::int32_t local_index = 0; local_index < num_local_vertices; local_index++)
		{
			uint32_t global_index = local_index + prefix_sum_num_vertices_per_pe[rank];
			
			r[local_index] = r_buffer[global_index];
			if (r_buffer[q_buffer[global_index]] > 0) //is dieses if nicht unnötig?
				r[local_index] = r_buffer[global_index] + r_buffer[q_buffer[global_index]];
			q[local_index] = q_buffer[q_buffer[global_index]];
		}
	}
	
	if (rank == 0)  std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - start_time << " of which is " << communication_time << " communication_time" << std::endl;
	
}

