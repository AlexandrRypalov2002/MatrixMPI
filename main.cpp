#include <cstring>

#include <mpi.h>

#include <matrix/Matrix.h>
#include <matrix/factory/MatrixFactory.h>

namespace
{

void count_submatrix_size(
	const math::matrix::Matrix & matrix,
	std::size_t & main_size,
	std::size_t & remain_size,
	int proc_num)
{
	std::size_t rows_num = matrix.rows();

	if (rows_num % proc_num == 0)
	{
		main_size = rows_num / proc_num;
		remain_size = main_size;

		return;
	}
}

void send_submatrix_data(double * matrix_src, std::size_t submatrix_rows, std::size_t submatrix_columns, int proc_id)
{
	MPI_Send(&submatrix_rows, 1, MPI_UNSIGNED_LONG, proc_id, 0, MPI_COMM_WORLD);
	MPI_Send(&submatrix_columns, 1, MPI_UNSIGNED_LONG, proc_id, 1, MPI_COMM_WORLD);
	MPI_Send(matrix_src, submatrix_rows * submatrix_columns, MPI_DOUBLE, proc_id, 2, MPI_COMM_WORLD);
}

void recv_submatrix_data(
	std::vector<double> & matrix_dst,
	std::size_t & submatrix_rows,
	std::size_t & submatrix_columns,
	int src_proc_id,
	MPI_Status & status)
{
	MPI_Recv(&submatrix_rows, 1, MPI_UNSIGNED_LONG, src_proc_id, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&submatrix_columns, 1, MPI_UNSIGNED_LONG, src_proc_id, 1, MPI_COMM_WORLD, &status);

	matrix_dst.reserve(submatrix_rows * submatrix_columns);
	double * dst = new double[submatrix_rows * submatrix_columns];
	MPI_Recv(dst, submatrix_rows * submatrix_columns, MPI_DOUBLE, src_proc_id, 2, MPI_COMM_WORLD, &status);

	for (std::size_t i = 0; i < submatrix_rows * submatrix_columns; i++)
	{
		matrix_dst.push_back(dst[i]);
	}

	delete[] dst;
}

} // namespace

int main(int argc, char ** argv)
{
	int proc_id;
	int proc_num;
	math::matrix::Matrix A; // первая матрица
	math::matrix::Matrix B; // вторая матрица
	math::matrix::Matrix A_sub; // подматрицы, на которые разбивается матрица А
	math::matrix::Matrix C_sub; // результат умножения A_sub на B
	math::matrix::Matrix C; // результат умножения матрицы А на матрицу В

	// данные для заполнения подматриц
	std::size_t remain_size{0};
	std::size_t submatrix_rows_num{0};
	std::size_t submatrix_columns_num{0};
	std::vector<double> submatrix_matrix;

	// вектор, куда процессы отправляют результат вычислений
	std::vector<double> result_matrix;
	double * result_matrix_arr;

	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

	if (proc_id == 0)
	{
		A = math::matrix::factory::MatrixFactory::create_1diag_matrix(1000);
	}

	B = math::matrix::factory::MatrixFactory::create_1diag_matrix(1000);

	MPI_Barrier(MPI_COMM_WORLD);

	if (proc_id == 0)
	{
		count_submatrix_size(A, submatrix_rows_num, remain_size, proc_num);
		submatrix_columns_num = A.columns();
		submatrix_matrix.reserve(submatrix_rows_num * submatrix_columns_num);

		for (std::size_t i = 0; i < submatrix_rows_num * submatrix_columns_num; i++)
		{
			submatrix_matrix.push_back(A.data()[i]);
		}

		A_sub = math::matrix::Matrix(submatrix_rows_num, submatrix_columns_num, submatrix_matrix);

		for (int i = 1; i < (proc_num - 1); i++)
		{
			send_submatrix_data(A.data() + (i * submatrix_rows_num * A.columns()), submatrix_rows_num, A.columns(), i);
		}

		send_submatrix_data(A.data() + (A.rows() - remain_size) * A.columns(), remain_size, A.columns(), proc_num - 1);
	}
	else
	{
		recv_submatrix_data(submatrix_matrix, submatrix_rows_num, submatrix_columns_num, 0, status);
		A_sub = math::matrix::Matrix(submatrix_rows_num, submatrix_columns_num, submatrix_matrix);
	}

	C_sub = A_sub * B;

	if (proc_id == 0)
	{
		C = math::matrix::Matrix(A.rows(), B.columns());
		double * C_data = C.data();
		double * C_sub_data = C_sub.data();

		for (std::size_t i = 0; i < C_sub.rows() * C_sub.columns(); i++)
		{
			C_data[i] = C_sub_data[i];
		}

		for (int i = 1; i < (proc_num - 1); i++)
		{
			MPI_Recv(
				C_data + i * (submatrix_rows_num * C.columns()),
				submatrix_rows_num * C.columns(),
				MPI_DOUBLE,
				i,
				3,
				MPI_COMM_WORLD,
				&status);
		}

		MPI_Recv(
			C_data + (C.rows() - remain_size) * C.columns(),
			remain_size * C.columns(),
			MPI_DOUBLE,
			proc_num - 1,
			3,
			MPI_COMM_WORLD,
			&status);

		std::ofstream output("");
		output << C;
	}
	else
	{
		MPI_Send(C_sub.data(), C_sub.rows() * C_sub.columns(), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
	}

	//MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
}
