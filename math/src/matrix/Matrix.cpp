#include <stdexcept>
#include <utility>

#include <matrix/Matrix.h>

namespace math
{
namespace matrix
{

Matrix::Matrix() = default;

Matrix::Matrix(const std::size_t rows, const std::size_t columns)
	: rows_{rows}
	, columns_{columns}
	, matrix_size_{rows * columns}
{
	memalloc();

	for (std::size_t i = 0; i < rows_; i++)
	{
		row_indexes_[i] = i * columns_;
	}
}

Matrix::Matrix(const std::size_t rows, const std::size_t columns, std::vector<double> matrix_src)
	: rows_{rows}
	, columns_{columns}
	, matrix_size_{rows * columns}
{
	if (matrix_size_ != matrix_src.size())
	{
		throw std::runtime_error{"Exception"};
	}

	memalloc();

	fill_indexes();

	for (std::size_t i = 0; i < matrix_src.size(); i++)
	{
		matrix_[i] = matrix_src[i];
	}
}

Matrix::Matrix(
	double * matrix_arr,
	std::size_t * row_indexes,
	const std::size_t rows,
	const std::size_t columns,
	const std::size_t matrix_size)
	: matrix_{matrix_arr}
	, row_indexes_{row_indexes}
	, rows_{rows}
	, columns_{columns}
	, matrix_size_{matrix_size}
{
}

Matrix::Matrix(const Matrix & matrix)
	: rows_{matrix.rows_}
	, columns_{matrix.columns_}
	, matrix_size_{matrix.matrix_size_}
{
	memalloc();
	memcopy(matrix.matrix_, matrix.row_indexes_);
}

Matrix::Matrix(Matrix && matrix)
{
	memmove(std::move(matrix));
}

Matrix::~Matrix()
{
	memclear();
}

Matrix & Matrix::operator=(const Matrix & matrix)
{
	memclear();

	rows_ = matrix.rows_;
	columns_ = matrix.columns_;
	matrix_size_ = matrix.matrix_size_;

	memalloc();
	memcopy(matrix.matrix_, matrix.row_indexes_);

	return *this;
}

std::ifstream & operator>>(std::ifstream & input, Matrix & matrix)
{
	matrix.memclear();

	input >> matrix.rows_;
	input >> matrix.columns_;
	matrix.matrix_size_ = matrix.rows_ * matrix.columns_;

	matrix.memalloc();
	matrix.fill_indexes();

	for (std::size_t i = 0; i < matrix.matrix_size_; i++)
	{
		input >> matrix.matrix_[i];
	}

	return input;
}

std::ofstream & operator<<(std::ofstream & output, Matrix & matrix)
{
	for (std::size_t i = 0; i < matrix.rows_; i++)
	{
		for (std::size_t j = 0; j < matrix.columns_; j++)
		{
			output << matrix.matrix_[matrix.row_indexes_[i] + j] << " ";
		}

		output << std::endl;
	}

	return output;
}

Matrix & Matrix::operator=(Matrix && matrix)
{
	memclear();
	memmove(std::move(matrix));

	return *this;
}

Matrix operator*(const Matrix & m1, const Matrix & m2)
{
	if (m1.columns_ != m2.rows_)
	{
		throw std::runtime_error{"Exception"};
	}

	std::size_t result_rows = m1.rows_;
	std::size_t result_columns = m2.columns_;
	std::size_t result_matrix_size = result_rows * result_columns;

	double * result_matrix = new double[result_matrix_size];
	std::size_t * result_rows_indexes = new std::size_t[result_rows];

	for (std::size_t i = 0; i < result_rows; i++)
	{
		result_rows_indexes[i] = i * result_columns;
	}

	for (std::size_t i = 0; i < result_rows; i++)
	{
		for (std::size_t j = 0; j < result_columns; j++)
		{
			std::size_t index = result_rows_indexes[i] + j;
			result_matrix[index] = 0;
			for (std::size_t k = 0; k < m2.rows_; k++)
			{
				result_matrix[index] += m1.matrix_[m1.row_indexes_[i] + k] * m2.matrix_[m2.row_indexes_[k] + j];
			}
		}
	}

	return Matrix(result_matrix, result_rows_indexes, result_rows, result_columns, result_matrix_size);
}

double Matrix::operator()(std::size_t row_index, std::size_t column_index)
{
	if (row_index > rows_ || column_index > columns_)
	{
		throw std::runtime_error{"Incorrect index"};
	}

	return matrix_[row_indexes_[row_index] + column_index];
}

double * Matrix::operator[](std::size_t column_index)
{
	if (column_index > columns_)
	{
		throw std::runtime_error{"Incorrect index"};
	}

	double * result = new double[rows_];

	for (std::size_t i = 0; i < rows_; i++)
	{
		result[i] = matrix_[row_indexes_[i] + column_index];
	}

	return result;
}

Matrix Matrix::submatrix(std::size_t start_row_index, std::size_t end_row_index) const
{
	if ((start_row_index >= end_row_index) || start_row_index > rows_ || end_row_index > rows_)
	{
		throw std::runtime_error{"Incorrect params"};
	}

	std::size_t result_rows{end_row_index - start_row_index};
	std::size_t result_columns{columns_};
	std::size_t result_matrix_size{result_rows * result_columns};

	double * result_matrix = new double[matrix_size_];
	std::size_t * result_rows_indexes = new std::size_t[result_rows];

	std::size_t iC{0};
	for (std::size_t i = start_row_index; i < end_row_index; i++)
	{
		result_matrix[iC++] = matrix_[i];
	}

	for (std::size_t i = 0; i < result_rows; i++)
	{
		result_rows_indexes[i] = i * result_columns;
	}

	return Matrix(result_matrix, result_rows_indexes, result_rows, result_columns, result_matrix_size);
}

void Matrix::memalloc()
{
	matrix_ = new double[matrix_size_];
	row_indexes_ = new std::size_t[rows_];
}

void Matrix::memmove(Matrix && matrix)
{
	rows_ = matrix.rows_;
	columns_ = matrix.columns_;
	matrix_size_ = matrix.matrix_size_;
	matrix_ = matrix.matrix_;
	row_indexes_ = matrix.row_indexes_;

	matrix.matrix_ = nullptr;
	matrix.row_indexes_ = nullptr;
}

void Matrix::memclear()
{
	if (matrix_)
	{
		delete[] matrix_;
		delete[] row_indexes_;
	}
}

void Matrix::memcopy(double * src, std::size_t * indexes_src)
{
	for (std::size_t i = 0; i < matrix_size_; i++)
	{
		matrix_[i] = src[i];
	}

	for (std::size_t i = 0; i < rows_; i++)
	{
		row_indexes_[i] = indexes_src[i];
	}
}

void Matrix::fill_indexes()
{
	for (std::size_t i = 0; i < rows_; i++)
	{
		row_indexes_[i] = i * columns_;
	}
}

} // namespace matrix
} // namespace math
