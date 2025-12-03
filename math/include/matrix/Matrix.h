#ifndef MATH_MATRIX_MATRIX_H
#define MATH_MATRIX_MATRIX_H

#include <cstdint>
#include <fstream>
#include <vector>

namespace math
{
namespace matrix
{

class Matrix
{
public:
	Matrix();

	Matrix(const std::size_t rows, const std::size_t columns);

	Matrix(const std::size_t rows, const std::size_t columns, std::vector<double> matrix_src);

	Matrix(const Matrix & other);

	Matrix(Matrix && other);

	~Matrix();

	Matrix & operator=(const Matrix & matrix);

	Matrix & operator=(Matrix && matrix);

	friend Matrix operator*(const Matrix & m1, const Matrix & m2);

	friend std::ifstream & operator>>(std::ifstream & input, Matrix & matrix);

	friend std::ofstream & operator<<(std::ofstream & output, Matrix & matrix_output);

	double operator()(std::size_t row_index, std::size_t column_index);

	double * operator[](std::size_t column_index);

	Matrix submatrix(std::size_t start_row_index, std::size_t end_row_index) const;

	inline double * data() noexcept
	{
		return matrix_;
	};

	inline std::size_t rows() const noexcept
	{
		return rows_;
	};

	inline std::size_t columns() const noexcept
	{
		return columns_;
	};

private:
	Matrix(
		double * matrix_arr,
		std::size_t * row_indexes,
		const std::size_t rows,
		const std::size_t columns,
		const std::size_t matrix_size);

	void memalloc();

	void memmove(Matrix && other);

	void memclear();

	void memcopy(double * matrix_src, std::size_t * indexes_src);

	void fill_indexes();

private:
	double * matrix_{nullptr};
	std::size_t * row_indexes_{nullptr};

	std::size_t rows_{0};
	std::size_t columns_{0};
	std::size_t matrix_size_{0};
};

} // namespace matrix
} // namespace math

#endif
