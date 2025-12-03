#include <matrix/Matrix.h>
#include <matrix/factory/MatrixFactory.h>

#include <vector>

namespace math
{
namespace matrix
{
namespace factory
{

Matrix MatrixFactory::create_1diag_matrix(std::size_t size)
{
	std::vector<double> result_matrix(size * size, 0);

	for (std::size_t i = 0; i < size; i++)
	{
		result_matrix[i + i * size] = 1;
	}

	return Matrix(size, size, result_matrix);
}

} // namespace factory
} // namespace matrix
} // namespace math
