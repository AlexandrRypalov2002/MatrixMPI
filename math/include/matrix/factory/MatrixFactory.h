#ifndef MATH_MATRIX_FACTORY_MATRIXFACTORY_H
#define MATH_MATRIX_FACTORY_MATRIXFACTORY_H

#include <cstdint>

namespace math
{
namespace matrix
{
class Matrix;
namespace factory
{

class MatrixFactory
{
public:
    /**
     * Создание квадратной матрицы с единицами на главной диагонали
     * @param size размер матрицы (будет создана размера size * size)
     */
    static Matrix create_1diag_matrix(std::size_t size);
};

} // namespace factory
} // namespace matrix
} // namespace math

#endif
