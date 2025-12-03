#include <matrix/Matrix.h>
#include <matrix/factory/MatrixFactory.h>

int main()
{
    std::ofstream output("");
    
    auto A = math::matrix::factory::MatrixFactory::create_1diag_matrix(1000);
    auto B = math::matrix::factory::MatrixFactory::create_1diag_matrix(1000);
    auto C = A * B;

    output << C;
}