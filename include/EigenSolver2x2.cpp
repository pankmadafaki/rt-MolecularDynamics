// EigenSolver2x2.cpp
#include "EigenSolver2x2.h"

namespace EigenSolver2x2 {

std::array<std::complex<float>, 2> eigenvalues(const std::array<std::array<float, 2>, 2>& matrix) {
    float a = matrix[0][0];
    float b = matrix[0][1];
    float c = matrix[1][0];
    float d = matrix[1][1];

    float trace = a + d;
    float determinant = a * d - b * c;

    float discriminant = trace * trace - 4 * determinant;

    std::complex<float> sqrt_discriminant;
    if (discriminant >= 0) {
        sqrt_discriminant = std::sqrt(discriminant);
    } else {
        sqrt_discriminant = std::sqrt(-discriminant) * std::complex<float>(0.0f, 1.0f);
    }

    std::array<std::complex<float>, 2> eigenValues;
    eigenValues[0] = (std::complex<float>(trace, 0.0f) + sqrt_discriminant) / 2.0f;
    eigenValues[1] = (std::complex<float>(trace, 0.0f) - sqrt_discriminant) / 2.0f;

    return eigenValues;
}

std::array<std::complex<float>, 2> eigenvector(const std::array<std::array<float, 2>, 2>& matrix, const std::complex<float>& eigenvalue) {
    std::complex<float> a = matrix[0][0];
    std::complex<float> b = matrix[0][1];
    std::complex<float> c = matrix[1][0];
    std::complex<float> d = matrix[1][1];

    std::complex<float> v1 = b;
    std::complex<float> v2 = eigenvalue - a;

    // If v1 and v2 are both close to zero, try the other column
    if (std::abs(v1) < 1e-6 && std::abs(v2) < 1e-6) {
        v1 = eigenvalue - d;
        v2 = c;
        if (std::abs(v1) < 1e-6 && std::abs(v2) < 1e-6) {
            throw std::runtime_error("No non-trivial eigenvector found for the given eigenvalue.");
        }
    }

    std::complex<float> norm = std::sqrt(std::norm(v1) + std::norm(v2));
    return {v1 / norm, v2 / norm};
}

} // namespace EigenSolver2x2
