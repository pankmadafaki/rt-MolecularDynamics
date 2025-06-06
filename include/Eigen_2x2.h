#ifndef EIGEN_2X2_H
#define EIGEN_2X2_H

#include <array>
#include <complex>
#include <cmath>
#include <tuple>
#include <limits>
#include <algorithm>
#include <vector>

namespace LinearAlgebra {

/**
 * @brief Calculates the eigenvalues and eigenvectors of a 2x2 real matrix and returns the eigenvector
 * associated with the eigenvalue having the largest magnitude.
 *
 * @param matrix A 2x2 matrix represented as a std::array of std::arrays of floats.
 * @return A std::tuple containing:
 * - A std::array of two std::complex<float> representing the eigenvalues.
 * - A std::array of two std::array<std::complex<float>, 2> representing the corresponding eigenvectors (as columns).
 * Eigenvectors are normalized if the eigenvalues are real. If eigenvalues are complex conjugates, the eigenvectors
 * will also be complex conjugates.
 * - A std::array<std::complex<float>, 2> representing the eigenvector associated with the eigenvalue having the
 * largest magnitude. If there's a tie in magnitude, the eigenvector corresponding to the first encountered
 * largest eigenvalue is returned.
 */
std::tuple<std::array<std::complex<float>, 2>, std::array<std::array<std::complex<float>, 2>, 2>, std::array<std::complex<float>, 2>, std::complex<float>>
calculateEigen(const std::array<std::array<float, 2>, 2>& matrix) {
    float a = matrix[0][0];
    float b = matrix[0][1];
    float c = matrix[1][0];
    float d = matrix[1][1];

    // Calculate the discriminant of the characteristic polynomial
    float discriminant = (a + d) * (a + d) - 4 * (a * d - b * c);

    std::array<std::complex<float>, 2> eigenvalues;
    std::array<std::array<std::complex<float>, 2>, 2> eigenvectors;

    if (discriminant >= 0) {
        // Real eigenvalues
        float sqrt_discriminant = std::sqrt(discriminant);
        eigenvalues[0] = 0.5f * (a + d + sqrt_discriminant);
        eigenvalues[1] = 0.5f * (a + d - sqrt_discriminant);

        // Calculate eigenvectors
        for (size_t i = 0; i < 2; ++i) {
            std::complex<float> lambda = eigenvalues[i];
            std::complex<float> v1 = b;
            std::complex<float> v2 = lambda - a;

            float norm = std::sqrt(std::norm(v1) + std::norm(v2));
            if (norm > std::numeric_limits<float>::epsilon()) {
                eigenvectors[i][0] = v1 / norm;
                eigenvectors[i][1] = v2 / norm;
            } else {
                // Handle the case where v1 and v2 are both close to zero (shouldn't happen for distinct eigenvalues)
                eigenvectors[i][0] = 1.0f;
                eigenvectors[i][1] = 0.0f;
            }
        }
    } else {
        // Complex conjugate eigenvalues
        float real_part = 0.5f * (a + d);
        float imaginary_part = 0.5f * std::sqrt(-discriminant);
        eigenvalues[0] = std::complex<float>(real_part, imaginary_part);
        eigenvalues[1] = std::complex<float>(real_part, -imaginary_part);

        // Calculate eigenvectors (will also be complex conjugates)
        std::complex<float> lambda0 = eigenvalues[0];
        std::complex<float> v1_0 = b;
        std::complex<float> v2_0 = lambda0 - a;
        eigenvectors[0][0] = v1_0;
        eigenvectors[0][1] = v2_0;

        std::complex<float> lambda1 = eigenvalues[1];
        std::complex<float> v1_1 = b;
        std::complex<float> v2_1 = lambda1 - a;
        eigenvectors[1][0] = v1_1;
        eigenvectors[1][1] = v2_1;
    }

    // Find the eigenvector associated with the largest eigenvalue (in magnitude)
    size_t largest_eigenvalue_index = 0;
    if (std::abs(eigenvalues[1]) > std::abs(eigenvalues[0])) {
        largest_eigenvalue_index = 1;
    }

    std::array<std::complex<float>, 2> largest_eigenvector = eigenvectors[largest_eigenvalue_index];
    std::complex<float> largest_eigenvalue = eigenvalues[largest_eigenvalue_index];

    return std::make_tuple(eigenvalues, eigenvectors, largest_eigenvector, largest_eigenvalue);
}

} // namespace LinearAlgebra

#endif // EIGEN_2X2_H
