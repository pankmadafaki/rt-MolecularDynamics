// EigenSolver2x2.h
#ifndef EIGEN_SOLVER_2X2_H
#define EIGEN_SOLVER_2X2_H

#include <array>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace EigenSolver2x2 {

/**
 * @brief Calculates the eigenvalues of a 2x2 matrix.
 *
 * Given a 2x2 matrix represented by a std::array of std::arrays of floats,
 * this function calculates its two eigenvalues (which may be complex).
 *
 * @param matrix The 2x2 matrix as a std::array<std::array<float, 2>, 2>.
 * @return A std::array of two std::complex<float> representing the eigenvalues.
 */
std::array<std::complex<float>, 2> eigenvalues(const std::array<std::array<float, 2>, 2>& matrix);

/**
 * @brief Calculates the eigenvectors corresponding to a given eigenvalue of a 2x2 matrix.
 *
 * Given a 2x2 matrix and one of its eigenvalues, this function calculates
 * the corresponding eigenvector. The eigenvector is returned as a normalized
 * std::array of two std::complex<float>.
 *
 * @param matrix The 2x2 matrix as a std::array<std::array<float, 2>, 2>.
 * @param eigenvalue The eigenvalue for which to find the eigenvector.
 * @return A std::array of two std::complex<float> representing the normalized eigenvector.
 * @throws std::runtime_error if no non-trivial eigenvector is found.
 */
std::array<std::complex<float>, 2> eigenvector(const std::array<std::array<float, 2>, 2>& matrix, const std::complex<float>& eigenvalue);

} // namespace EigenSolver2x2

#endif // EIGEN_SOLVER_2X2_H
