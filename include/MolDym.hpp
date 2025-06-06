#ifndef MOLDYM_HPP
#define MOLDYM_HPP

#include <cmath>
#include <math.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <iomanip>
#include <vector>


struct Vector2d {
    double x, y;

    Vector2d() : x(0), y(0) {}
    Vector2d(double x, double y) : x(x), y(y) {}

    Vector2d operator-(const Vector2d& other) const {
        return Vector2d(x - other.x, y - other.y);
    }

    Vector2d operator+(const Vector2d& other) const {
        return Vector2d(x + other.x, y + other.y);
    }


    Vector2d operator*(double scalar) const {
        return Vector2d(x * scalar, y * scalar);
    }
    
    Vector2d operator/(double scalar) const {
        return Vector2d(x / scalar, y / scalar);
    }

    double norm() const {
        return std::sqrt(x * x + y * y);
    }

    Vector2d normalized() const {
        double n = norm();
        if (n == 0.0) return Vector2d(0, 0); // Avoid division by zero
        return Vector2d(x / n, y / n);
    }
};
struct Molecule {
    Vector2d position;
    Vector2d velocity;
    Vector2d acceleration;
};

void saveVectorToFile(const std::vector<double>& data, const std::string& filename);

Vector2d force_WF(const Vector2d& r_i, const Vector2d& r_j, double sigma, double rc, double epsilon); 

double potential_WF(const Vector2d& r_i, const Vector2d& r_j, double sigma, double rc, double epsilon, float dt);

std::vector<Molecule> initVelocities(float temperature, int N);

double meanSquareDisplacement(std::vector<Molecule>& molecules, std::vector<Molecule>& initial_molecules, int n_mol);

std::tuple<double, double, std::vector<double>> velocityVerlet(std::vector<Molecule>, int n_mol, double sigma, double rc0, double epsilon, float delta, float T, int numBins, float epsilon0);

#endif

