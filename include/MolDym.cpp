#include "MolDym.hpp"

#define SCREENWIDTH 800
#define SCREENHEIGHT 800



void saveVectorToFile(const std::vector<double>& data, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    for (double value : data) {
        outFile << value << '\n';  // Write one value per line
    }

    outFile.close();
    std::cout << "Vector saved to " << filename << "\n";
}

Vector2d force_WF(const Vector2d& r_i, const Vector2d& r_j, double sigma, double rc, double epsilon) {
    Vector2d r_ij = r_i - r_j;
    double rij = r_ij.norm();
    if (abs(rij)>1000.0) {
      std::cout << rij << " abs rij \n";
      
    }
    if (rij >= rc || rij == 0.0) {
        //std::cout << "no calc\n";
        return Vector2d(0.0, 0.0);
    }

    // Compute intermediate quantities
    double rij2 = rij * rij;
    double rij3 = rij2 * rij;

    double sigma2 = sigma * sigma;
    double rc2 = rc * rc;

    double A = (sigma2 / rij2) - 1.0;
    double C = (rc2 / rij2) - 1.0;
    double B = C * C;

    double dA = -2.0 * sigma2 / rij3;
    double dB = -4.0 * rc2 / rij3 * C;
    

    // Derivative of the potential with respect to rij
    double dU_dr = epsilon * (dA * B + A * dB);
    //std::cout << dU_dr << " inside force " << rij << " " << rij2 << "\n";
      

    
    // Compute force vector: -dU/dr * unit vector
    return r_ij.normalized() * (-dU_dr);
}

double potential_WF(const Vector2d& r_i, const Vector2d& r_j, double sigma, double rc, double epsilon) {
    Vector2d r_ij = r_i - r_j;
    double rij = r_ij.norm();

    if (rij >= rc || rij == 0.0) {
        return 0.0;
    }

    double sigma2 = sigma * sigma;
    double rc2 = rc * rc;
    double rij2 = rij * rij;

    double term1 = (sigma2 / rij2) - 1.0;
    double term2 = (rc2 / rij2) - 1.0;

    double potential = epsilon * term1 * term2 * term2;
    return potential;
}



std::vector<Molecule> initVelocities(float temperature, int N) {
    std::vector<Molecule> particles;
    particles.reserve(N);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float> velocity_dist(0.0f, temperature);

    // Determine how many particles per side we need (assuming a square grid)
    int particlesPerSide = std::ceil(std::sqrt(N));

    // Compute spacing a such that we fit particlesPerSide particles with spacing a and buffer a/2 on each side
    float a = SCREENWIDTH / particlesPerSide;

    // Bounds check: make sure a is reasonable
    if (a <= 0.0f) {
        throw std::runtime_error("Computed spacing 'a' is non-positive.");
    }

    // Start placing particles
    int count = 0;
    for (int row = 0; row < particlesPerSide && count < N; ++row) {
        for (int col = 0; col < particlesPerSide && count < N; ++col) {
            float x = a * col + a / 2.0f;
            float y = a * row + a / 2.0f;

            particles.push_back({
                {x, y},
                {velocity_dist(gen), velocity_dist(gen)},
                {0.0f, 0.0f}
            });

            ++count;
        }
    }

    return particles;
}
double meanSquareDisplacement(std::vector<Molecule>& molecules, std::vector<Molecule>& initial_molecules, int n_mol) {
  double msd = 0.0;
  for (int i = 0; i < n_mol; i++) {
    Vector2d displ = molecules[i].position - initial_molecules[i].position;
    msd =+ displ.x * displ.x + displ.y * displ.y;
  }
  return msd/n_mol;
}
std::tuple<double, double, std::vector<double>> velocityVerlet(std::vector<Molecule>& molecules, int n_mol, double sigma, double rc0, double epsilon, float dt, float delta, float T, int numBins, float epsilon0) {
  std::vector<Molecule> new_molecules = molecules;
  double U_tot = 0.0;
  double Ek_tot = 0.0;
  double m = 1.0e-08;
  double Nf = 2*n_mol-2-1;
  std::random_device rd;  
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine
  std::uniform_real_distribution<> dis(-delta, delta);
  const double k_B = 1.380649e-23; // Boltzmann constant 
  double dr = SCREENWIDTH / (2.0 * numBins);  // max radius = boxLength / 2
  std::vector<double> g(numBins, 0.0);
  std::uniform_real_distribution<> dist_uni(0.0, 1.0); // Uniform distribution between 0.0 and 1.0
  
  const float box_size = float(SCREENWIDTH);
  float z = exp(-dis(gen));
  for (int i = 0; i < n_mol; i++) {
          
          Vector2d force_ij;
          Ek_tot =+ 0.5*m*(std::pow(molecules[i].velocity.x, 2) +  std::pow(molecules[i].velocity.y, 2));
          for (int j = 0; j < n_mol; j++) {
              if (i != j) {
                Vector2d r_ij = molecules[i].position - molecules[j].position;
                U_tot += potential_WF(molecules[i].position, molecules[j].position, sigma, rc0, epsilon);

                if (r_ij.x > 0.5 * box_size) r_ij.x -= box_size;
                else if (r_ij.x < -0.5 * box_size) r_ij.x += box_size;

                if (r_ij.y > 0.5 * box_size) r_ij.y -= box_size;
                else if (r_ij.y < -0.5 * box_size) r_ij.y += box_size;
                double rij = r_ij.norm();
                int bin = static_cast<int>(rij / dr);
                if (bin < numBins) {
                  if (abs(rij - bin * dr) < epsilon0) {
                    double delta = 1.0;                    
                  } else {
                    double delta = 0.0;
                  }
                  double shellArea = 2 * M_PI * (bin * dr);
  
                  g[bin] += delta / shellArea;
                }
                force_ij = force_ij + force_WF(Vector2d(0, 0), r_ij, sigma, rc0, epsilon);
              } 
          }
          double normFactor = (2.0 * SCREENWIDTH * SCREENWIDTH) / (static_cast<double>(n_mol) * n_mol);
          
          
          Vector2d new_pos = molecules[i].position + molecules[i].velocity * dt + molecules[i].acceleration * (dt * dt * 0.5);
          
          if (new_pos.x >= box_size) new_pos.x -= box_size;
          else if (new_pos.x < 0) new_pos.x += box_size;

          if (new_pos.y >= box_size) new_pos.y -= box_size;
          else if (new_pos.y < 0) new_pos.y += box_size;

          Vector2d new_acc = force_ij;
          Vector2d new_vel = molecules[i].velocity + (molecules[i].acceleration + new_acc) * (dt * 0.5);
          // Rescale initVelocities

          new_molecules[i].position = new_pos;
          new_molecules[i].velocity = new_vel;
          new_molecules[i].acceleration = new_acc;
        }
        double p = std::pow(z,2*(n_mol -1))*exp(-(Ek_tot*(z*z-1))/(T));
        if (p > 1.0) {
          p = 1.0;          
        }
        if (dist_uni(gen) < p) {
          molecules = new_molecules;
          for (auto& molecule : molecules) {
            molecule.velocity = molecule.velocity*z ;
          }
        } else {
          molecules = new_molecules;
        }
        return std::make_tuple(U_tot, Ek_tot/(k_B*Nf), g);
}

