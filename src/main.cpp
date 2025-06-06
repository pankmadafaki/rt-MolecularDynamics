#include "../include/raylib-cpp.hpp"
#include <cmath>
#include <math.h>
// #include <cstddef>
// #include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <random>
#include <iomanip>
#include "../include/Eigen_2x2.h"
#include <sciplot/Plot2D.hpp>
#include <sciplot/sciplot.hpp>
#include <vector>
#include "../include/MolDym.cpp"
#define RAYGUI_IMPLEMENTATION
#include "../include/raygui.h"
#include "../include/raylib.h"

#define SCREENWIDTH 800
#define SCREENHEIGHT 800


int main() {
    const int n_mol = 25*25; //841
    const float radius = 12.0f;
    const float box_size = float(SCREENWIDTH);
    float sigma = 1.1; // 1.1
    float rc0 = 50.0; // 60.0
    float epsilon = 0.3; // 0.3
    float epsilon0 = 20.0;
    float T = 0.0;
    int counter = 0;
    int numbins = 300;
    double Nf = 2*n_mol - 2 -1;
    std::vector<double> ek_list;
    std::vector<double> ep_list;
    std::vector<double> g(numbins, 0.0);
    std::vector<double> msd_list;
    // Raylib stuff for font and appearance
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(SCREENWIDTH + 400, SCREENHEIGHT, "Velocity-Verlet");
    GuiLoadStyle("include/dark.rgs");
    Font fontTtf = LoadFontEx("include/PixelOperator.ttf", 32, 0, 250); //    Load
    // Raylib stuff for font and appearance
    SetTargetFPS(120);
    float dt = 0.15;
    // Create some circles with initial velocities
    std::vector<Molecule> molecules = initVelocities(0.1, n_mol);
    std::vector<Molecule> initial_molecules = molecules;
    while (!WindowShouldClose()) {
        
        // Update circle positions
        auto [ep, ek, g0] = velocityVerlet(molecules, n_mol, sigma, rc0, epsilon, dt, 0.001, T, numbins, epsilon0);
        double msd = meanSquareDisplacement(molecules, initial_molecules, n_mol);
        ep_list.push_back(ep);
        ek_list.push_back(ek/Nf);
        msd_list.push_back(msd);
        counter ++;

        for (int i = 0; i < numbins; i++) {
          g[i] =+ g0[i];
        }


        // Drawing
        BeginDrawing();
        ClearBackground(RAYWHITE);
        
        
        for (auto &c : molecules) {
            Vector2 basePos = {c.position.x, c.position.y};

            // Define offset positions to simulate periodic boundary copies
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    Vector2 offsetPos = {
                        basePos.x + dx * SCREENWIDTH,
                        basePos.y + dy * SCREENHEIGHT
                    };

                    // Only draw if the copy is within the visible screen area (optional)
                    if (offsetPos.x + radius >= 0 && offsetPos.x - radius <= SCREENWIDTH &&
                        offsetPos.y + radius >= 0 && offsetPos.y - radius <= SCREENHEIGHT) {
                        DrawCircleV(offsetPos, radius, DARKGRAY);
                    }
                }
            }
        }
        DrawRectangle(SCREENWIDTH, 0, SCREENWIDTH - 400, SCREENHEIGHT,
                    GRAY);
        
        GuiSliderBar((Rectangle){SCREENWIDTH + 90, 40, 220, 20}, "dt", NULL, &dt, 0.02, 0.3);
        
        GuiSliderBar((Rectangle){SCREENWIDTH + 90, 70, 220, 20}, "sigma", NULL, &sigma, 0.1, 1.5);
        GuiSliderBar((Rectangle){SCREENWIDTH + 90, 90, 220, 20}, "rc_0", NULL, &rc0, 0.0, 600.0);
        GuiSliderBar((Rectangle){SCREENWIDTH + 90, 110, 220, 20}, "epsilon", NULL, &epsilon, 0.0, 3.0);
        GuiSliderBar((Rectangle){SCREENWIDTH + 90, 130, 220, 20}, "Temperature", NULL, &T, 0.0, 600.0);

        DrawTextEx(fontTtf, TextFormat("%.2f", dt), {SCREENWIDTH + 90 + 220, 40}, 20, 0.2, WHITE);
        DrawTextEx(fontTtf, TextFormat("%.2f", sigma), {SCREENWIDTH + 90 + 220, 70}, 20, 0.2, WHITE);
        DrawTextEx(fontTtf, TextFormat("%.2f", rc0), {SCREENWIDTH + 90 + 220, 90}, 20, 0.2, WHITE);
        DrawTextEx(fontTtf, TextFormat("%.2f", epsilon), {SCREENWIDTH + 90 + 220, 110}, 20, 0.2, WHITE);
        DrawTextEx(fontTtf, TextFormat("%.2f", T), {SCREENWIDTH + 90 + 220, 130}, 20, 0.2, WHITE);

        DrawFPS(0,0);
        EndDrawing();
        // Drawing - END
        if (counter == 10000) {
          break;
        }
    }
    CloseWindow(); // Cleanup
    for (auto &b: g) {
      b = b/counter;      
    }
    saveVectorToFile(ep_list, "e_p_600K.txt");
    saveVectorToFile(g, "g_gas.txt");
    saveVectorToFile(ek_list, "e_k_600K.txt");
    saveVectorToFile(msd_list, "/home/lr0n/documents/Uni/CS/Ex4/msd/msd.txt");
    return 0;
}

