// who needs 1D-  Diddys version

#include <array>
#include <vector>
#include <iostream>
//#include <glm/glm.hpp>
#include <cmath>
#include <fstream>

//using Vec = glm::vec3;



std::vector<float> linspace(float start, float end, int n) {
    std::vector<float> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = start + (end - start)*i/(n-1.0f);
    }
    return result;
}




int main() {
    float L = 1.0f;
    float sigma = 0.1f;
    float c = 1.0f;
    float t_final = 1.0f;

    /// nx,ny stuff:
    float nx = 100;
    float ny = 100;
    float dx = L / (nx - 1.0f);
    float dy = L / (ny - 1.0f);

    float r = 0.5f;
    float rx = r;
    float ry = r;
    
    float dt = r*(dx/c);

    int nt = static_cast<int>(t_final / dt) + 1;

    std::vector<float> x = linspace(0.0f, L, static_cast<int>(nx));
    std::vector<float> y = linspace(0.0f, L, static_cast<int>(ny));
    std::vector<float> t = linspace(0.0f, t_final, nt);


    //Meshgrid  : 2v vectros
    // initalizatoin"
    std::vector<std::vector<float>> X(static_cast<int>(nx), std::vector<float>(static_cast<int>(ny)));
    std::vector<std::vector<float>> Y(static_cast<int>(nx), std::vector<float>(static_cast<int>(ny)));
    // filling:
    for (int i = 0; i < nx; ++i) {
        for (int k = 0; k < ny; ++k) {
            X[i][k] = x[i];
            Y[i][k] = y[k];
        }
    }

    // u array: 3D vector (time, x, y):
    std::vector<std::vector<std::vector<float>>> u(
        nt, 
        std::vector<std::vector<float>>(static_cast<int>(nx), 
                                        std::vector<float>(static_cast<int>(ny), 0.0f))
    );

    // iniial condition: Gaussian
    for (int i = 0; i < nx; ++i) {
        for(int k = 0; k < ny; ++k) {
            float dx_term = (X[i][k] - L/2.0f) * (X[i][k] - L/2.0f);
            float dy_term = (Y[i][k] - L/2.0f) * (Y[i][k] - L/2.0f);
            u[0][i][k] = std::exp(-(dx_term + dy_term)/(2.0f*sigma*sigma));
        }
    }

    // inital conditions of 0s:
    for (int i=0; i<nx; ++i){
        u[0][i][0] = 0.0f;
        u[0][i][ny-1] = 0.0f;
    }
    for (int k=0; k<nx; ++k){
        u[0][0][k] = 0.0f;
        u[0][nx-1][k] = 0.0f;
    }

    // first iteration
    for (int i=1; i < nx-1; ++i) {
        for (int k=1; k<ny-1; ++k) {
            u[1][i][k] = u[0][i][k] + (rx*rx)/2.0f * (u[0][i+1][k] - 2.0f*u[0][i][k] + u[0][i-1][k]) + (ry*ry)/2.0f * (u[0][i][k+1] - 2.0f*u[0][i][k] + u[0][i][k-1]);
        }
    }


    // inital conditions of 0s:
    for (int i=0; i<nx; ++i){
        u[1][i][0] = 0.0f;
        u[1][i][ny-1] = 0.0f;
    }
    for (int k=0; k<nx; ++k){
        u[1][0][k] = 0.0f;
        u[1][nx-1][k] = 0.0f;
    }

    // main loop (j = time, i = x, k = y):
    for (int j=1; j<nt-1; ++j) {
        for (int i=1; i < nx-1; ++i) {
            for (int k=1; k<ny-1; ++k) {
                u[j+1][i][k] = 2*u[j][i][k] - u[j-1][i][k] + (rx*rx)*(u[j][i+1][k] - 2*u[j][i][k] + u[j][i-1][k]) + (ry*ry)*(u[j][i][k+1] - 2*u[j][i][k] + u[j][i][k-1]);
            }
        }


        for (int i=0; i<nx; ++i){
            u[j+1][i][0] = 0.0f;
            u[j+1][i][ny-1] = 0.0f;
        }
        for (int k=0; k<nx; ++k){
            u[j+1][0][k] = 0.0f;
            u[j+1][nx-1][k] = 0.0f;
        }
    }

    
// Output initial condition
std::ofstream file_t0("wave_t0.csv");
for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < ny; ++k) {
        file_t0 << u[0][i][k];
        if (k < ny-1) file_t0 << ",";
    }
    file_t0 << "\n";
}
file_t0.close();

// Output final time
std::ofstream file_final("wave_final.csv");
for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < ny; ++k) {
        file_final << u[nt-1][i][k];
        if (k < ny-1) file_final << ",";
    }
    file_final << "\n";
}
file_final.close();

std::cout << "Saved wave_t0.csv and wave_final.csv" << std::endl;

    return 0;

}