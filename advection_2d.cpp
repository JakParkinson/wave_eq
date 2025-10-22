#include <array>
#include <vector>
#include <iostream>
//#include <glm/glm.hpp>
#include <cmath>
#include <fstream>



std::vector<float> linspace(float start, float end, int n) {
    std::vector<float> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = start + (end - start)*i/(n-1.0f);
    }
    return result;
}


int main (){
    // parameters:
    float L = 4.0f;
    float cx = 1.5f;
    float cy = 3.5f;
    // IC:
    float sigma = 0.5f;
    float amplitude = 5.0f;
    float t_final = 3.0f;
    // discretization:
    int nx = 100;
    int ny = 100;
    float dx = L/ static_cast<float>(nx - 1);
    float dy = L/ static_cast<float>(ny - 1);
    float dt = std::min(dx/cx, dy/cy);
    int nt = static_cast<int>(t_final / dt) + 1;
    
    std::vector<float> x = linspace(0.0f, L, nx);
    std::vector<float> y = linspace(0.0f, L, ny);
    std::vector<float> t = linspace(0.0f, t_final, nt);


    std::vector<std::vector<float>> X(nx, std::vector<float>(ny));
    std::vector<std::vector<float>> Y(nx, std::vector<float>(ny));


    for (int i = 0; i < nx; ++i) {
        for (int k = 0; k < ny; ++k) {
            X[i][k] = x[i];
            Y[i][k] = y[k];
        }
    }

    std::vector<std::vector<std::vector<float>>> u(
        nt, 
        std::vector<std::vector<float>>(nx, 
                                        std::vector<float>(ny, 0.0f))
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
    for (int k=0; k<ny; ++k){
        u[0][0][k] = 0.0f;
        u[0][nx-1][k] = 0.0f;
    }

    // main loop (j=time,i=x,k=y)
    for (int j=0; j<nt-1; ++j) {
        // k arrays:
        std::vector<std::vector<float>> k1(nx, std::vector<float>(ny,0.0f));
        std::vector<std::vector<float>> k2(nx, std::vector<float>(ny,0.0f));
        std::vector<std::vector<float>> k3(nx, std::vector<float>(ny,0.0f));
        std::vector<std::vector<float>> k4(nx, std::vector<float>(ny,0.0f));

        // k1:
        for(int i=1; i<nx-1; ++i) {
            for (int k=1; k<ny-1; ++k) {
                k1[i][k] = -cx*(u[j][i][k] - u[j][i-1][k])/dx - cy*(u[j][i][k] - u[j][i][k-1])/dy;
            }
        }

        //k2s:
        std::vector<std::vector<float>> u_temp(nx, std::vector<float>(ny, 0.0f));
        for (int i=0;i<nx;++i){
            for (int k=0;k<ny;++k){
                u_temp[i][k]=u[j][i][k] + (dt/2.0f)*k1[i][k];
            }
        }
        // BC (dirichlet):
        for (int i =0; i < nx; ++i) {
            u_temp[i][0] = 0.0f;
            u_temp[i][ny-1] = 0.0f;
        }
        for (int k =0; k < ny; ++k) {
            u_temp[0][k] = 0.0f;
            u_temp[ny-1][k] = 0.0f;
        }
        for (int i = 1; i<nx-1; ++i){
            for (int k=1; k<ny-1; ++k){
                k2[i][k] = -cx*(u_temp[i][k] - u_temp[i-1][k])/dx - cy*(u_temp[i][k] - u_temp[i][k-1])/dy;
            }
        }

        //k3:
        for (int i=0;i<nx;++i){
            for (int k=0;k<ny;++k){
                u_temp[i][k]=u[j][i][k] + (dt/2.0f)*k2[i][k];
            }
        }
        // BC (dirichlet):
        for (int i =0; i < nx; ++i) {
            u_temp[i][0] = 0.0f;
            u_temp[i][ny-1] = 0.0f;
        }
        for (int k =0; k < ny; ++k) {
            u_temp[0][k] = 0.0f;
            u_temp[ny-1][k] = 0.0f;
        }
        for (int i = 1; i<nx-1; ++i){
            for (int k=1; k<ny-1; ++k){
                k3[i][k] = -cx*(u_temp[i][k] - u_temp[i-1][k])/dx - cy*(u_temp[i][k] - u_temp[i][k-1])/dy;
            }
        }

        //k4:
        for (int i=0;i<nx;++i){
            for (int k=0;k<ny;++k){
                u_temp[i][k]=u[j][i][k] + dt*k3[i][k];
            }
        }
        // BC (dirichlet):
        for (int i =0; i < nx; ++i) {
            u_temp[i][0] = 0.0f;
            u_temp[i][ny-1] = 0.0f;
        }
        for (int k =0; k < ny; ++k) {
            u_temp[0][k] = 0.0f;
            u_temp[ny-1][k] = 0.0f;
        }
        for (int i = 1; i<nx-1; ++i){
            for (int k=1; k<ny-1; ++k){
                k4[i][k] = -cx*(u_temp[i][k] - u_temp[i-1][k])/dx - cy*(u_temp[i][k] - u_temp[i][k-1])/dy;
            }
        }

        for (int i=0;i<nx;++i){
            for (int k=0;k<ny;++k){
                u[j+1][i][k] = u[j][i][k] + (dt/6.0f)*(k1[i][k]+2.0f*k2[i][k]+2.0f*k3[i][k]+k4[i][k]);
            }
        }

        for (int i = 0; i < nx; ++i) {
            u[j+1][i][0] = 0.0f;
            u[j+1][i][ny-1] = 0.0f;
        }
        for (int k = 0; k < ny; ++k) {
            u[j+1][0][k] = 0.0f;
            u[j+1][nx-1][k] = 0.0f;
        }

    }



    // Output initial condition
    std::ofstream file_t0("adv_t0.csv");
    for (int i = 0; i < nx; ++i) {
        for (int k = 0; k < ny; ++k) {
            file_t0 << u[0][i][k];
            if (k < ny-1) file_t0 << ",";
        }
        file_t0 << "\n";
    }
    file_t0.close();

    // Output final time
    std::ofstream file_final("adv_final.csv");
    for (int i = 0; i < nx; ++i) {
        for (int k = 0; k < ny; ++k) {
            file_final << u[nt-1][i][k];
            if (k < ny-1) file_final << ",";
        }
        file_final << "\n";
    }
    file_final.close();

    std::cout << "Saved adv_t0.csv and adv_final.csv" << std::endl;

    return 0;
    
}