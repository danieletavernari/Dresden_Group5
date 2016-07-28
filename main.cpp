/* 
 * File:   main.cpp
 * Author: daniele
 *
 * Created on 5 dicembre 2015, 17.54
 * 
 * NOTE
 * libreria "DUNE" per risolvere pde
 * if (abs(expr1 –expr2) < 0.000001) ...  // Ok !
 * std::list could be better than std::vector. 
 * Boost vector should not be better (it actually depends on compilers-machines (try both)
 * preallocate vectors with v.reserve()
 * 
 */
#include <iostream>
#include <math.h>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <numeric>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>
#include <fftw3.h>

using namespace std;
using namespace boost::numeric::ublas;
const double pi = boost::math::constants::pi<double>();

const double Kap = 1;
const double Kpa = 1;
const double Alpha = 1;
const double Beta = 1;

const double ini_time = -44;
const double ini_pos = -11.076;

const double Pea = 5*2.5*10^4/(0.28*4.4*10^3); // V0*Vcyto/(Da*Omega)
const double Pep = 5*2.5*10^4/(0.15*4.4*10^3); // V0*Vcyto/(Dp*Omega)
const double Daoffa = 3.24*10^(-3)*(2.5*10^4)^2/(0.28*(4.4*10^3)^2); // koffa*Vcyto^2/(Da*Omega^2)
const double Daoffp = 7.19*10^(-3)*(2.5*10^4)^2/(0.15*(4.4*10^3)^2); // koffp*Vcyto^2/(Dp*Omega^2)
const double Daona = 6.29*10^(-3)*(2.5*10^4)/(0.28*4.4*10^3); // kona*Vcyto/(Da*Omega)
const double Daonp = 7.682*10^(-2)*(2.5*10^4)/(0.15*4.4*10^3); // konp*Vcyto/(Dp*Omega)
const double Daap = Kap*(2.5*10^4)^2/(0.28*(4.4*10^3)^2); // kap*Vcyto^2/(Da*Omega^2)
const double Dapa = Kpa*(2.5*10^4)^2/(0.15*(4.4*10^3)^2); // kpa*Vcyto^2/(Dp*Omega^2)
const double Lalpha = (9.8*10^4*(4.4*10^3)^2/((2.5*10^4)^2))^Alpha; // (Np*Omega^2/(Vcyto^2))^Alpha
const double Lbeta = (2.4*10^5*(4.4*10^3)^2/((2.5*10^4)^2))^Beta; // (Na*Omega^2/(Vcyto^2))^Beta
const double Ca = 6.29*10^(-3)*(2.5*10^4)^3/(0.28*(4.4*10^3)^4); // kona*Vcyto^3/(Da*Omega^4)
const double Cp = 7.682*10^(-2)*(2.5*10^4)^3/(0.15*(4.4*10^3)^4); // konp*Vcyto^3/(Dp*Omega^4)

void DFsolver(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 2*sin(ini_pos/3.2)*exp(-(((ini_time+step*dt)/10-20)^4)/(10^5));
    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(((ini_time+step*dt)/10-20)^4)/(10^5));
    u_new[0] = ( (1/Pea)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx^2) - (Daoffa-Daona+Daap*Lalpha*(w[0])^Alpha)*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*2*dt + u_old[0];
    w_new[0] = ( (1/Pep)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx^2) - (Daoffp-Daonp+Dapa*Lbeta*(u[0])^Beta)*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*2*dt + w_old[0];
    
    // Second to (Last-1) nodes
    for(i=1;i<Last;i++)
    {
        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(((ini_time+step*dt)/10-20)^4)/(10^5));
        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(((ini_time+step*dt)/10-20)^4)/(10^5));
        u_new[i] = ( (1/Pea)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx^2) - (Daoffa-Daona+Daap*Lalpha*(w[i])^Alpha)*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + u_old[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx^2) - (Daoffp-Daonp+Dapa*Lbeta*(u[i])^Beta)*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + w_old[i]; 
    }
    
    // Last node
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new; 
}

void DFsolver_firststep(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 2*sin(ini_pos/3.2)*exp(-((ini_time/10-20)^4)/(10^5));
    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-((ini_time/10-20)^4)/(10^5));
    u_new[0] = ( (1/Pea)*( (u[1]-2*u[0]+u[Last-1])/(dx^2) - (Daoffa-Daona+Daap*Lalpha*(w[0])^Alpha)*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx )*dt + u[0];
    w_new[0] = ( (1/Pep)*( (w[1]-2*w[0]+w[Last-1])/(dx^2) - (Daoffp-Daonp+Dapa*Lbeta*(u[0])^Beta)*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];
    
    // Second to (Last-1) nodes
    for(i=1;i<Last;i++)
    {
        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-((ini_time/10-20)^4)/(10^5));
        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-((ini_time/10-20)^4)/(10^5));
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*u[i]+u[i-1])/(dx^2) - (Daoffa-Daona+Daap*Lalpha*(w[i])^Alpha)*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*w[i]+w[i-1])/(dx^2) - (Daoffp-Daonp+Dapa*Lbeta*(u[i])^Beta)*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*dt + w[i]; 
    }
    
    // Last node
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void Iterate(int maxsteps, double accuracy, int max_unchanged_steps, int step_interval_writetofile, vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, 
        void (*method)(vector<double>& u_old, vector<double>& w_old,vector<double>& u, vector<double>& w, double dt, double dx, int step)) {
    double elapsed;
    clock_t t0 = clock();
    int step = 1;
    vector<double>& stability_errors(max_unchanged_steps,1);
    double unchanged_error;
    ofstream ofs;
    ofs.open("error_data_temp", ofstream::out);
    
    while (true) 
    {
        
        method(u_old, w_old, u, w, dt, dx, step);
        
        if (step % step_interval_writetofile == 0)
        {
            Writetofile(u,w,step,dt,dx);
            Plot(step, dt, dx);
        }
        
        double error = Relative_error(u_old, w_old, u, w);
        stability_errors.insert(stability_errors.begin(),error);
        stability_errors.pop_back();
        double unchanged_error = std::accumulate(stability_errors.begin(), stability_errors.end(), 0.0);
        unchanged_error /= stability_errors.size();
        if (unchanged_error < accuracy)
        {
            cout<<"Concentrations almost unchanged over the last "<<max_unchanged_steps<<" steps. Stability reached."<<endl;
            break;
        }
        
        step++;
        if (step>=maxsteps)
        {
            cout<<"Max number of steps ("<<maxsteps<<") reached. No stability."<<endl;
            break;
        }

        ofs << step << '\t' << error << '\n';
    }
    ofs.close();
    cout << "Total number of steps = " << step << " steps."<<endl;

    clock_t t1 = clock();
    elapsed = double(t1 - t0) / CLOCKS_PER_SEC;
    cout << "Computational time = " << elapsed << " seconds." << endl;
    
    Writetofile(u,w,step,dt,dx);
    Plot(step, dt, dx);
}

double Relative_error(vector<double> u_old, vector<double> w_old, vector<double> u_new, vector<double> w_new) {

    double error = 0;
    
    vector<double> udiff = u_old;
    std::transform (udiff.begin(), udiff.end(), u_new.begin(), udiff.begin(), std::minus<int>());
    error += sqrt(inner_product(udiff.begin(), udiff.end(), udiff.begin(), 0.0))<<endl;
    
    vector<double> wdiff = w_old;
    std::transform (wdiff.begin(), wdiff.end(), w_new.begin(), wdiff.begin(), std::minus<int>());
    error += sqrt(inner_product(wdiff.begin(), wdiff.end(), wdiff.begin(), 0.0))<<endl;
    
    return error;
}

void Writetofile(vector<double> u, vector<double> w, int step, double dt, double dx) {
    
    stringstream file_namestr;
    file_namestr << "Concentration_A_attime_" << step*dt << "_withdx_" << dx <<".data";
    string file_name = file_namestr.str();

    ofstream ofs;
    ofs.open(file_name.c_str(), ofstream::out);
    for(vector<double>::const_iterator i = u.begin(); i != u.end(); ++i) {ofs << *i << '\n';}
    ofs << '\n';
    ofs.close();
    
    cout<<"Concentrations of A over the 1D membrane have been saved to"<<file_name<<endl;
    
    stringstream file_namestr;
    file_namestr << "Concentration_P_attime_" << step*dt << "_withdx_" << dx <<".data";
    string file_name = file_namestr.str();

    ofstream ofs;
    ofs.open(file_name.c_str(), ofstream::out);
    for(vector<double>::const_iterator i = w.begin(); i != w.end(); ++i) {ofs << *i << '\n';}
    ofs << '\n';
    ofs.close();
    
    cout<<"Concentrations of P over the 1D membrane have been saved to"<<file_name<<endl;
}

void Plot(int step, double dt, double dx){
    stringstream file_namestr;
    file_namestr << "Concentration_A_attime_" << step*dt << "_withdx_" << dx <<".data";
    string file_name = file_namestr.str();
    
    stringstream gnuscript_namestr;
    gnuscript_namestr << "gnuscript_temp";
    string gnuscript_name = gnuscript_namestr.str();
    
    ofstream ofs2;
    ofs2.open(gnuscript_name.c_str(), ofstream::out);
    ofs2 << "set terminal png\n" 
           "set output 'C_A_"<<step*dt<<"_"<<dx<<".png'\n"
           "set title \"Concentration of A along the 1D membrane\"\n"
           "set xlabel \"Position\"\n"
           "set ylabel \"Concentration\"\n"
           "plot '"<<file_name<<"'\n"
           "set terminal x11\n"
           "set output\n"
           "replot\n"
           "pause -1 \"Press to continue\"";
    ofs2.close();
    stringstream terminal_line2;
    terminal_line2 << "gnuplot" << " " << gnuscript_name;
    string command2 = terminal_line2.str();
    system(command2.c_str());
    cout<<"An image of the plot has been saved."<<endl;
    
    stringstream file_namestr;
    file_namestr << "Concentration_P_attime_" << step*dt << "_withdx_" << dx <<".data";
    string file_name = file_namestr.str();
    
    stringstream gnuscript_namestr;
    gnuscript_namestr << "gnuscript_temp";
    string gnuscript_name = gnuscript_namestr.str();
    
    ofstream ofs2;
    ofs2.open(gnuscript_name.c_str(), ofstream::out);
    ofs2 << "set terminal png\n" 
           "set output 'C_P_"<<step*dt<<"_"<<dx<<".png'\n"
           "set title \"Concentration of P along the 1D membrane\"\n"
           "set xlabel \"Position\"\n"
           "set ylabel \"Concentration\"\n"
           "plot '"<<file_name<<"'\n"
           "set terminal x11\n"
           "set output\n"
           "replot\n"
           "pause -1 \"Press to continue\"";
    ofs2.close();
    stringstream terminal_line2;
    terminal_line2 << "gnuplot" << " " << gnuscript_name;
    string command2 = terminal_line2.str();
    system(command2.c_str());
    cout<<"An image of the plot has been saved."<<endl;
}

int main(int argc, char** argv) {
    double dt = 0.01;
    double final_time = 440;
    double maxtime = final_time-ini_time;
    int maxsteps = floor(maxtime/dt);
    
    double dx = 0.006;
    double final_pos = 11.076;
    double maxlength = final_pos-ini_pos;
    int Npoints = floor(maxlength/dx);
    
    double accuracy = 0.001;
    int max_unchanged_steps = 1000;
    int step_interval_writetofile = 100; // every how many step you should write to file
    
    vector<double> u(Npoints,0.5); // concentration of A at each point of the domain, initialized to 0.5
    vector<double> w(Npoints,0.5); // concentration of P at each point of the domain, initialized to 0.5
    vector<double> u_old = u;
    vector<double> w_old = w;
    DFsolver_firststep(u, w, dt, dx);
    Iterate(maxsteps, accuracy, max_unchanged_steps, step_interval_writetofile, u_old, w_old, u, w, dt, dx, DFsolver);
    
    return 0;
}