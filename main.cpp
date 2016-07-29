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
//using namespace boost::numeric::ublas;
const double pi = boost::math::constants::pi<double>();

const double Kap = 0.001; // also: (1,2), (0.19,0.3)
const double Kpa = 0.006;
const double Alpha = 1;
const double Beta = 2;

const double ini_time = -44;
const double ini_pos = -16*pi/5; //-11.08; -10.0531;

const double Pea = 5*2.5*10000/(0.28*4.4*1000); // V0*Vcyto/(Da*Omega) 101.461
const double Pep = 5*2.5*10000/(0.15*4.4*1000); // V0*Vcyto/(Dp*Omega) 189.394
const double Daoffa = 3.24*0.001*pow(2.5*10000,2)/(0.28*pow(4.4*1000,2)); // koffa*Vcyto^2/(Da*Omega^2) 0.373561
const double Daoffp = 7.19*0.001*pow(2.5*10000,2)/(0.15*pow(4.4*1000,2)); // koffp*Vcyto^2/(Dp*Omega^2) 1.54743
const double Daona = 6.29*0.001*(2.5*10000)/(0.28*4.4*1000); // kona*Vcyto/(Da*Omega) 0.127638
const double Daonp = 7.682*0.001*(2.5*10000)/(0.15*4.4*1000); // konp*Vcyto/(Dp*Omega) 0.290985
const double Daap = Kap*pow(2.5*10000,2)/(0.28*pow(4.4*1000,2)); // kap*Vcyto^2/(Da*Omega^2) 115.297
const double Dapa = Kpa*pow(2.5*10000,2)/(0.15*pow(4.4*1000,2)); // kpa*Vcyto^2/(Dp*Omega^2) 215.22
const double Lalpha = pow(9.8*10000*pow(4.4*1000,2)/(pow(2.5*10000,2)),Alpha); // (Np*Omega^2/(Vcyto^2))^Alpha 3035.65
const double Lbeta = pow(2.4*100000*pow(4.4*1000,2)/(pow(2.5*10000,2)),Beta); // (Na*Omega^2/(Vcyto^2))^Beta 7434.24
const double Ca = 6.29*0.001*pow(2.5*10000,3)/(0.28*pow(4.4*1000,4)); // kona*Vcyto^3/(Da*Omega^4) 0.000936487
const double Cp = 7.682*0.01*pow(2.5*10000,3)/(0.15*pow(4.4*1000,4)); // konp*Vcyto^3/(Dp*Omega^4) 0.0213497

double Vel(double x){
    double Vel = 1.5*sin(x/3.2);
    return Vel;
}

void InitializeGaussian(vector<double>& u, vector<double>& w, double dx) {
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    
    // gaussian c.i.
    for(int i=0;i<=Last;i++)
    {
        u_new[i] = u[i]*exp(-pow(ini_pos+dx*i,2))/sqrt(2*pi);
    }
    u = u_new;
}

void Writetofile(vector<double>& u, vector<double>& w, int step, double dt, double dx) {
    
    stringstream file_namestrA;
    file_namestrA << "ResultsA/Concentration_A_atstep_" << step << "_withdx_" << dx <<".data";
    string file_nameA = file_namestrA.str();

    ofstream ofsA;
    ofsA.open(file_nameA.c_str(), ofstream::out);
    for(vector<double>::const_iterator i = u.begin(); i != u.end(); ++i) {ofsA << *i << '\n';}
    ofsA << '\n';
    ofsA.close();
    
    cout<<"Concentrations of A over the 1D membrane have been saved to "<<file_nameA<<endl;
    
    stringstream file_namestrP;
    file_namestrP << "ResultsP/Concentration_P_atstep_" << step << "_withdx_" << dx <<".data";
    string file_nameP = file_namestrP.str();

    ofstream ofsP;
    ofsP.open(file_nameP.c_str(), ofstream::out);
    for(vector<double>::const_iterator i = w.begin(); i != w.end(); ++i) {ofsP << *i << '\n';}
    ofsP << '\n';
    ofsP.close();
    
    cout<<"Concentrations of P over the 1D membrane have been saved to "<<file_nameP<<endl;
}

void Plot(int step, double dt, double dx){
    stringstream file_namestrA;
    file_namestrA << "ResultsA/Concentration_A_atstep_" << step << "_withdx_" << dx <<".data";
    string file_nameA = file_namestrA.str();
    
    stringstream gnuscript_namestrA;
    gnuscript_namestrA << "Scripts_and_error/gnuscript_temp";
    string gnuscript_nameA = gnuscript_namestrA.str();
    
    ofstream ofsA;
    ofsA.open(gnuscript_nameA.c_str(), ofstream::out);
    ofsA << "set terminal png\n" 
           "set output 'ImagesA/C_A_"<<step<<"_"<<dx<<".png'\n"
           "set title \"Concentration of A along the 1D membrane\"\n"
           "set xlabel \"Position\"\n"
           "set ylabel \"Concentration\"\n"
//           "set yrange [0:1]\n"
           "plot '"<<file_nameA<<"' pt 7 ps 1\n";
    ofsA.close();
    stringstream terminal_line2A;
    terminal_line2A << "gnuplot" << " " << gnuscript_nameA;
    string command2A = terminal_line2A.str();
    system(command2A.c_str());
    cout<<"An image of the plot has been saved."<<endl;
    
    stringstream file_namestrP;
    file_namestrP << "ResultsP/Concentration_P_atstep_" << step << "_withdx_" << dx <<".data";
    string file_nameP = file_namestrP.str();
    
    stringstream gnuscript_namestrP;
    gnuscript_namestrP << "Scripts_and_error/gnuscript_temp";
    string gnuscript_nameP = gnuscript_namestrP.str();
    
    ofstream ofsP;
    ofsP.open(gnuscript_nameP.c_str(), ofstream::out);
    ofsP << "set terminal png\n" 
           "set output 'ImagesP/C_P_"<<step<<"_"<<dx<<".png'\n"
           "set title \"Concentration of P along the 1D membrane\"\n"
           "set xlabel \"Position\"\n"
           "set ylabel \"Concentration\"\n"
           "plot '"<<file_nameP<<"' pt 7 ps 1\n";
    ofsP.close();
    stringstream terminal_line2P;
    terminal_line2P << "gnuplot" << " " << gnuscript_nameP;
    string command2P = terminal_line2P.str();
    system(command2P.c_str());
    cout<<"An image of the plot has been saved."<<endl;
}

void DFsolver(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
//    Vx = 2*sin(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
//    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(1000000));
    Vx = 1.5*sin(ini_pos/3.2);
    dV_dx = (1.5/3.2)*cos(ini_pos/3.2);
    
    u_new[0] = ( (1/Pea)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*2*dt + u_old[0];
    w_new[0] = ( (1/Pep)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*2*dt + w_old[0];

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
//        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
//        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        Vx = 1.5*sin((ini_pos+i*dx)/3.2);
        dV_dx = (1.5/3.2)*cos((ini_pos+i*dx)/3.2);
        u_new[i] = ( (1/Pea)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + u_old[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*2*dt + w_old[i];
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
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
//    Vx = 2*sin(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4)/(100000)));
    Vx = 1.5*sin(ini_pos/3.2);
    dV_dx = (1.5/3.2)*cos(ini_pos/3.2);

    cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos<<endl;
    
    u_new[0] = ( (1/Pea)*( (u[1]-2*u[0]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx )*dt + u[0];
    w_new[0] = ( (1/Pep)*( (w[1]-2*w[0]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];
    

    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 1.5*sin((ini_pos+i*dx)/3.2);
        dV_dx = (1.5/3.2)*cos((ini_pos+i*dx)/3.2);
//        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
        cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos+i*dx<<endl;
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*u[i]+u[i-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*w[i]+w[i-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*dt + w[i]; 

    }
    
    // Last node 
//    Vx = 2*sin((ini_pos+Last*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    dV_dx = 0.625*cos((ini_pos+Last*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos+Last*dx<<endl;
        
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[Last]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[Last]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void DFsolver_AvgCon(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    double u_avg = std::accumulate(u.begin(), u.end(), 0.0);
    u_avg /= u.size();
    double w_avg = std::accumulate(w.begin(), w.end(), 0.0);
    w_avg /= w.size();
    
    // First node
//    Vx = 2*sin(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
//    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(1000000));
    Vx = 1.5*sin(ini_pos/3.2);
    dV_dx = (1.5/3.2)*cos(ini_pos/3.2);
    
    u_new[0] = ( (1/Pea)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca-Daona*u_avg) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*2*dt + u_old[0];
    w_new[0] = ( (1/Pep)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp-Daonp*w_avg) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*2*dt + w_old[0];

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
//        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
//        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        Vx = 1.5*sin((ini_pos+i*dx)/3.2);
        dV_dx = (1.5/3.2)*cos((ini_pos+i*dx)/3.2);
        u_new[i] = ( (1/Pea)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca-Daona*u_avg) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + u_old[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp-Daonp*w_avg) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*2*dt + w_old[i];
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

void DFsolver_firststep_AvgCon(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    double u_avg = std::accumulate(u.begin(), u.end(), 0.0);
    u_avg /= u.size();
    double w_avg = std::accumulate(w.begin(), w.end(), 0.0);
    w_avg /= w.size();
    
    // First node
//    Vx = 2*sin(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4)/(100000)));
    Vx = 1.5*sin(ini_pos/3.2);
    dV_dx = (1.5/3.2)*cos(ini_pos/3.2);

    cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos<<endl;
    
    u_new[0] = ( (1/Pea)*( (u[1]-2*u[0]+u[Last-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca-Daona*u_avg) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx )*dt + u[0];
    w_new[0] = ( (1/Pep)*( (w[1]-2*w[0]+w[Last-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp-Daonp*w_avg) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];
    

    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 1.5*sin((ini_pos+i*dx)/3.2);
        dV_dx = (1.5/3.2)*cos((ini_pos+i*dx)/3.2);
//        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
        cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos+i*dx<<endl;
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*u[i]+u[i-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca-Daona*u_avg) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*w[i]+w[i-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp-Daonp*w_avg) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*dt + w[i]; 

    }
    
    // Last node 
//    Vx = 2*sin((ini_pos+Last*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    dV_dx = 0.625*cos((ini_pos+Last*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos+Last*dx<<endl;
        
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[Last]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[Last]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void DFsolver_AvgConUpw(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    double u_avg = std::accumulate(u.begin(), u.end(), 0.0);
    u_avg /= u.size();
    double w_avg = std::accumulate(w.begin(), w.end(), 0.0);
    w_avg /= w.size();
    
    // First node
//    Vx = 2*sin(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
//    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(1000000));
    Vx = 1.5*sin(ini_pos/3.2);
    dV_dx = (1.5/3.2)*cos(ini_pos/3.2);
    
    u_new[0] = ( (1/Pea)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca-Daona*u_avg) - ((Vel(ini_pos)>=0.0)*(u[0]*Vel(ini_pos)-u[Last-1]*Vel(ini_pos+(Last-1)*dx))/dx+(Vel(ini_pos)<0.0)*(u[1]*Vel(ini_pos+dx)-u[0]*Vel(ini_pos))/dx))*2*dt + u_old[0];
    w_new[0] = ( (1/Pep)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp-Daonp*w_avg) - ((Vel(ini_pos)>=0.0)*(w[0]*Vel(ini_pos)-w[Last-1]*Vel(ini_pos+(Last-1)*dx))/dx+(Vel(ini_pos)<0.0)*(w[1]*Vel(ini_pos+dx)-w[0]*Vel(ini_pos))/dx))*2*dt + w_old[0];

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
//        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
//        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        Vx = 1.5*sin((ini_pos+i*dx)/3.2);
        dV_dx = (1.5/3.2)*cos((ini_pos+i*dx)/3.2);
        u_new[i] = ( (1/Pea)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca-Daona*u_avg) - ((Vel(ini_pos+i*dx)>=0.0)*(u[i]*Vel(ini_pos+i*dx)-u[i-1]*Vel(ini_pos+(i-1)*dx))/dx+(Vel(ini_pos+i*dx)<0.0)*(u[i+1]*Vel(ini_pos+(i+1)*dx)-u[i]*Vel(ini_pos+i*dx))/dx))*2*dt + u_old[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp-Daonp*w_avg) - ((Vel(ini_pos+i*dx)>=0.0)*(w[i]*Vel(ini_pos+i*dx)-w[i-1]*Vel(ini_pos+(i-1)*dx))/dx+(Vel(ini_pos+i*dx)<0.0)*(w[i+1]*Vel(ini_pos+(i+1)*dx)-w[i]*Vel(ini_pos+i*dx))/dx))*2*dt + w_old[i];
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

void DFsolver_firststep_AvgConUpw(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    double u_avg = std::accumulate(u.begin(), u.end(), 0.0);
    u_avg /= u.size();
    double w_avg = std::accumulate(w.begin(), w.end(), 0.0);
    w_avg /= w.size();
    
    // First node
//    Vx = 2*sin(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4)/(100000)));
    Vx = 1.5*sin(ini_pos/3.2);
    dV_dx = (1.5/3.2)*cos(ini_pos/3.2);
    
    Vx = 1.5*sin(ini_pos/3.2);

    cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos<<endl;
    
    u_new[0] = ( (1/Pea)*( (u[1]-2*u[0]+u[Last-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca-Daona*u_avg) - ((Vel(ini_pos)>=0.0)*(u[0]*Vel(ini_pos)-u[Last-1]*Vel(ini_pos+(Last-1)*dx))/dx+(Vel(ini_pos)<0.0)*(u[1]*Vel(ini_pos+dx)-u[0]*Vel(ini_pos))/dx) )*dt + u[0];
    w_new[0] = ( (1/Pep)*( (w[1]-2*w[0]+w[Last-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp-Daonp*w_avg) - ((Vel(ini_pos)>=0.0)*(w[0]*Vel(ini_pos)-w[Last-1]*Vel(ini_pos+(Last-1)*dx))/dx+(Vel(ini_pos)<0.0)*(w[1]*Vel(ini_pos+dx)-w[0]*Vel(ini_pos))/dx) )*dt + w[0];
    

    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 1.5*sin((ini_pos+i*dx)/3.2);
        dV_dx = (1.5/3.2)*cos((ini_pos+i*dx)/3.2);
//        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
        cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos+i*dx<<endl;
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*u[i]+u[i-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca-Daona*u_avg) - ((Vel(ini_pos+i*dx)>=0.0)*(u[i]*Vel(ini_pos+i*dx)-u[i-1]*Vel(ini_pos+(i-1)*dx))/dx+(Vel(ini_pos+i*dx)<0.0)*(u[i+1]*Vel(ini_pos+(i+1)*dx)-u[i]*Vel(ini_pos+i*dx))/dx))*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*w[i]+w[i-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp-Daonp*w_avg) - ((Vel(ini_pos+i*dx)>=0.0)*(w[i]*Vel(ini_pos+i*dx)-w[i-1]*Vel(ini_pos+(i-1)*dx))/dx+(Vel(ini_pos+i*dx)<0.0)*(w[i+1]*Vel(ini_pos+(i+1)*dx)-w[i]*Vel(ini_pos+i*dx))/dx))*dt + w[i]; 

    }
    
    // Last node 
//    Vx = 2*sin((ini_pos+Last*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    dV_dx = 0.625*cos((ini_pos+Last*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos+Last*dx<<endl;
        
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[Last]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[Last]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void DFsolver_AvgConLax(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    double u_avg = std::accumulate(u.begin(), u.end(), 0.0);
    u_avg /= u.size();
    double w_avg = std::accumulate(w.begin(), w.end(), 0.0);
    w_avg /= w.size();
    
    // First node
//    Vx = 2*sin(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
//    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(1000000));
    Vx = 1.5*sin(ini_pos/3.2);
    dV_dx = (1.5/3.2)*cos(ini_pos/3.2);
    
    u_new[0] = ( (1/Pea)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca-Daona*u_avg) - Vx*((u[1]-u[Last-1])/2/dx +(u[1]-2*u[0]+u[Last-1])/2/dx/dx*dt) - u[0]*dV_dx )*dt + (u_old[1]+u_old[Last-1])/2;
    w_new[0] = ( (1/Pep)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp-Daonp*w_avg) - Vx*((w[1]-w[Last-1])/2/dx +(w[1]-2*w[0]+w[Last-1])/2/dx/dx*dt) - w[0]*dV_dx )*dt + (w_old[1]+w_old[Last-1])/2;

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
//        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
//        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        Vx = 1.5*sin((ini_pos+i*dx)/3.2);
        dV_dx = (1.5/3.2)*cos((ini_pos+i*dx)/3.2);
        u_new[i] = ( (1/Pea)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca-Daona*u_avg) - Vx*((u[i+1]-u[i-1])/2/dx +(u[i+1]-2*u[i]+u[i-1])/2/dx/dx*dt)  - u[i]*dV_dx )*dt + (u_old[i+1]+u_old[i-1])/2;
        w_new[i] = ( (1/Pep)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp-Daonp*w_avg) - Vx*((w[i+1]-w[i-1])/2/dx +(w[i+1]-2*w[i]+w[i-1])/2/dx/dx*dt)  - w[i]*dV_dx )*dt + (w_old[i+1]+w_old[i-1])/2;
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

void DFsolver_firststep_AvgConLax(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    double u_avg = std::accumulate(u.begin(), u.end(), 0.0);
    u_avg /= u.size();
    double w_avg = std::accumulate(w.begin(), w.end(), 0.0);
    w_avg /= w.size();
    
    // First node
//    Vx = 2*sin(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4)/(100000)));
    Vx = 1.5*sin(ini_pos/3.2);
    dV_dx = (1.5/3.2)*cos(ini_pos/3.2);
    
    Vx = 1.5*sin(ini_pos/3.2);

    cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos<<endl;
    
    u_new[0] = ( (1/Pea)*( (u[1]-2*u[0]+u[Last-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca-Daona*u_avg) - ((Vel(ini_pos)>=0.0)*(u[0]*Vel(ini_pos)-u[Last-1]*Vel(ini_pos+(Last-1)*dx))/dx+(Vel(ini_pos)<0.0)*(u[1]*Vel(ini_pos+dx)-u[0]*Vel(ini_pos))/dx) )*dt + (u[1]+u[Last-1])/2;
    w_new[0] = ( (1/Pep)*( (w[1]-2*w[0]+w[Last-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp-Daonp*w_avg) - ((Vel(ini_pos)>=0.0)*(w[0]*Vel(ini_pos)-w[Last-1]*Vel(ini_pos+(Last-1)*dx))/dx+(Vel(ini_pos)<0.0)*(w[1]*Vel(ini_pos+dx)-w[0]*Vel(ini_pos))/dx) )*dt + w[0];
    

    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 1.5*sin((ini_pos+i*dx)/3.2);
        dV_dx = (1.5/3.2)*cos((ini_pos+i*dx)/3.2);
//        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
        cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos+i*dx<<endl;
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*u[i]+u[i-1])/(dx*dx) - (Daoffa+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca-Daona*u_avg) - ((Vel(ini_pos+i*dx)>=0.0)*(u[i]*Vel(ini_pos+i*dx)-u[i-1]*Vel(ini_pos+(i-1)*dx))/dx+(Vel(ini_pos+i*dx)<0.0)*(u[i+1]*Vel(ini_pos+(i+1)*dx)-u[i]*Vel(ini_pos+i*dx))/dx))*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*w[i]+w[i-1])/(dx*dx) - (Daoffp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp-Daonp*w_avg) - ((Vel(ini_pos+i*dx)>=0.0)*(w[i]*Vel(ini_pos+i*dx)-w[i-1]*Vel(ini_pos+(i-1)*dx))/dx+(Vel(ini_pos+i*dx)<0.0)*(w[i+1]*Vel(ini_pos+(i+1)*dx)-w[i]*Vel(ini_pos+i*dx))/dx))*dt + w[i]; 

    }
    
    // Last node 
//    Vx = 2*sin((ini_pos+Last*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    dV_dx = 0.625*cos((ini_pos+Last*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
//    cout<<"Vx = "<<Vx<<" dV = "<<dV_dx<<" and "<<ini_pos+Last*dx<<endl;
        
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[Last]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[Last]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void DFsolver_norea(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 2*sin(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(1000000));
    
    u_new[0] = ( (1/Pea)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*2*dt + u_old[0];
    w_new[0] = ( (1/Pep)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*2*dt + w_old[0];

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        u_new[i] = ( (1/Pea)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx*dx) ) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + u_old[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx*dx) ) - Vx*(w[i+1]-w[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + w_old[i];
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

void DFsolver_firststep_norea(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 2*sin(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4)/(100000)));    
    
    u_new[0] = ( (1/Pea)*( (u[1]-2*u[0]+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx )*dt + u[0];
    w_new[0] = ( (1/Pep)*( (w[1]-2*w[0]+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];
    
    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*u[i]+u[i-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*w[i]+w[i-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*dt + w[i]; 
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void DFsolver_noadv(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 0;
    dV_dx = 0;
    
    u_new[0] = ( (1/Pea)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*2*dt + u_old[0];
    w_new[0] = ( (1/Pep)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*2*dt + w_old[0];

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 0;
        dV_dx = 0;
        u_new[i] = ( (1/Pea)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + u_old[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + w_old[i];
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
//    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

void DFsolver_firststep_noadv(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 0;
    dV_dx = 0;    
    
    u_new[0] = ( (1/Pea)*( (u[1]-2*u[0]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx )*dt + u[0];
    w_new[0] = ( (1/Pep)*( (w[1]-2*w[0]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];

    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 0;
        dV_dx = 0;
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*u[i]+u[i-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*w[i]+w[i-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*dt + w[i]; 

    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[Last]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[Last]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void DFsolver_onlydiff(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 0;
    dV_dx = 0;
    
    u_new[0] = ( (1/Pea)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*2*dt + u_old[0];
    w_new[0] = ( (1/Pep)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*2*dt + w_old[0];

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 0;
        dV_dx = 0;
        u_new[i] = ( (1/Pea)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx*dx) ) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + u_old[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx*dx) ) - Vx*(w[i+1]-w[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + w_old[i];
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
//    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

void DFsolver_firststep_onlydiff(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 0;
    dV_dx = 0;    
    
    u_new[0] = ( (1/Pea)*( (u[1]-2*u[0]+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx )*dt + u[0];
    w_new[0] = ( (1/Pep)*( (w[1]-2*w[0]+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];

    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 0;
        dV_dx = 0;
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*u[i]+u[i-1])/(dx*dx) ) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*w[i]+w[i-1])/(dx*dx) ) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*dt + w[i]; 
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[Last]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[Last]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void DFsolver_onlyrea(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 0;
    dV_dx = 0;
    
    u_new[0] = ( (1/Pea)*( - (Daoffa+Daona+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*2*dt + u_old[0];
    w_new[0] = ( (1/Pep)*( - (Daoffp+Daonp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*2*dt + w_old[0];

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 0;
        dV_dx = 0;
        u_new[i] = ( (1/Pea)*( - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + u_old[i];
        w_new[i] = ( (1/Pep)*( - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + w_old[i];
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
//    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

void DFsolver_firststep_onlyrea(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 0;
    dV_dx = 0;    
    
    u_new[0] = ( (1/Pea)*( - (Daoffa+Daona+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx )*dt + u[0];
    w_new[0] = ( (1/Pep)*( - (Daoffp+Daonp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];
    
    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 0;
        dV_dx = 0;
        u_new[i] = ( (1/Pea)*( - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (1/Pep)*( - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*dt + w[i]; 
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[Last]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[Last]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void DFsolver_onlyadv(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 2*sin(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(1000000));
    
    u_new[0] = ( (0)*( (u[1]-(u[0]+u_old[0])+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*2*dt + u_old[0];
    w_new[0] = ( (0)*( (w[1]-(w[0]+w_old[0])+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*2*dt + w_old[0];
    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        u_new[i] = ( (0)*( (u[i+1]-(u[i]+u_old[i])+u[i-1])/(dx*dx) ) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + u_old[i];
        w_new[i] = ( (0)*( (w[i+1]-(w[i]+w_old[i])+w[i-1])/(dx*dx) ) - Vx*(w[i+1]-w[i-1])/(2*dx) -u[i]*dV_dx)*2*dt + w_old[i];
        
    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-(u[Last]+u_old[Last])+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-(w[Last]+w_old[Last])+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

void DFsolver_firststep_onlyadv(vector<double>& u, vector<double>& w, double dt, double dx){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 2*sin(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow(ini_time/10-20,4)/(100000)));    
    
    u_new[0] = ( (0)*( (u[1]-2*u[0]+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx )*dt + u[0];
    w_new[0] = ( (0)*( (w[1]-2*w[0]+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];
    
    
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow(ini_time/10-20,4))/(100000));
        u_new[i] = ( (0)*( (u[i+1]-2*u[i]+u[i-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (0)*( (w[i+1]-2*w[i]+w[i-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -w[i]*dV_dx)*dt + w[i]; 

    }
    
    // Last node 
    u_new[Last] = u_new[0];
    w_new[Last] = w_new[0];
    
//    // Last node 
//    u_new[Last] = ( (1/Pea)*( (u[1]-2*u[Last]+u[Last-1])/(dx*dx) ) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx )*dt + u[Last];
//    w_new[Last] = ( (1/Pep)*( (w[1]-2*w[Last]+w[Last-1])/(dx*dx) ) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*dt + w[Last];
    
    // Updating u and w
    u = u_new;
    w = w_new;
}

void Simplersolver(vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, int step){
    
    vector<double> u_new(u.size());
    vector<double> w_new(w.size());
    int Last = u.size()-1;
    double Vx;
    double dV_dx;
    
    // First node
    Vx = 2*sin(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
    dV_dx = 0.625*cos(ini_pos/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(1000000));
    
    u_new[0] = ( (1/Pea)*( (u[1]-2*(u[0])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[0],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[0]*dV_dx)*dt + u[0];
    w_new[0] = ( (1/Pep)*( (w[1]-2*(w[0])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[0],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[0]*dV_dx)*dt + w[0];

    
//    // PROBE! ANYTHING WRONG HERE?
//    cout<<"HERE? "<<u_new[0]<<endl;
//    cout<<"HERE? "<<w_new[0]<<endl;
//    exit(0);
//    // END OF PROBE
    // Second to (Last-1) nodes
    for(int i=1;i<Last;i++)
    {
        Vx = 2*sin((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        dV_dx = 0.625*cos((ini_pos+i*dx)/3.2)*exp(-(pow((ini_time+step*dt)/10-20,4))/(100000));
        u_new[i] = ( (1/Pea)*( (u[i+1]-2*(u[i])+u[i-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[i],Alpha))*u[i]+Ca) - Vx*(u[i+1]-u[i-1])/(2*dx) -u[i]*dV_dx)*dt + u[i];
        w_new[i] = ( (1/Pep)*( (w[i+1]-2*(w[i])+w[i-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[i],Beta))*w[i]+Cp) - Vx*(w[i+1]-w[i-1])/(2*dx) -u[i]*dV_dx)*dt + w[i];
        
    }
    
//    // Last node 
//    u_new[Last] = u_new[0];
//    w_new[Last] = w_new[0];
    
    // Last node 
    u_new[Last] = ( (1/Pea)*( (u[1]-2*(u[Last])+u[Last-1])/(dx*dx) - (Daoffa+Daona+Daap*Lalpha*pow(w[Last],Alpha))*u[0]+Ca) - Vx*(u[1]-u[Last-1])/(2*dx) -u[Last]*dV_dx)*2*dt + u_old[Last];
    w_new[Last] = ( (1/Pep)*( (w[1]-2*(w[Last])+w[Last-1])/(dx*dx) - (Daoffp+Daonp+Dapa*Lbeta*pow(u[Last],Beta))*w[0]+Cp) - Vx*(w[1]-w[Last-1])/(2*dx) -w[Last]*dV_dx)*2*dt + w_old[Last];
    
    // Updating u_old, w_old, u, w
    u_old = u;
    w_old = w;
    u = u_new;
    w = w_new;
}

double Relative_error(vector<double> u_old, vector<double> w_old, vector<double> u_new, vector<double> w_new) {

    double error = 0;
    
    vector<double> udiff = u_old;
    std::transform (udiff.begin(), udiff.end(), u_new.begin(), udiff.begin(), std::minus<int>());
    error += sqrt(inner_product(udiff.begin(), udiff.end(), udiff.begin(), 0.0));
    
    vector<double> wdiff = w_old;
    std::transform (wdiff.begin(), wdiff.end(), w_new.begin(), wdiff.begin(), std::minus<int>());
    error += sqrt(inner_product(wdiff.begin(), wdiff.end(), wdiff.begin(), 0.0));
    
    return error;
}

void Iterate(int maxsteps, int step_interval_writetofile, vector<double>& u_old, vector<double>& w_old, vector<double>& u, vector<double>& w, double dt, double dx, 
        void (*method)(vector<double>& u_old, vector<double>& w_old,vector<double>& u, vector<double>& w, double dt, double dx, int step)) {
    double elapsed;
    clock_t t0 = clock();
    int step = 1;
//    vector<double> stability_errors(max_unchanged_steps,1.0);
//    double unchanged_error;
//    ofstream ofs;
//    ofs.open("Scripts_and_error/error_data_temp", ofstream::out);
    Writetofile(u,w,step,dt,dx);
    Plot(step, dt, dx);
    while (true) 
    {
        
        method(u_old, w_old, u, w, dt, dx, step);
        
        if (step % step_interval_writetofile == 0)
        {
            Writetofile(u,w,step,dt,dx);
            Plot(step, dt, dx);
        }
        
//        double error = Relative_error(u_old, w_old, u, w);
//        stability_errors.insert(stability_errors.begin(),error);
//        stability_errors.pop_back();
//        double unchanged_error = std::accumulate(stability_errors.begin(), stability_errors.end(), 0.0);
//        unchanged_error /= stability_errors.size();
//        if (unchanged_error < accuracy)
//        {
//            cout<<"Concentrations almost unchanged over the last "<<max_unchanged_steps<<" steps. Stability reached."<<endl;
//            break;
//        }
        
        step++;
        if (step>=maxsteps)
        {
            cout<<"Max number of steps ("<<maxsteps<<") reached. No stability."<<endl;
            break;
        }

//        ofs << step << '\t' << error << '\n';
    }
//    ofs.close();
    cout << "Total number of steps = " << step << " steps."<<endl;

    clock_t t1 = clock();
    elapsed = double(t1 - t0) / CLOCKS_PER_SEC;
    cout << "Computational time = " << elapsed << " seconds." << endl;
    
    Writetofile(u,w,step,dt,dx);
    Plot(step, dt, dx);
}

int main(int argc, char** argv) {
    double dt = 0.000001;
    double final_time = 440;
    double maxtime = final_time-ini_time;
    int maxsteps = floor(maxtime/dt); // 100000;//
    
    double final_pos = 16*pi/5;
    double maxlength = final_pos-ini_pos;
    int Npoints = 200;
    double dx = maxlength/Npoints;
    
//    double accuracy = 0.00006;
//    int max_unchanged_steps = 100000;
    int step_interval_writetofile = 10000; // every how many step you should write to file
    
    vector<double> u(Npoints+1,1.045*0.001); // concentration of A at each point of the domain, initialized to 1.045*0.001
    vector<double> w(Npoints+1,0.235*0.001); // concentration of P at each point of the domain, initialized to 0.235*0.001
    //InitializeGaussian(u, w, dx);
//    InitializeGaussian(u, w, dx);
    vector<double> u_old = u;
    vector<double> w_old = w;
    //DFsolver_firststep_AvgConUpw(u, w, dt, dx);
    Iterate(maxsteps, step_interval_writetofile, u_old, w_old, u, w, dt, dx, DFsolver_AvgConLax);
    
    return 0;
}