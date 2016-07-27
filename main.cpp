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
const double Pea = 0.0; // V0*Vcyto/(Da*Omega)
const double Pep = 0.0; // V0*Vcyto/(Dp*Omega)
const double Vx = 0.0; // Vx
const double Daoffa = 0.0; // koffa*Vcyto^2/(Da*Omega^2)
const double Daoffp = 0.0; // koffp*Vcyto^2/(Dp*Omega^2)
const double Daona = 0.0; // kona*Vcyto/(Da*Omega)
const double Daonp = 0.0; // konp*Vcyto/(Da*Omega)
const double Daap = 0.0; // kap*Vcyto^2/(Da*Omega^2)
const double Dapa = 0.0; // kpa*Vcyto^2/(Dp*Omega^2)
const double Alpha = 0.0;
const double Beta = 0.0;
const double Lalpha = 0.0; // (Np*Omega^2/(Vcyto^2))^Alpha
const double Lbeta = 0.0; // (Na*Omega^2/(Vcyto^2))^Beta
const double Ca = 0.0; // kona*Vcyto^3/(Da*Omega^4)
const double Cp = 0.0; // konp*Vcyto^3/(Dp*Omega^4)

void Iterate(int maxsteps, double accuracy, int max_unchanged_steps, int step_interval_writetofile, vector<double>& u, vector<double>& w, double dt, double dx, 
        void (*method)(vector<double>& u, vector<double>& w, double dt, double dx)) {
    double elapsed;
    clock_t t0 = clock();
    int step = 1;
    vector<double>& stability_errors(max_unchanged_steps,1);
    double unchanged_error;
    ofstream ofs;
    ofs.open("error_data_temp", ofstream::out);
    
    while (true) 
    {
        vector<double> u_old = u; 
        vector<double> w_old = w;
        
        method(u,w,dt,dx);
        
        if (step % step_interval_writetofile == 0) Writetofile(u,w,step,dt,dx);
        
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
    ofs.close();
    
    cout<<"Concentrations of A over the 1D membrane have been saved to"<<file_name<<endl;
    
    stringstream file_namestr;
    file_namestr << "Concentration_P_attime_" << step*dt << "_withdx_" << dx <<".data";
    string file_name = file_namestr.str();

    ofstream ofs;
    ofs.open(file_name.c_str(), ofstream::out);
    for(vector<double>::const_iterator i = w.begin(); i != w.end(); ++i) {ofs << *i << '\n';}
    ofs.close();
    
    cout<<"Concentrations of P over the 1D membrane have been saved to"<<file_name<<endl;
}

void Plot(string method){
    stringstream file_namestr;
    file_namestr << "poisson_" << method << ".data";
    string file_name = file_namestr.str();
    
    stringstream gnu_script_namestr;
    gnu_script_namestr << "gnuscript_" << method;
    string gnu_script_name = gnu_script_namestr.str();
    
    ofstream ofs;
    ofs.open(gnu_script_name.c_str(), ofstream::out);
    ofs << "set terminal png\n" 
           "set output 'plot_"<<method<<".png'\n"
           "set title \"Superficie 3D del potenziale V (metodo "<<method<<")\"\n"
           "set xlabel \"x\"\n"
           "set ylabel \"y\"\n"
           "set zlabel \"V(x,y)\"\n"
           "$map2 << EOD\n"; 
    ifstream ifs(file_name.c_str(), ifstream::in);
    string line;
    while (getline(ifs, line)) ofs<<line<<"\n";
    
    ofs << "EOD\n"
           "set ticslevel 0\n"
           "set pm3d at s hidden3d 100\n"
           "set palette rgbformulae 21,22,23\n"  
//           "set logscale cb\n"
//           "set style fill transparent solid 0.9\n"
           "splot '$map2' using 1:2:3 with lines notitle\n"
           "set terminal x11\n"
           "set output\n"
           "replot\n"
           "pause -1 \"Premere INVIO per continuare\"";
    ofs.close();
    
    stringstream terminal_line;
    terminal_line << "gnuplot" << " " << gnu_script_name;
    string command = terminal_line.str();
    system(command.c_str());
    cout<<"Un'immagine del plot è stata salvata nel file plot_"<<method<<".png"<<endl;
    
    if (method=="fft") return;
    
    stringstream gnu_error_namestr;
    gnu_error_namestr << "gnu_error_" << method;
    string gnu_error_name = gnu_error_namestr.str();
    
    ofstream ofs2;
    ofs2.open(gnu_error_name.c_str(), ofstream::out);
    ofs2 << "set terminal png\n" 
           "set output 'error_"<<method<<".png'\n"
           "set title \"Errore vs Iterazioni (metodo "<<method<<")\"\n"
           "set xlabel \"Iterazione\"\n"
           "set ylabel \"Errore (scala logaritmica)\"\n"
           "set logscale y\n"
           "plot 'error_data_temp'\n"
           "set terminal x11\n"
           "set output\n"
           "replot\n"
           "pause -1 \"Premere INVIO per continuare\"";
    ofs2.close();
    stringstream terminal_line2;
    terminal_line2 << "gnuplot" << " " << gnu_error_name;
    string command2 = terminal_line2.str();
    system(command2.c_str());
    cout<<"Un'immagine dell'andamento dell'errore è stata salvata nel file error_"<<method<<".png"<<endl;
}

int main(int argc, char** argv) {
    double dt = 0.01;
    double maxtime = 20;
    int maxsteps = floor(maxtime/dt);
    double dx = 0.01;
    double maxlength = 20;
    int Npoints = floor(maxlength/dx);
    double accuracy = 0.001;
    int max_unchanged_steps = 1000;
    int step_interval_writetofile = 100; // every how many step you should write to file
    vector<double> u(Npoints,0.5); // concentration of A at each point of the domain, initialized to 0.5
    vector<double> w(Npoints,0.5); // concentration of P at each point of the domain, initialized to 0.5
    
    
    return 0;
}