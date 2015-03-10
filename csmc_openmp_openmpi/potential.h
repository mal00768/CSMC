/** \file potential.h
 \brief potential class
 \ingroup gp_nr
 */

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include <cmath>
#include <iostream>
#include <pot/discrete.cxx>
using pot::DiscreteSolutionEigenSystem;

//this is a struct that creates a square well potential of finite depth
struct FiniteSquareWellPotential {

    double twicemass; //2*m/\hbar^2
    double R; //radius of the potential in fm.
    double V0; //depth of the potential in MeV;
   
    double operator()(double r) {
        if(r < R) { return V0; }
        else { return 0; }
    }
    
};

//this is a struct that creates a wood saxon potential of finite depth
struct WoodsSaxonPotential {

    double twicemass; //2*m/\hbar^2
    double R; //radius of the potential in fm.
    double a;//diffusiveness
    double V0; //depth of the potential in MeV;
   
    double operator()(double r) {
        return V0/(1 + exp((r - R)/a));
    }

};

double factorial(int n){
    
    //Returns the value n! as a floating-point number.
    static std::vector<double> a(171, 0);
    static bool init = true;
    if(init){
        init = false;
        a[0] = 1.0;
        for(int i = 1; i < 171; ++i) { a[i] = i*a[i - 1]; }
    }
    if (n < 0 || n > 170) { throw("factrl out of range"); }
    return a[n];
    
}

double plegendre(const int l, const int m, const double x) {
    static const double PI = 3.141592653589793;
    int i, ll;
    double fact, oldfact, pll, pmm, pmmp1, omx2;
    if (m < 0 || m > l || abs(x) > 1.0)
        throw("Bad arguments in routine plgndr");
    pmm = 1.0;
    if (m > 0) {
        omx2 = (1.0 - x)*(1.0 + x);
        fact = 1.0;
        for(i = 1;i <= m; ++i){
            pmm *= omx2*fact/(fact + 1.0);
            fact += 2.0; }
    }
    pmm = sqrt((2*m + 1)*pmm/(4.0*PI));
    
    if (m & 1)
        pmm = -pmm;
    if (l == m)
        return pmm;
    else {
        pmmp1 = x*sqrt(2.0*m + 3.0)*pmm;
        if (l == (m + 1))
            return pmmp1;
        else {
            oldfact = sqrt(2.0*m + 3.0);
            for(ll = m + 2; ll <= l; ++ll) {
                fact = sqrt((4.0*ll*ll - 1.0)/(ll*ll - m*m));
                pll = (x*pmmp1 - pmm/oldfact)*fact;
                oldfact = fact;
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}

double plgndr(const int l, const int m, const double x) {
    const double PI = 3.141592653589793238;
    if (m < 0 || m > l || abs(x) > 1.0)
        throw("Bad arguments in routine plgndr");
    double prod = 1.0;
    for (int j = l - m + 1; j <= l + m; ++j)
        prod *= j;
    return sqrt(4.0*PI*prod/(2*l + 1))*plegendre(l, m, x);
}

struct Quadrature{
    
    int n; //Current level of refinement.
    virtual double next() = 0;
    
};

template<class T>
struct Trapzd : Quadrature {
    
    double a,b,s;
    T& func;
    Trapzd() {};
    Trapzd(T &funcc, const double aa, const double bb) : func(funcc), a(aa), b(bb) { n = 0; }
    
    double next() {
        double x, tnm, sum, del;
        int it, j;
        ++n;
        if (n == 1) {
            return (s = 0.5*(b - a)*(func(a) + func(b)));
        } else {
            for (it = 1, j = 1; j < n - 1; ++j) it <<= 1;
            tnm = it;
            del = (b - a)/tnm;
            x = a + 0.5*del;
            for(sum = 0.0, j = 0; j < it; ++j, x += del) sum += func(x);
            s = 0.5*(s + (b - a)*sum/tnm);
            return s;
        }
    }
    
};

template<class T>
double qsimp(T& func, const double a, const double b, const double eps = 1.0e-10) {
    const int JMAX = 20;
    double s, st, ost = 0.0, os = 0.0;
    Trapzd<T> t(func, a, b);
    for(int j = 0; j < JMAX; ++j) {
        st = t.next();
        s = (4.0*st - ost)/3.0;
        if(j > 5)
            if (abs(s-os) < eps*abs(os) || (s == 0.0 && os == 0.0)) return s;
        os = s;
        ost = st;
    }
    throw("Too many steps in routine qsimp");
    
}

struct AngularPart {
    
    AngularPart(double ll1, double mm1, double ll2, double mm2) : l1(ll1), m1(mm1), l2(ll2), m2(mm2) {}
    
    double prefactor(double l, double m) {
        const double PI = 3.141592653589793238;
        double a = (2*l + 1)/(4*PI);
        a *= factorial(l - m)/factorial(l + m);
        return sqrt(a);
    }
    
    double operator()(double theta) {
        double p = prefactor(l1, m1)*prefactor(l1, m1)*prefactor(l2, m2)*prefactor(l2, m2);
        
        if(m1 >= 0 && m2 >= 0){
            return p*plgndr(l1, m1, cos(theta))*plgndr(l1, m1, cos(theta))*plgndr(l2, m2, cos(theta))*plgndr(l2, m2, cos(theta))*sin(theta);
        }
        else if(m1 < 0 && m2 >= 0){
            m1 = abs(m1);
            return p*plgndr(l1, m1, cos(theta))*(factorial(l1 - m1)/factorial(l1 + m1))*plgndr(l1, m1, cos(theta))*plgndr(l2, m2, cos(theta))*plgndr(l2, m2, cos(theta))*(factorial(l1 - m2)/factorial(l1 + m2))*sin(theta);
        }
        else if(m1 >= 0 && m2 < 0){
            m2 = abs(m2);
            return p*plgndr(l1, m1, cos(theta))*(factorial(l1 - m1)/factorial(l1 + m1))*plgndr(l1, m1, cos(theta))*plgndr(l2, m2, cos(theta))*plgndr(l2, m2, cos(theta))*(factorial(l1 - m2)/factorial(l1 + m2))*sin(theta);
        }
        else if(m1 < 0 && m2 < 0){
            m1 = abs(m1);
            m2 = abs(m2);
            return p*plgndr(l1, m1, cos(theta))*(factorial(l1 - m1)/factorial(l1 + m1))*plgndr(l1, m1, cos(theta))*plgndr(l2, m2, cos(theta))*plgndr(l2, m2, cos(theta))*(factorial(l1 - m2)/factorial(l1 + m2))*sin(theta);
        }
        else{
            throw("m1 or m2 is not a number in AngularPart");
        }
        
    }
    
private:
    double l1, l2;
    double m1, m2;
    
};


class Potential {
    
public:
    
    Potential(int, double, double, double, double, double, double, double);
    Potential(int, double, double, double, double, double, double, double, double);
    ~Potential();
    void SolveSquareWell();
    void SolveWoodSaxon();
    void PairingMatrixElements();
    
    std::vector<double> Energies; //vector to store energies for states
    std::vector< std::vector<double> > Gij; //store all possible values
    
private:
    
    int nx;
    int numlevels;
    double xmin, xmax;
    double m;
    double twicemass;
    double R;
    double a;
    double V0;
    double hc;
    double v;
    
    double** wf;
    
    double f(int);
    
};

Potential::Potential(int nlevels, double xminn, double xmaxx, double v0, double v1, double r, double mm, double n)
: numlevels(nlevels), xmin(xminn), xmax(xmaxx), V0(v0), v(v1), R(r), m(mm), nx(n) {
   
    //create a double array (matrix of size nx times nx+5)
    wf = new double* [nx];
    for(int i = 0; i < nx; ++i) { wf[i] = new double [nx + 5]; }
    
    Energies.resize(numlevels, 0);
    std::vector<double> tmpG(numlevels, 0);
    for(int i = 0; i < numlevels; ++i){
        Gij.push_back(tmpG);
    }
    
    hc = 1.0;
    
}

Potential::Potential(int nlevels, double xminn, double xmaxx, double v0, double v1, double r, double aa, double mm, double n)
: numlevels(nlevels), xmin(xminn), xmax(xmaxx), V0(v0), v(v1), R(r), a(aa), m(mm), nx(n) {
    
    //create a double array (matrix of size nx times nx+5)
    wf = new double* [nx];
    for(int i = 0; i < nx; ++i) { wf[i] = new double [nx + 5]; }
    
    Energies.resize(numlevels, 0);
    std::vector<double> tmpG(numlevels, 0);
    for(int i = 0; i < numlevels; ++i){
        Gij.push_back(tmpG);
    }
    
    hc = 197.327;
    
}

Potential::~Potential() {
    
    //reclaim the memory
    for (int i = 0; i < nx; ++i) { delete [] wf[i]; }
    delete [] wf;
    
    Energies.clear();
    Gij.clear();
    
}

void Potential::SolveSquareWell() {
    
    FiniteSquareWellPotential V;
    
    V.R = R;
    V.V0 = V0;
    V.twicemass = 2.0*m/(hc*hc); //define mass
    
    //solve the problem
    DiscreteSolutionEigenSystem(nx, wf, xmin, xmax, V, false, false, 1.0/V.twicemass);
  
    for(int i = 0; i < numlevels; ++i) {
        Energies[i] = wf[i][0];
    }
    
    PairingMatrixElements();
    
}

void Potential::SolveWoodSaxon() {
    
    WoodsSaxonPotential V; //create variable of potential type
    
    V.R = R;
    V.a = a;
    V.V0 = V0;
    V.twicemass = 2.0*m/(hc*hc); //define mass
    
    //solve the problem
    DiscreteSolutionEigenSystem(nx, wf, xmin, xmax, V, false, false, 1.0/V.twicemass);
    
    for(int i = 0; i < numlevels; ++i) {
        Energies[i] = wf[i][0];
    }
    
    PairingMatrixElements();
    
}

void Potential::PairingMatrixElements() {
    
    std::vector<double> temp(nx, 0); //save each vector seperately and push back
    //std::vector<double> rho; //in case ρ(r) is known
    std::vector< std::vector<double> > WaveFunctions; //store u(r) for each state
    
    for(int i = 0; i < numlevels; ++i){
        for(int j = 0; j < nx; ++j){
            temp[j] = wf[i][j + 1];
        }
        WaveFunctions.push_back(temp);
    }
    
    double dx = xmax/double(nx - 1);
    double sum = 0;
    
    for(int i = 0; i < numlevels; ++i){
        for(int j = i; j < numlevels; ++j){
            sum = 0;
            for(int k = 0; k < nx; ++k){
                sum += WaveFunctions[i][k]*WaveFunctions[i][k]*WaveFunctions[j][k]*WaveFunctions[j][k]*dx/(k*k*dx*dx+1.e-13);
            }
            
            AngularPart A(0, 0, 0, 0); //only doing l = 0
            double Ang_answer = 2*(3.141592653589793238)*qsimp(A, 0, 3.141592653589793238);
            Gij[i][j] = v*sum*Ang_answer;
            Gij[j][i] = Gij[i][j];
        }
    }
    
}

double Potential::f(int i) {
   // return sqrt((1/((1 + exp((Energies[i] - a)/b))*(1 + exp((-Energies[i] - a)/b)))));
}

#endif
