/*! \file discrete.cxx
 \brief Discrete potential
 \ingroup gr_pot
 */
#ifndef __DISCRETE__CXX__
#define __DISCRETE__CXX__
#include <tav/TriDiag.cxx>
#include <tav/sort.cxx>
namespace pot{
    /**!
     \brief Function to discretize potential
     Note note every part has inverse mass.
     solve equation -imass*u''+V(x)u=Eu
     If there n-variables, there are n+1 intervals and n+2 actual points. Two edges are determined by the boundary condition.
     For example with Dirichlet (zero) last points u(-1)=0 and u(n)=0 so we have n non-zero variables and these points correspond to xmin and xmax
     thus h=(xmax-xmin)/(nx+1)
     The difference in boundary condition is the value of the first diagonal, it is -2 for all in Dirichlet(value) case and -1 in Neumann(derivative)
     Another difference is that for derivative we set value between points so for Dirichlet x[0]=h, but for Neumann x[0]=h/2
     */
    template <typename funct>
    int DiscretizePotential(
                            int nx, ///<number of discrete points to be
                            double *diag, ///< [o] diagonal values
                            double *offd, ///< [o] off-diagonal
                            double xmin, ///<left side of box
                            double xmax, ///<right side of box
                            funct &VV, ///<potential itself
                            bool bleft=false, ///<left boundary condition (false- u[0]=0; true - u'[0]=0)
                            bool bright=false, ///<right boundary condition
                            double invmass=1.0 ///< inverse mass imass=hbar^2/2m
    ) {
        //
        //cerr<<"run discritizer values\n";
        //double *diag=new double [nx];
        double mpx=nx+1.0; //number of step-long intervals
        if (bleft) mpx-=0.5;
        if (bright) mpx-=0.5; //derivative is determined in the middle 1/2
        double h=(xmax-xmin)/(mpx);
        //  cerr<<xmax<<" "<<xmin<<" step="<<h<<endl;
        double h2=h*h;
        double od=-invmass/(h*h); //-1.0/h2
        double dd=-od*2.0; // 2.0/h2
        double x; //position
        if (bleft) x=xmin+0.5*h; else x=xmin+h; //starting position
        for (int ii=0;ii<nx;ii++) {
            offd[ii]=od;
            diag[ii]=dd+VV(x);
            //  cerr<<x<<" "<<diag[ii]<<" "<<offd[ii]<<endl;
            x+=h ;  //  advance step
        }
        if (bleft) diag[0]+=od;
        if (bright) diag[nx-1]+=od;
        /*
         cerr<<"output discrete values\n";
         for (int i=0;i<nx;i++) cerr<<diag[i]<<" "<<offd[i]<<endl;
         cerr<<"end output discrete \n";
         */
        return 0;
    }
    
    /*! Function to diagonalize potential
     */
    template<class cmatrix, class cvector, class funct>
    int DiscreteSolutionT(int nx, cmatrix &a, cvector &diag, double xmin, double xmax, funct &VV,
                          bool bleft=false, ///left boundary condition (false- u[0]=0; true - u'[0]=0)
                          bool bright=false, ///right boundary condition
                          double invmass=1.0 ///< inverse mass imass=hbar^2/2m
    ) {
        //first and last points u(0)=0 and u(n)=0 so we have n-1 variables
        //double *diag=new double [nx];
        //typedef typeof(a[0][0]) cnumber;
        typedef double cnumber;
        cnumber *offd=new cnumber [nx];
        DiscretizePotential(nx,diag,offd,xmin,xmax,VV,bleft,bright,invmass);
        tav::TqliT(nx,diag,offd,a);
        // delete [] diag;
        delete [] offd;
        return 0;
    }
    
    template<class cmatrix, class cvector, class funct>
    int DiscreteSolutionEivenvalues(int nx, cvector &diag, double xmin, double xmax, funct &VV,
                                    bool bleft=false, ///left boundary condition (false- u[0]=0; true - u'[0]=0)
                                    bool bright=false, ///right boundary condition
                                    double invmass=1.0 ///< inverse mass imass=hbar^2/2m
    ) {
        //first and last points u(0)=0 and u(n)=0 so we have n-1 variables
        //double *diag=new double [nx];
        //typedef typeof(diag[0]) cnumber;
        typedef double cnumber;
        cnumber *offd=new cnumber [nx];
        DiscretizePotential(nx,diag,offd,xmin,xmax,VV,bleft,bright,invmass);
        tav::TqliT(nx,diag,offd);
        // delete [] diag;
        delete [] offd;
        tav::QuickSort(nx, diag);
        return 0;
    }
    
    
    /*full packaged routine to compute eigensystem a[nx][nx+5] a[i][-]-eigenvalue a[i][nx+1]=xmin
     Storage a[i][0]=energy, a[i][1]=boundary a[i][2]...a[i][nx+1] -variable vector a[i][nx+2]-right boundary
     */
    template<class cmatrix, class funct>
    int DiscreteSolutionEigenSystem(
                                    int nx, ///<[i] dimensionality of eigensystem
                                    cmatrix &a, ///<output matrix of eigenstates, a[i][0]=eigenvalue a[i][2...nx+1]-normalized vector a[nx+2]-boundary a[nx+3]=xmin, a[nx+4]=xmax
                                    double xmin, ///<[i] xmin
                                    double xmax, ///< [i] xmax
                                    funct &VV, ///< potential
                                    bool bleft=false, ///left boundary condition (false- u[0]=0; true - u'[0]=0)
                                    bool bright=false, ///right boundary condition
                                    double invmass=1.0 ///< inverse mass imass=hbar^2/2m
    ) {
        //first and last points u(0)=0 and u(n)=0 so we have n-1 variables
        //double *diag=new double [nx];
        //typedef typeof(a[0][0]) cnumber;
        typedef double cnumber;
        cnumber *diag=new cnumber [nx];
        cnumber **vf=new cnumber* [nx];
        for (int i=0;i<nx;i++) {vf[i]=new cnumber [nx];
            for (int j=0;j<nx;j++) vf[i][j]=0.0;
            vf[i][i]=1.0;
        }
        DiscreteSolutionT(nx,vf, diag, xmin, xmax, VV,bleft,bright,invmass);
        int *index=new int [nx];
        for (int i=0;i<nx;i++) index[i]=i;
        tav::QuickSort(nx, diag, index);
        //normalization coefficient is to be computed so that function is propery normalized between given xmin and xmax
        //Normalization is always 1/\sqrt(h) because edges add as 1/2
        //DETERMINE step
        double mpx=nx+1.0; //number of step-long intervals
        if (bleft) mpx-=0.5;
        if (bright) mpx-=0.5; //derivative is determined in the middle 1/2
        double h=(xmax-xmin)/(mpx);
        
        double sqh=sqrt(h);
        for (int i=0;i<nx;i++) {
            a[i][0]=diag[i]; //save eigenvalue
            
            //Setting up the wave funciton
            double *v=a[i]+2; //for convenience set pointer to starting of the wf
            //copy vector
            for (int j=0;j<nx;j++) v[j]=vf[index[i]][j]/sqh;
            if (bleft) a[i][1]=a[i][2]; else a[i][1]=0.0; //set boundary condition
            if (bright) a[i][nx+2]=a[i][nx+1]; else a[i][nx+2]=0.0;
            //note if derivative boundary was used left/right points are shifted
            if (bleft) a[i][nx+3]=xmin-0.5*h; else a[i][nx+3]=xmin;
            if (bright) a[i][nx+4]=xmax+0.5*h; else a[i][nx+4]=xmax;
            delete [] vf[i];
        }
        delete [] diag;
        delete [] vf;
        delete [] index;
        return 0;
    }
    /*!\brief return wave function
     return value of the wf at position x, use linear extrapolation
     Important input is number of points nx, note that number of eigenvalues ne=nx-2
     number of intervals is nx-1
     storage format xmin=wf[nx], xmax=wf[nx+1]
     note that x[n]=(n+1)h+xmin
     */
    double DiscreteWaveFunction(double x, int nx, double *wf) {
        //double h=(wf[nx+1]-wf[nx])/(nx-1.0);
        if ((x<wf[nx])||(x>wf[nx+1])) CriticalError("x="<<x<<" is outside of range ["<<wf[nx]<<","<<wf[nx+1]<<"]\n");
        double jd=(nx-1.0)*(x-wf[nx])/(wf[nx+1]-wf[nx]);
        int j=int(jd);
        //if ((j<=0)||(j>(nx))) return 0.0;
        jd-=double(j);
        return (1.0-jd)*wf[j]+jd*wf[j+1];
        //if (j==(nx-1)) return jd*wf[j]; //this should never happen bu just in case?
        //else return jd*wf[j]+(1.0-jd)*wf[j+1];
    }
    
    /*full packaged routine to compute eigensystem a[nx][nx+3] a[i][-]-eigenvalue a[i][nx+1]=xmin
     a[2nx][2nx+6]  32123
     we add extra point to symmetric solution
     */
    
    template<class cmatrix, class funct>
    int DiscreteSymmetricSolutionEigenSystem(
                                             int n2, ///<[i] dimensionality of eigensystem
                                             cmatrix &a, ///<output matrix of eigenstates, a[i][0]=eigenvalue a[i][1...nx]-normalized vector a[nx+1]=xmin, a[nx+2]=xmax
                                             double xmax, ///< [i] xmax
                                             funct &VV, ///< potential
                                             double invmass=1.0 ///< inverse mass imass=hbar^2/2m
    ) {
        int nx=n2/2;
        double xmin=0.0; //potential starts from zero
        //  tav::Tensor<2,double> wf(nx,nx+5);
        double **wf =new double* [nx];
        for (int i=0;i<nx;i++) wf[i]=new double [nx+5];
        
        double sqrt2=sqrt(2.0);
        DiscreteSolutionEigenSystem(nx,wf,xmin, xmax, VV,true,false,invmass); //with zero boundary condition
        //process storage of symmetric system
        for (int i=0;i<nx;i++) {
            int iv=2*i;
            a[iv][0]=wf[i][0]; //save energy
            //if (i<10) cerr<<wf[i][0]<<endl;
            //nx+2 points in wf, with symmetry this is 2nx+3 points (one overlap at zero)
            for (int j=1;j<=nx+2;j++) a[iv][j]=wf[i][nx+3-j]/sqrt2; //finish with a[][nx+2]
            for (int j=1;j<=nx+2;j++) a[iv][nx+j]=wf[i][j]/sqrt2; //correctly one should do this !!! no triple overlap at zero
            // for (int j=1;j<=nx+2;j++) a[iv][nx+j+1]=wf[i][j]; //note j=0 here nx+2 coinsides with above and is zero last 2nx+3
            // a[iv][2*nx+3]=0.0; //fill extra point with zero
            double h=xmax/(nx+0.5); //step
            a[iv][2*nx+4]=-xmax; 
            //double h=xmax/(nx+0.5); //step
            a[iv][2*nx+5]= xmax+h; //because we add extra blank point
            //detete  if we do not use tensor
            for (int i=0;i<nx;i++) delete [] wf[i];
            delete [] wf;
        }
        //asymmetric   
        DiscreteSolutionEigenSystem(nx,wf,xmin, xmax, VV,false,false,invmass); //with zero boundary condition
        //process storage of symmetric system
        //for (int jj=0;jj<nx+2;jj++) cerr<<wf[0][jj+1]<<endl;
        for (int i=0;i<nx;i++) {
            int iv=2*i+1;
            a[iv][0]=wf[i][0]; //save energy
            //if (i<10) cerr<<wf[i][0]<<endl;
            //nx+2 points in wf, with symmetry this is 2nx+3 points (one overlap at zero)
            for (int j=1;j<=nx+2;j++) a[iv][j]=-wf[i][nx+3-j]/sqrt2; //finish with a[][nx+2]
            for (int j=1;j<=nx+2;j++) a[iv][nx+j+1]=wf[i][j]/sqrt2; //note j=0 here nx+2 coinsides with above and is zero last 2nx+3
            a[iv][2*nx+4]=-xmax;
            a[iv][2*nx+5]= xmax;
        }
        
        
        return 0;
    }
    
    
    /*
     double PotentialInterpolationWF(int n, double x) {
     double step=(xmax-xmin)/(nx+1.);
     int j=int((x-xmin)/step);
     if ((j<0)||(j>(nx-1))) return 0.0;
     //return wf[n][j+1];
     double y, dy;
     int imin=0; 
     if (j>2) imin=j-2;
     if (j>nx-2) imin=nx-4;
     //cout<<"j= "<<j<<" jmin="<<imin<<endl;
     //for (int q=0;q<4;q++) cout<<(&(wf[n][imin+1]))[q]<<" "<<wf[n][imin+q+1]<<endl;
     //recall that array starts from 1
     return tav::EquidistantPolynomialInterpolation(4,xmin+imin*step,step,(&(wf[n][imin+1])), x, y,dy);
     //return wf[n][j+1];
     }
     */
    
} //end of namespace pot


#endif //__DISCRETE__CXX__
