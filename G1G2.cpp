//
// Computes G'(omega), G"(omega), tandelta(omega), eta'(omega) and eta"(omega)
//          from G(t) (stress autocorrelation function, shear relaxation modulus).
//
// If you publish work based on this code, please cite:
// David, A., De Nicola, A., Tartaglino, U., Milano, G. and Raos, G., 2019.
// Viscoelasticity of Short Polymer Liquids from Atomistic Simulations.
// Journal of The Electrochemical Society, 166(9), p.B3246.
// https://iopscience.iop.org/article/10.1149/2.0371909jes
//
// To compile: g++ -O G1G2.cpp -o G1G2.exe
// To execute: ./G1G2.exe OMEGA_MIN OMEGA_MAX NT_MIN NT_MAX DELTA_T > OUTPUT.dat
//             where:
//             OMEGA_MIN and OMEGA_MAX specify the range of frequencies (in rad/ns),
//             NT_MIN and NT_MAX specify the time range for averaging the integrals (no. of points),
//             DELTA_T is the timestep for the input G(t) (in ns).
//
// Other important "input" parameters that are not specified on the command line.
// The filenames with the components of G_ij(t) (input):
char filexx[]="pxx_subtrace.in.acf", \
     filexy[]="pxy.in.acf", \
     filexz[]="pxz.in.acf", \
     fileyy[]="pyy_subtrace.in.acf", \
     fileyz[]="pyz.in.acf", \
     filezz[]="pzz_subtrace.in.acf";
// Max. number of G(t) points to be read from files.
const int nmax=3000000;
double vec[nmax], Gt[nmax];
// Max. number of points to be computed in the frequency domain
const int nmaxo=1000;
double G1[nmaxo], G2[nmaxo], eta1[nmaxo], eta2[nmaxo];
// The filename for the output data:
//char fileout[]="G1G2.dat";

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string.h>

using namespace std;

// Writes one vector
void writevec(double *vec, int n, char *filename){
   ofstream outfile(filename);
   for (int i=0; i<n; i++) {
//      outfile << i << " " << vec[i] << endl;
      outfile << vec[i] << endl;
   }
   outfile.close();
}


// Writes three vectors on three parallel columns
void writevec3(double *vec0, double *vec1, double *vec2, int n, char *filename){
   ofstream outfile(filename);
   for (int i=0; i<n; i++) {
      outfile << i << " " << vec0[i] << " " << vec1[i] << " " << vec2[i] << endl;
   }
   outfile.close();
}


// Reads one vector
int readvec(double *vec, char *filename){
   string line;
   double t;
   ifstream infile(filename);
   if (!infile) {
      cout << "Cannot open " << filename << endl;
      return -1;
   }
   int n=0;
   while (getline(infile, line)) n++;
//   if (n>nmax) {
//      cout << "Too many data: n=" << n << " > " << nmax << endl;
//      return -2;
//   }
   infile.clear();
   infile.seekg(0);
   n = min(n,nmax);
   for (int i=0; i<n; i++) {
      getline(infile, line);
      istringstream iss(line);
//      iss >> t >> vec[i];
      iss >> vec[i];
   }
   infile.close();
   return n;
}


// Sums vec2 to vec1; vec2 scan be scaled by f before the summation.
void sumvec(double *vec1, double *vec2, int n, double f=1.){
   for (int i=0; i<n; i++) vec1[i] += f*vec2[i];
}


// This computes the average of a vector.
double avrg(double *vec, int n){
   int i;
   double a=0.0;
   for (i=0; i<n; a+=vec[i], i++);
   a /= n;
   return a;
}


// This computes the variance of a vector.
// It also modifies VEC by subtracting its AVRG!
double sigma2(double *vec, int n){
   int i;
   double s2=0.0;
   double a = avrg(vec,n);
   for (i=0; i<n; i++){
      vec[i] -= a;
      s2 += vec[i]*vec[i];
   }
   s2 /= n;
   return s2;
}


// The MAIN program
int main(int argc, char *argv[]) {
   int i, j, n, n0, n1, n2;
   double omega, t, temp1, temp2;
   char outfile[80];

   if (argc != 6) {
      cout << "Input error: enter OMEGA_MIN OMEGA_MAX NT_MIN NT_MAX DELTA_T on command line" << endl;
      return argc;
   }
   double omin = atof(argv[1]);
   double omax = atof(argv[2]);
   int ntmin = atoi(argv[3]);
   int ntmax = atoi(argv[4]);
   double dt = atof(argv[5]);
    
// fomega determines the "step" (on log scale) for the spacing of frequencies:
   double fomega = pow(omax/omin,1./(nmaxo-1.));
    
// Read and accumulate data with autocorrelation functions for all stress components.    
// Diagonal components (weight=1)
   n = readvec(vec,filexx);
   cout << n << " data have been read for XX from " << filexx << endl;
   n0 = n;
   sumvec(Gt,vec,n);

   n = readvec(vec,fileyy);
   cout << n << " data have been read for YY from " << fileyy << endl;
   if (n!=n0) cout << "ERROR ON NO. OF DATA!" << endl;
   sumvec(Gt,vec,n);

   n = readvec(vec,filezz);
   cout << n << " data have been read for ZZ from " << filezz << endl;
   if (n!=n0) cout << "ERROR ON NO. OF DATA!" << endl;
   sumvec(Gt,vec,n);

// Off-diagonal components (weight=2)
   n = readvec(vec,filexy);
   cout << n << " data have been read for XY from " << filexy << endl;
   if (n!=n0) cout << "ERROR ON NO. OF DATA!" << endl;
   sumvec(Gt,vec,n,2.);

   n = readvec(vec,filexz);
   cout << n << " data have been read for XZ from " << filexz << endl;
   if (n!=n0) cout << "ERROR ON NO. OF DATA!" << endl;
   sumvec(Gt,vec,n,2.);

   n = readvec(vec,fileyz);
   cout << n << " data have been read for YZ from " << fileyz << endl;
   if (n!=n0) cout << "ERROR ON NO. OF DATA!" << endl;
   sumvec(Gt,vec,n,2.);

// Normalize the sum
   for (i=0; i<n0; i++) {
       Gt[i] /= 10.;
   }

   n=0;
// The real calculation starts here
   n1 = min(n0,ntmin);
   n2 = min(n0,ntmax);
   cout << "Time integration from " << n1 << " (points) x " << \
        dt << " (ns/point) = " << n1*dt <<" ns" << endl;
   cout << "                up to " << n2 << " (points) x " << \
        dt << " (ns/point) = " << n2*dt <<" ns" << endl;
   cout << endl;
   cout << "#  omega        G1         G2       tandelta       eta1         eta2 " << endl;
// Loop over all frequencies  
   omega = omin;
   for (i=0; i<nmaxo; i++){
//     Loop over all times and integrate by trapezoidal rule.
//     Average the integrals from n1 to n2.
       temp1 = Gt[0]/2;
       temp2 = 0.;
       eta1[i] = 0.;
       eta2[i] = 0.;
       for (j=1; j<n2; j++) {
           t = j*dt;
           temp1 += Gt[j]*cos(omega*t);
           temp2 += Gt[j]*sin(omega*t);
           if (j>=n1) {
               eta1[i] += temp1;
               eta2[i] += temp2;
           }
       }
       eta1[i] *= dt/(n2-n1);
       eta2[i] *= dt/(n2-n1);
       G1[i] = eta2[i]*omega;
       G2[i] = eta1[i]*omega;
       cout << omega <<" "<< G1[i] <<" "<< G2[i] <<" "<< \
               eta1[i]/eta2[i] <<" "<< eta1[i] <<" "<< eta2[i] << endl;
//     Change omega at the end of the cycle
       omega *= fomega;
   }
  
   return 0;
}
