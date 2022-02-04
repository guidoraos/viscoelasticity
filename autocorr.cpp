//
// Computes autocorrelation functions by recursive coarse-graining in time.
//
// If you publish work based on this code, please cite:
// David, A., De Nicola, A., Tartaglino, U., Milano, G. and Raos, G., 2019.
// Viscoelasticity of Short Polymer Liquids from Atomistic Simulations.
// Journal of The Electrochemical Society, 166(9), p.B3246.
// https://iopscience.iop.org/article/10.1149/2.0371909jes 
//
// To compile:    g++ autocorr.cpp -o autocorr
// To execute:    ./autocorr FILENAME NBLOCK
//                where NBLOCK determines the size of block averages;
//                if NBLOCK is omitted, it is assumed NBLOCK=1 (i.e., no blocking).
//                The autocorrelations functions are writen to FILENAME.acf
//

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string.h>

using namespace std;

// The NBLOCK parameter (from the input) determines the amount of
// coarse-graining.
int nblock;

// The FACTOR parameter (below) determines the delay in the application
// of coarse-graining. A higher value of FACTOR implies a more accurate
// (and expensive) calculation.
const int factor = 100;

// This vector size should be enough for most practical purposes.
const int nmax=100000000;
double vec[nmax], vacf[nmax];
//double vacf0[nmax], vacf1[nmax];


// This performs block averaging
int blockave(double *vec, int n){
   int i, j, m=0;
   double t;
   for (i=0; i<n-nblock+1; i+=nblock){
      t = 0.0;
      for (j=i; j<i+nblock; j++)
         t += vec[j];
      vec[m] = t/nblock;
      m++;
   }
   return m;
}


// This performs linear interpolation
void interpol(double *vec, int m, double v1, double v2){
   double f;
   for (int i=0; i<m; i++){
      f = ((double)(m-i)) / ((double)(m));
      vec[i] = f*v1 + (1.0-f)*v2;
   }
}


// Another linear interpolation
void interpol2(double *vec, int m, double v1, double v2, double dl){
   double f;
   for (int i=0; i<m; i++){
      f = (dl-i) / dl;
      vec[i] = f*v1 + (1.0-f)*v2;
   }
}


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
   if (n>nmax) {
      cout << "Too many data: n=" << n << " > " << nmax << endl;
      return -2;
   }
   infile.clear();
   infile.seekg(0);
   for (int i=0; i<n; i++) {
      getline(infile, line);
      istringstream iss(line);
//      iss >> t >> vec[i];
      iss >> vec[i];
   }
   infile.close();
   return n;
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


// This computes the exact ACF (cost is O(N^2))
void acf0(double *vec, int n, double *acf){
   int i, j;
   for (i=0; i<n; acf[i]=0.0, i++);
   cout << endl << "ACF0 running, now at i=... " << endl;
   for (i=0; i<n; i++){
      if (i%1000==0) cout << i << endl;
      for (j=i; j<n; j++)
         acf[j-i] += vec[i]*vec[j];
   }
   for (i=0; i<n; i++){
      acf[i] /= (n-i);
   }
}


// This computes ACF1 (coarse-grained) from ACF0 (exact)
void acf0acf1(double *acf0, int n, double *acf1){
   int i, j, k, ngen, mm, m2;
   double f1, f2;
   int delta = factor*nblock;
   for (i=0; i<delta; i++) acf1[i]=acf0[i];
   f2 = acf0[delta-1];
   m2 = 1;
   k = delta;
   j = delta-1;
   cout << endl;
   for (ngen=1, mm=nblock; k<n ; ngen++, mm*=nblock){
      cout << "ACF0ACF1 at generation no. " << ngen;
      cout << ", block size " << mm;
      cout << ", now at position " << k << endl;
      for (i=0; i<delta; i++){
         f1 = f2;
         f2 = avrg(acf0+k, mm);
         interpol(acf1+j, (mm+m2)/2, f1, f2);
         k += mm;
         j += (mm+m2)/2;
         m2 = mm;
         if (j>n) break;
      }
   }
}


// This computes the coarse-grained ACF (cost is O(N)).
void acf1(double *vec, int n, double *acf){
   int i, j, kk, mm, ngen, mt;
   int jstart, jstop, imax;
   double temp[10000], ttemp;   // delta should not exceed 10000
   int count[10000];
   int nn=n;
   int delta=factor*nblock;
   for (i=0; i<n; acf[i]=0.0, i++);
   cout << endl;
   for (ngen=0, mm=1, kk=0, jstop=0; kk<nn; ngen++, mm*=nblock){
      cout << "ACF at generation no. " << ngen;
      cout << ", block size " << mm;
      cout << ", now at position " << kk << endl;
      for (i=0; i<delta; i++){
         temp[i] = 0.0;
         count[i] = 0;
      }
      jstart = jstop/nblock;
      jstop = jstart+delta;
      for (i=0; i<n; i++){
         for (j=i+jstart; j<i+jstop; j++) {
            if (j==n) break;
            temp[j-i-jstart] += vec[i]*vec[j];
            count[j-i-jstart] += 1;
         }
      }
      for (i=0; i<delta; i++){
         if (count[i]==0) break;
         temp[i] /= count[i];
      }
      imax = i-1;
      if (ngen==0){
         for (i=0; i<imax; i++){
//            cout << kk << " ";
            acf[kk]=temp[i];
            kk++;
         }
//         cout << endl;
         ttemp = temp[imax];
         mt = 1;
      }
      else{
//         cout << kk << " ";
         interpol(acf+kk, mt, ttemp, temp[0]);
         kk += mt;
         for (i=0; i<imax; i++){
//            cout << kk << " ";
            interpol(acf+kk, mm, temp[i], temp[i+1]);
            kk += mm;
         }
//         cout << endl;
         ttemp=temp[imax];
         mt = mm;
      }
      n = blockave(vec, n);
   }
}


// The MAIN program
int main(int argc, char *argv[]) {
   int i, n;
   char outfile[80];

   if (argc < 2 || argc > 3) {
      cout << "Input error: enter FILENAME (NBLOCK) on command line" << endl;
      return argc;
   }
   char* filename = argv[1];
   if (argc==3)
      nblock = atoi(argv[2]);
   else
      nblock = 1;
   n = readvec(vec,filename);
   cout << endl;
   cout << n << " data lines have been read." << endl;

//// Uncomment this if you want to compute the average and variance of VEC,
//// and subtract the average from the elements of VEC.
//   cout << "Average: " << avrg(vec,n) << endl;
//   cout << "Sigma2:  " << sigma2(vec,n) << endl;
  
//// Compute the coarse-grained ACF (or the exact one if NBLOCK==1).
   if (nblock>1)
      acf1(vec, n, vacf);
   else
      acf0(vec, n, vacf);
   strcpy(outfile,filename);
   strcat(outfile,".acf");
   writevec(vacf, n, outfile);

//// Uncomment to compute the exact ACF, and then coarse-grain for comparison.
//// Uncomment also  "//double vacf0[nmax], vacf1[nmax];" at the top of file.
//   n = readvec(vec,filename);
//   acf0(vec, n, vacf0);
//   acf0acf1(vacf0, n, vacf1);
//   strcpy(outfile,filename);
//   strcat(outfile,".3acf");
//   writevec3(vacf, vacf0, vacf1, n, outfile);

   return 0;
}
