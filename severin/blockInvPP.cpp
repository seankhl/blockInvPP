#include <iostream>
#include <omp.h>
#include <Eigen/SparseCore>
#include <math.h>

using namespace Eigen;
using namespace std;

int main()
{
    const int chunk = 4; // must be >= bandwidth
    const int bandwidth = 2; // must be = to the actual bandwidth of the input matrix
    int numLevels = 2; // must be log_2 (numProcs)
    const int numProcs = 4; // must be equal to the actual #procs on the platform
    typedef Matrix<double, bandwidth+1, bandwidth+1> MatrixBd;
    MatrixBd* h = new MatrixBd[numProcs];
    MatrixBd* hprdcts = new MatrixBd[numProcs];
    typedef Matrix<double, bandwidth, 1> VectorBd;
    VectorBd* z = new VectorBd[numProcs];
    typedef Matrix<double, chunk-bandwidth, 1> VectorUd;
    VectorUd* u = new VectorUd[numProcs];
    VectorUd* y = new VectorUd[numProcs];
    typedef Matrix<double, chunk-bandwidth, bandwidth> MatrixMd;
    MatrixMd* M = new MatrixMd[numProcs];
    
    
    // Parallel Phase 1:
    // Each processor initializes it's matrices and vector and
    // does the forward substitution.
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        
        if (tid == 0)
        {
            SparseMatrix<double> L(chunk,chunk);
            VectorXd b = VectorXd::Zero(chunk);
            VectorXd c = VectorXd::Zero(chunk);
            
            for (int i=0; i<chunk; i++)
            {
                L.insert(i,i) = 1;
                b(i) = 1;
            }
            
            for (int i=0; i<chunk-1; i++) {
                L.insert(i+1,i) = -1;
            }
            
            c = L.triangularView<Eigen::Lower>().solve(b);
            for (int k=0; k<bandwidth; k++) {
                h[tid](k,bandwidth) = c(chunk-bandwidth+k);
            }
            h[tid](bandwidth,bandwidth) = 1;
            hprdcts[tid] = h[tid];
            
            y[tid] = c.head(chunk-bandwidth);
            
        } else {
            
            SparseMatrix<double> L(chunk,chunk);
            MatrixXd R = MatrixXd::Zero(chunk,chunk);
            MatrixXd G = MatrixXd::Zero(chunk,chunk);
            VectorXd b = VectorXd::Zero(chunk);
            VectorXd c = VectorXd::Zero(chunk);
            
            for (int i=0; i<chunk; i++)
            {
                L.insert(i,i) = 1;
                b(i) = 1;
            }
            
            for (int i=0; i<chunk-1; i++) {
                L.insert(i+1,i) = -1;
            }
            
            R(0,chunk-1) = -1;
            
            G = L.triangularView<Eigen::Lower>().solve(R);
            c = L.triangularView<Eigen::Lower>().solve(b);
            for (int k=0; k<bandwidth; k++) {
                h[tid](k,bandwidth) = c(chunk-bandwidth+k);
                for (int l=0; l<bandwidth; l++) {
                    h[tid](k,l) = - G(chunk-bandwidth+k,chunk-bandwidth+l);
                }
            }
            h[tid](bandwidth,bandwidth) = 1;
            if (tid != numProcs-1) {
                hprdcts[tid] = h[tid];
            } else {
                hprdcts[tid] = MatrixBd::Identity();
            }
            
            u[tid] = c.head(chunk-bandwidth);
            for (int k=0; k<chunk-bandwidth; k++) {
                for (int l=0; l<bandwidth; l++) {
                    M[tid](k,l) = G(k,chunk-bandwidth+l);
                }
            }
        }
    }
    
    // Parallel Phase 2: Parallel Prefix
    // 1st : Climbing up the tree
    for (int level=1; level<numLevels; level++) {
        int iStart = pow(2,(double) level)-1;
        int iOffset = pow(2,(double) level);
        int iOffsetPrev = iOffset / 2;
#pragma omp for schedule(static,1)
        for (int pos=iStart; pos<numProcs-1; pos+=iOffset) {
            hprdcts[pos] = hprdcts[pos] * hprdcts[pos-iOffsetPrev];
        }
    }
    // 2nd : Climbing down the tree
    for (int level=numLevels; level>0; level--) {
        int iOffset = pow(2,(double) level);
        int iOffsetPrev = iOffset / 2;
#pragma omp for schedule(static,1)
        for (int pos=numProcs-1; pos>0; pos-=iOffset) {
            MatrixBd temp = hprdcts[pos-iOffsetPrev];
            hprdcts[pos-iOffsetPrev] = hprdcts[pos];
            hprdcts[pos] = temp * hprdcts[pos];
        }
    }
    // 3rd : Doing one last multiplication step
#pragma omp for schedule(static,1)
    for (int i=0; i<numProcs; i++) {
        hprdcts[i] = h[i] * hprdcts[i];
    }
    
    // Parallel Phase 3: Computing the z_i
#pragma omp for schedule(static,1)
    for (int i=0; i<numProcs; i++) {
        VectorXd vecForExtr = VectorXd::Zero(bandwidth+1);
        vecForExtr(bandwidth) = 1;
        VectorXd temp = hprdcts[i] * vecForExtr;
        z[i] = temp.head(bandwidth);
    }
    
    // Parallel Phase 4: Computing the y_i
#pragma omp for schedule(static,1)
    for (int i=1; i<numProcs; i++) {
        y[i] = u[i] - M[i] * z[i-1];
    }
    
    for (int a=0; a<numProcs; a++) {
        for (int b=0; b<chunk-bandwidth; b++) {
            printf("%1f",y[a](b));
            printf("\n");
        }
        for (int c=0; c<bandwidth; c++) {
            printf("%1f",z[a](c));
            printf("\n");
        }
    }
    
    return 0;
}

