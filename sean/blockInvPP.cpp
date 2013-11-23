#include "Matrix.h"
#include "Barrier.h"
#include <math.h>
#include <stdio.h>
#include <thread>
#include <chrono>

using namespace std;

// return n square empty matrices by pointer
vector<Matrix> *allocateNMatrices(size_t n) {
    vector<Matrix> *ret = new vector<Matrix>;
    for (size_t i=0; i<n; i++) {
        ret->push_back(Matrix(0,0));
    }
    return ret;
}

// performs forward substitution on one column a matrix
Matrix forwardSubCol(Matrix &coefs, vector<double> &eqls) {
	Matrix xs = Matrix(eqls.size(), 1);
	for (size_t y=0; y<coefs.rows.size(); y++) {
		double sum = 0;
		for (size_t x=0; x<y; x++) {
            // this is multiplying the 
			sum += coefs.rows[y][x] * xs.rows[x][0];
		}
		xs.rows[y][0] = eqls[y] - sum / coefs.rows[y][y];
	}
	return xs;
}

// performs forward substitution on a matrix
Matrix forwardSub(Matrix &coefs, Matrix &eqls) {
	Matrix xs = Matrix(coefs.rows.size(), 0);
	for (size_t x=0; x<eqls.rows[0].size(); x++) {
        // get xth column and perform forward substitution on that column
		vector<double> eqlsCol = eqls.sliceCol(x);
		Matrix xsCol = forwardSubCol(coefs, eqlsCol);
		xs.appendCol(xsCol);
	}
	return xs;
}

vector<Matrix> sliceLs(Matrix coefs, size_t blockSize) {
	vector<Matrix> ls;
	for (size_t i=0; i<(size_t)coefs.rows.size() / blockSize; i++) {
		ls.push_back(coefs.sliceBlock(i * blockSize, i * blockSize, blockSize, blockSize));
	}
	return ls;
}

vector<Matrix> sliceRs(Matrix coefs, size_t blockSize) {
	vector<Matrix> rs;
	for (size_t i=1; i<(size_t)coefs.rows.size() / blockSize; i++) {
		rs.push_back(coefs.sliceBlock(i * blockSize, (i - 1) * blockSize, blockSize, blockSize));
	}
	return rs;
}

vector<Matrix> calcGs(vector<Matrix> &rs, vector<Matrix> &ls) {
	vector<Matrix> gs;
	for (size_t i=0; i<rs.size(); i++) {
		gs.push_back(forwardSub(ls[i+1], rs[1]));
	}
	return gs;

}

vector<Matrix> calcBis(vector<Matrix> &ls, vector<Matrix> &ansBlocks) {
	vector<Matrix> bis;
	for (size_t i=0; i<ls.size(); i++) {
		bis.push_back(forwardSub(ls[i], ansBlocks[i]));
	}
	return bis;
}

vector<Matrix> makeGHats(vector<Matrix> &gs, size_t blockSize, size_t bandwidth) {
	vector<Matrix> gHats;
	size_t firstZero = blockSize - bandwidth;
	for (size_t i=0; i<gs.size(); i++) {
		Matrix gHat = Matrix(blockSize, bandwidth);
		for (size_t y=0; y<blockSize; y++) {
			for (size_t x=0; x<bandwidth; x++) {
				gHat.rows[y][x] = gs[i].rows[y][firstZero + x];
			}

		}
		gHats.push_back(gHat);
	}
	return gHats;
}

vector<SplitMatrix> splitByBand(vector<Matrix> &mtx, size_t blockSize, size_t bandwidth) {
	vector<SplitMatrix> split;
	for (size_t i=0; i<mtx.size(); i++) {
        Matrix top = mtx[i].sliceRows(0, blockSize - bandwidth);
        Matrix bot = mtx[i].sliceRows(blockSize - bandwidth, blockSize);
		split.push_back(SplitMatrix(top, bot)); 
	}
	return split;
}

vector<Matrix> solveZsBlockInv(vector<SplitMatrix> &UVs, vector<SplitMatrix> &MHs) {
	vector<Matrix> zs;
	zs.push_back(UVs[0].bottom);
	for (size_t i=1; i<UVs.size(); i++) {
		zs.push_back(UVs[i].bottom - MHs[i-1].bottom * zs[i-1]);
	}
	return zs;
}

vector<Matrix> assemblePrefixComponents(vector<SplitMatrix> &MHs, vector<SplitMatrix> &UVs) {
	vector<Matrix> comps;
	Matrix first = Matrix(UVs[0].bottom.rows.size() + 1, MHs[0].bottom.rows[0].size() + UVs[0].bottom.rows[0].size());
	first.populateCol(first.rows[0].size() - 1, 1);
	first.pasteIn(UVs[0].bottom, 0, MHs[0].bottom.rows[0].size());
	comps.push_back(first);
	for (size_t i=0; i<MHs.size(); i++) {
		Matrix compNew = Matrix(UVs[0].bottom.rows.size() + 1, MHs[0].bottom.rows[0].size() + UVs[0].bottom.rows[0].size());
		compNew.populateCol(compNew.rows[0].size() - 1, 1);
		compNew.pasteIn(MHs[i].bottom * -1, 0, 0);
		compNew.pasteIn(UVs[i + 1].bottom, 0, MHs[0].bottom.rows[0].size());
		comps.push_back(compNew);
	}
	return comps;
}

void recursiveSolvePrefix(vector<Matrix> &xs) {

    size_t n = xs.size();
    
    for (size_t stepSize=2; stepSize<=n; stepSize*=2) {
        size_t start = stepSize/2 - 1;
        size_t lookForward = stepSize/2;
        for (size_t j=start; j<n; j+=stepSize) {
            xs[j+lookForward] = xs[j+lookForward] * xs[j];
        }
    }

    for (size_t stepSize=n/2+0.5; stepSize>=2; stepSize/=2) {
        size_t start = stepSize-1;
        size_t lookForward = stepSize/2;
        for (size_t i=start; i<n-1; i+=stepSize) {
            xs[i+lookForward] = xs[i+lookForward] * xs[i];
        }
    }
    
}

// performs the prefix sum algorithm in parallel
// &xs: reference to a vector of matrices that will store the solutions
//      done in-place with shared memory
// &b: a barrier object that allows for the processors to painlessly synchronize
// tid: the thread id of the processor
// P: the number of processors
//
// i'm not 100% happy with the style (passing in b and P is not ideal)
// but for now this works!
void recursiveSolvePrefix_par(vector<Matrix> &xs, ThBarrier &b, size_t tid, size_t P) 
{
    size_t n = xs.size();
    
    // number of matrices in the xs vector delegated to each processor
    size_t chunk = n/P;
    
    // descent: multiply adjacent pairs, always storing to the right
    for (size_t stepSize=2; stepSize<=n; stepSize*=2) {
        // only take the processors that are necessary
        //
        // if there are only 4 matrices of interest and 8 processors,
        // take the first 4 processors and scale them to the correct location
        if (chunk*tid % stepSize == 0) {
            size_t start = stepSize/2 - 1 + chunk*tid;
            size_t lookForward = stepSize/2;
            // perform all multiplications within this processor's chunk
            for (size_t j=start; j<start+chunk; j+=stepSize) {
                xs[j+lookForward] = xs[j+lookForward] * xs[j];
            }
        }
        // wait for all threads to complete this iteration before continuing
        b.sync();
    }

    // ascent: multiply adjacent pairs, always storing to the right
    // this time, though, the step is effectively half what it was before
    for (size_t stepSize=n/2+0.5; stepSize>=2; stepSize/=2) {
        // again, only take the processors we actually need...
        if (chunk*tid % stepSize == 0) {
            size_t start = stepSize-1 + chunk*tid;
            size_t lookForward = stepSize/2;
            for (size_t i=start; i<start+chunk && i<n-1; i+=stepSize) {
                xs[i+lookForward] = xs[i+lookForward] * xs[i];
            }
        }
        b.sync();
    }
    
}

void recursiveSolvePrefix_threads(vector<Matrix> &continComponents) {
    size_t P = 4;
    size_t tid = 0;
    vector<thread> threads(P);
    ThBarrier b(P);
    while (tid < P) {
        threads[tid] = thread(recursiveSolvePrefix_par, std::ref(continComponents), std::ref(b), tid, P);
        ++tid;
    }
    tid = 0;
    while (tid < P) {
        threads[tid].join();
        ++tid;
    }
}

vector<Matrix> solveZsPrefix(vector<SplitMatrix> MHs, vector<SplitMatrix> UVs) {
	vector<Matrix> continComponents = assemblePrefixComponents(MHs, UVs);
    for (size_t i=0; i<continComponents.size(); i++) {
        //printf("%lu: %f %f \n", i, continComponents[i].rows[0][2], continComponents[i].rows[1][2]);
        continComponents[i].print();
    }



	recursiveSolvePrefix_threads(continComponents);
	//recursiveSolvePrefix(continComponents);
    for (size_t i=0; i<continComponents.size(); i++) {
        //printf("%lu: %f %f \n", i, continComponents[i].rows[0][2], continComponents[i].rows[1][2]);
        continComponents[i].print();
    }

	vector<Matrix> hProds;
	Matrix first = Matrix(UVs[0].bottom.rows.size() + 1, MHs[0].bottom.rows[0].size() + UVs[0].bottom.rows[0].size());
	first.populateCol(first.rows[0].size() - 1, 1);
	first.pasteIn(UVs[0].bottom, 0, MHs[0].bottom.rows[0].size());
	hProds.push_back(first);
	for (size_t i=0; i<MHs.size(); i++) {
		Matrix hNew = Matrix(UVs[0].bottom.rows.size() + 1, MHs[0].bottom.rows[0].size() + UVs[0].bottom.rows[0].size());
		hNew.populateCol(hNew.rows[0].size() - 1, 1);
		hNew.pasteIn(MHs[i].bottom * -1, 0, 0);
		hNew.pasteIn(UVs[i + 1].bottom, 0, MHs[0].bottom.rows[0].size());
		for (size_t j=hProds.size() - 1; j>=0; j--) {
			hNew = hNew * hProds[j];
		}
		hProds.push_back(hNew);
	}
	vector<Matrix> foo;
	return foo;
}

vector<Matrix> solveYs(vector<SplitMatrix> &UVs, vector<SplitMatrix> &MHs, vector<Matrix> &zs) {
	vector<Matrix> ys;
	ys.push_back(UVs[0].top);
	for (size_t i=1; i<UVs.size(); i++) {
		//ys.push_back(UVs[i].top - MHs[i-1].top * zs[i-1]);
	}
	return ys;
}

Matrix solveXs(vector<Matrix> &ls, vector<Matrix> &bis, vector<Matrix> &ans, vector<Matrix> &gs, size_t bandwidth) {
	size_t blockSize = ls[0].rows.size();
	vector<Matrix> gHats = makeGHats(gs, blockSize, bandwidth);
	vector<SplitMatrix> MHs = splitByBand(gHats, blockSize, bandwidth);
	vector<SplitMatrix> UVs = splitByBand(bis, blockSize, bandwidth);
	vector<Matrix> zs = solveZsPrefix(MHs, UVs);
    /*
	//vector<Matrix> zs = solveZsBlockInv(UVs, MHs);
	vector<Matrix> ys = solveYs(UVs, MHs, zs);
	Matrix xs = Matrix(zs.size() * zs[0].rows.size() + ys.size() * ys[0].rows.size(), 1);
	size_t index = 0;
	for (size_t i=0; i<zs.size(); i++) {
		Matrix *yGroup = &ys[i];
		Matrix *zGroup = &zs[i];

		for (size_t row=0; row<yGroup->rows.size(); row++) {
			xs.rows[index][0] = yGroup->rows[row][0];
			index++;
		}
		for (size_t row=0; row<zGroup->rows.size(); row++) {
			xs.rows[index][0] = zGroup->rows[row][0];
			index++;
		}
	}
	return xs;
    */
	return gHats[0];
}

int main(int argc, char *argv[])
{
    // numProc MUST be a power of 2 :/ or else it doesn't work...
	size_t numProc = 8;
	size_t blockSize = 5;
	size_t mtxSize = numProc * blockSize;
	Matrix coefs = Matrix(mtxSize, mtxSize);
    coefs.populateDiagonal(0, 0, 1);
	coefs.populateDiagonal(1, 0, -1);
    //coefs.print();
	//coefs.populateDiagonal(2, 0, -2);
	vector<Matrix> rs = sliceRs(coefs, blockSize);
	vector<Matrix> ls = sliceLs(coefs, blockSize);
	vector<Matrix> gs = calcGs(rs, ls);
	size_t bandwidth = 2;
	Matrix ans = Matrix(mtxSize, 1);
	ans.populateCol(0, 1);
	vector<Matrix> ansBlocks = ans.asRowBlocks(blockSize);
	vector<Matrix> bis = calcBis(ls, ansBlocks);
	Matrix xs = solveXs(ls, bis, ansBlocks, gs, bandwidth);
	return 0;
}

