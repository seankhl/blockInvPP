
#include <vector>
using namespace std;

class Matrix {
  public:
	vector<vector<double>> rows;
	
    Matrix(int nRow, int nCol);
	Matrix();
	void populateRow(int row, double val);
	void populateCol(int col, double val);
    void populateDiagonal(int row, int col, double val); 

	void appendCol(Matrix &col);
	//void pasteIn(Matrix &paste, int row, int col);
	void pasteIn(Matrix paste, int row, int col);
	
	Matrix sliceRow(int row);
    Matrix sliceCol(int col);
	Matrix sliceRows(int min, int max); //up to but not including b
	Matrix sliceBlock(int row, int col, int numRows, int numCols);
	vector<Matrix> asRowBlocks(int numBlocks);

	Matrix operator+ (Matrix &m);
	Matrix operator- (Matrix &m);
	Matrix operator- (Matrix m);
	Matrix operator* (Matrix &m);
	Matrix operator* (double x);
	void   operator*=(Matrix &m);
};

class SplitMatrix {
  public:
	SplitMatrix(Matrix &top_, Matrix &bottom_);
	SplitMatrix();
	Matrix top;
	Matrix bottom;
};

