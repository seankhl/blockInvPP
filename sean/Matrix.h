
#include <vector>
#include <cstdio>
using namespace std;

class Matrix 
{
  public:
	vector<vector<double>> rows;
	
    Matrix(size_t nRow, size_t nCol);
	
	void populateRow(size_t row, double val);
	void populateCol(size_t col, double val);
    void populateDiagonal(size_t row, size_t col, double val); 

	void appendCol(Matrix &col);
	//void pasteIn(Matrix &paste, size_t row, size_t col);
	void pasteIn(Matrix paste, size_t row, size_t col);
	
	vector<double> sliceRow(size_t row);
	Matrix sliceRows(size_t min, size_t max); //up to but not including b
    vector<double> sliceCol(size_t col);
	Matrix sliceBlock(size_t row,     size_t col, 
                      size_t numRows, size_t numCols);
	vector<Matrix> asRowBlocks(size_t numBlocks);

	Matrix operator+ (Matrix &m);
	Matrix operator- (Matrix &m);
	Matrix operator- (Matrix m);
	Matrix operator* (Matrix &m);
	Matrix operator* (double x);

    void print();
};

class SplitMatrix {
  public:
	SplitMatrix(Matrix &top_, Matrix &bottom_);
	Matrix top;
	Matrix bottom;
};

