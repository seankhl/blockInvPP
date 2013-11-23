#include "Matrix.h"

/* constructor */
Matrix::Matrix(int nRow, int nCol) 
{
	for (int y=0; y<nRow; y++) {
		vector<double> row;
		row.reserve(nCol);
		for (int x=0; x<nCol; x++) {
			row.push_back(0);
		}
		rows.push_back(row);
	}
}

SplitMatrix::SplitMatrix(Matrix &top_, Matrix &bottom_) 
  : top(top_), bottom(bottom_) 
{
	this->top = top;
	this->bottom = bottom;
}

/************************************************************************/

void Matrix::populateRow(int nRow, double val) 
{
	for (int x=0; x<rows[0].size(); x++) {
		rows[nRow][x] = val;
	}
}

void Matrix::populateCol(int nCol, double val) 
{
	for (int y=0; y<rows.size(); y++) {
		rows[y][nCol] = val;
	}
}

void Matrix::populateDiagonal(int row, int col, double val) 
{
	int y = row;
	int x = col;
	while (y < rows.size() && x < rows[0].size()) {
		rows[y][x] = val;
		y++;
		x++;
	}
}

/************************************************************************/

void Matrix::appendCol(Matrix &col) 
{
	for (int y=0; y<rows.size(); y++) {
		rows[y].push_back(col.rows[y][0]);
	}
}

void Matrix::pasteIn(Matrix paste, int row, int col)
{
	for (int y=0; y<paste.rows.size(); y++) {
		for (int x=0; x<paste.rows[0].size(); x++) {
			rows[row+y][col+x] = paste.rows[y][x];
		}
	}
}

/************************************************************************/

Matrix Matrix::sliceRow(int nRow) 
{
	Matrix destRow = Matrix(1, rows[0].size());
	destRow.rows[0] = rows[nRow];
	return destRow;
}

Matrix Matrix::sliceCol(int nCol) 
{
	Matrix destCol = Matrix(rows.size(), 1);
	for (int y=0; y<rows.size(); y++) {
		destCol.rows[y][0] = rows[y][nCol];
	}
	return destCol;
}

Matrix Matrix::sliceRows(int min, int max) 
{
	Matrix sliced = Matrix(max - min, rows[0].size());
	for (int i=min; i<max; i++) {
		sliced.rows[i - min] = rows[i];
	}
	return sliced;
}

Matrix Matrix::sliceBlock(int row, int col, int nRows, int nCols) 
{
	Matrix block = Matrix(nRows, nCols);
	for (int y=0; y<nRows; y++) {
		for (int x=0; x<nCols; x++) {
			block.rows[y][x] = rows[row + y][col + x];
		}
	}
	return block;
}

vector<Matrix> Matrix::asRowBlocks(int blockSize) 
{
	vector<Matrix> blocks;
	for (int i=0; i<rows.size()/blockSize; i++) {
		Matrix cpRows = Matrix(blockSize, rows[0].size());
		for (int j=0; j<blockSize; j++) {
			cpRows.rows[j] = rows[i * blockSize + j];
		}
		blocks.push_back(cpRows);
	}
	return blocks;
}

/************************************************************************/

Matrix Matrix::operator+ (Matrix &m) 
{
	Matrix res = Matrix(rows.size(), rows[0].size());
	for (int y=0; y<rows.size(); y++) {
		for (int x=0; x<rows[0].size();  x++) {
			res.rows[y][x] = rows[y][x] + m.rows[y][x];
		}
	}
	return res;
}

Matrix Matrix::operator- (Matrix &m) 
{
	Matrix res = Matrix(rows.size(), rows[0].size());
	for (int y=0; y<rows.size(); y++) {
		for (int x=0; x<rows[0].size();  x++) {
			res.rows[y][x] = rows[y][x] - m.rows[y][x];
		}
	}
	return res;
}
Matrix Matrix::operator- (Matrix m) 
{
	Matrix res = Matrix(rows.size(), rows[0].size());
	for (int y=0; y<rows.size(); y++) {
		for (int x=0; x<rows[0].size();  x++) {
			res.rows[y][x] = rows[y][x] - m.rows[y][x];
		}
	}
	return res;
}


Matrix Matrix::operator* (Matrix &m) 
{
	Matrix res = Matrix(rows.size(), m.rows[0].size());

	for (int y=0; y<res.rows.size(); y++) {
		for (int x=0; x<res.rows[0].size(); x++) {
			double sum = 0;
			for (int my=0; my<m.rows.size(); my++) {
				sum += m.rows[my][x] * rows[y][my];
			}
			res.rows[y][x] = sum;
		}
	}
	return res;
}

Matrix Matrix::operator* (double x) 
{
	Matrix res = Matrix(rows.size(), rows[0].size());
	for (int row=0; row<rows.size(); row++) {
		for (int col=0; col<rows[0].size(); col++) {
			res.rows[row][col] = rows[row][col] * x;
		}
	}
	return res;
}


