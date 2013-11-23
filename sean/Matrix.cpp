#include "Matrix.h"

/* constructor */
Matrix::Matrix(size_t nRow, size_t nCol) 
{
	for (size_t y=0; y<nRow; y++) {
		vector<double> row;
		row.reserve(nCol);
		for (size_t x=0; x<nCol; x++) {
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

void Matrix::populateRow(size_t nRow, double val) 
{
	for (size_t x=0; x<rows[0].size(); x++) {
		rows[nRow][x] = val;
	}
}

void Matrix::populateCol(size_t nCol, double val) 
{
	for (size_t y=0; y<rows.size(); y++) {
		rows[y][nCol] = val;
	}
}

void Matrix::populateDiagonal(size_t row, size_t col, double val) 
{
	size_t y = row;
	size_t x = col;
	while (y < rows.size() && x < rows[0].size()) {
		rows[y][x] = val;
		y++;
		x++;
	}
}

/************************************************************************/

void Matrix::appendCol(Matrix &col) 
{
	for (size_t y=0; y<rows.size(); y++) {
		rows[y].push_back(col.rows[y][0]);
	}
}

void Matrix::pasteIn(Matrix paste, size_t row, size_t col)
{
	for (size_t y=0; y<paste.rows.size(); y++) {
		for (size_t x=0; x<paste.rows[0].size(); x++) {
			rows[row+y][col+x] = paste.rows[y][x];
		}
	}
}

/************************************************************************/

// get row nRow
vector<double> Matrix::sliceRow(size_t nRow) 
{
	vector<double> destRow = vector<double>(rows[0].size());
	destRow = rows[nRow];
	return destRow;
}

// get rows min through max (in/ex-clusive, respectively)
Matrix Matrix::sliceRows(size_t min, size_t max) 
{
	Matrix sliced = Matrix(max - min, rows[0].size());
	for (size_t i=min; i<max; i++) {
		sliced.rows[i - min] = rows[i];
	}
	return sliced;
}

// get col nCol
vector<double> Matrix::sliceCol(size_t nCol) 
{
	vector<double> destCol = vector<double>(rows.size(), 1);
	for (size_t y=0; y<rows.size(); y++) {
		destCol[y] = rows[y][nCol];
	}
	return destCol;
}

// get submatrix from row to nRows and col to nCols
Matrix Matrix::sliceBlock(size_t row,   size_t col, 
                          size_t nRows, size_t nCols) 
{
	Matrix block = Matrix(nRows, nCols);
	for (size_t y=0; y<nRows; y++) {
		for (size_t x=0; x<nCols; x++) {
			block.rows[y][x] = rows[row + y][col + x];
		}
	}
	return block;
}

vector<Matrix> Matrix::asRowBlocks(size_t blockSize) 
{
	vector<Matrix> blocks;
	for (size_t i=0; i<rows.size()/blockSize; i++) {
		Matrix cpRows = Matrix(blockSize, rows[0].size());
		for (size_t j=0; j<blockSize; j++) {
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
	for (size_t y=0; y<rows.size(); y++) {
		for (size_t x=0; x<rows[0].size();  x++) {
			res.rows[y][x] = rows[y][x] + m.rows[y][x];
		}
	}
	return res;
}

Matrix Matrix::operator- (Matrix &m) 
{
	Matrix res = Matrix(rows.size(), rows[0].size());
	for (size_t y=0; y<rows.size(); y++) {
		for (size_t x=0; x<rows[0].size();  x++) {
			res.rows[y][x] = rows[y][x] - m.rows[y][x];
		}
	}
	return res;
}
Matrix Matrix::operator- (Matrix m) 
{
	Matrix res = Matrix(rows.size(), rows[0].size());
	for (size_t y=0; y<rows.size(); y++) {
		for (size_t x=0; x<rows[0].size();  x++) {
			res.rows[y][x] = rows[y][x] - m.rows[y][x];
		}
	}
	return res;
}


Matrix Matrix::operator* (Matrix &m) 
{
	Matrix res = Matrix(rows.size(), m.rows[0].size());

	for (size_t y=0; y<res.rows.size(); y++) {
		for (size_t x=0; x<res.rows[0].size(); x++) {
			double sum = 0;
			for (size_t my=0; my<m.rows.size(); my++) {
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
	for (size_t row=0; row<rows.size(); row++) {
		for (size_t col=0; col<rows[0].size(); col++) {
			res.rows[row][col] = rows[row][col] * x;
		}
	}
	return res;
}

void Matrix::print()
{
	for (size_t row=0; row<rows.size(); row++) {
        if (row == 0) {
            printf("[ ");
        }
        else {
            printf("  ");
        }
		for (size_t col=0; col<rows[0].size(); col++) {
			printf("% f ",rows[row][col]);
		}
        if (row == rows.size()-1) {
            printf("]\n");
        }
        else {
            printf("\n");
        }
	}
}

