#ifndef _MYMATRIX_H_
#define _MYMATRIX_H_
 
typedef struct
{
  int start;
  int end;
} SubArrayIndices;

typedef struct
{
  int *Data;
  int Rows;
  int Columns;
}Matrix;

typedef struct
{
  Matrix *Matrix;
  int StartRow;
  int StartColumn;
  int EndRow;
  int EndColumn;
}SubMatrix;

// FuncDecs
void DisplayMatrix(Matrix matrix, int startRow, int startColumn, int endRow,int endColumn);
int FindMaxSubArray(SubArrayIndices *inputArray, int *array, int arrayLength, int start, int end);
int FindMaxMiddleSubArray(SubArrayIndices *inputArray, int *array, int arrayLength, int start, int end);
void GraftMatrix(Matrix matrix, int row, int column, Matrix smallerMatrix);
void MatrixAdd(Matrix *out, Matrix A, Matrix B);
void MatrixSubtract(Matrix *out, Matrix A, Matrix B);
void SplitMatrix(SubMatrix matrix, SubMatrix *matrix11, SubMatrix *matrix12,SubMatrix *matrix21, SubMatrix *matrix22);
void MatrixMultiply(Matrix *out, Matrix A, Matrix B);
void MultipleMatrixAdd(Matrix *out,int numArgs, ...);
Matrix *MatrixNegate(Matrix *out, Matrix A);
void PadMatrix(Matrix *out, Matrix in);
void RecursiveMatrixMultiply(Matrix *out, SubMatrix A, SubMatrix B);
void StrassenMatrixMultiply(Matrix *out, Matrix A, Matrix B);
void StrassenSubMultiply(Matrix *out, SubMatrix A, SubMatrix B);
Matrix *SubMatrixToMatrix(Matrix *out, SubMatrix A);
void UnpadMatrix(Matrix *out, Matrix in);
int ValueFromMatrix(Matrix matrix, int row, int column);
int ValueFromSubMatrix(SubMatrix subMatrix, int row, int column);
void ValueToMatrix(Matrix matrix, int row, int column, int value);
void ValueToSubMatrix(SubMatrix subMatrix, int row, int column, int value);

#define CREATEMATRIX(matrixName, rows, columns); Matrix matrixName; int matrixName##Data[(rows)][(columns)]; matrixName.Data = &matrixName##Data[0][0]; matrixName.Rows = (rows); matrixName.Columns = (columns)

#define CREATESUBMATRIX(subMatrixName, dataName, startRow, startColumn, endRow, endColumn); SubMatrix subMatrixName = {.Matrix = &dataName, .StartRow = (startRow), .StartColumn = (startColumn), .EndRow = (endRow), .EndColumn = (endColumn)}
#endif

