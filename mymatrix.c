#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include "mylib.h"
#include "mymatrix.h"

// Displays elements of the matrix structure
void DisplayMatrix(Matrix matrix, int startRow, int startColumn, int endRow, int endColumn)
{
  int i; 
  for(i = startRow; i <= endRow; i++)
    {      
      puts("\n");
      int j;
      for(j = startColumn; j <= endColumn; j++)
	{
	  printf("%d\t",ValueFromMatrix(matrix, i, j));
	}
    }
}

// Display the submatrix
void DisplaySubMatrix(SubMatrix subMatrix, int startRow, int startColumn, int endRow, int endColumn)
{
  int i; 
  for(i = startRow; i <= endRow; i++)
    {      
      puts("\n");
      int j;
      for(j = startColumn; j <= endColumn; j++)
	{
	  printf("%d\t",ValueFromSubMatrix(subMatrix, i, j));
	}
    }
}

// Find maximum sequence in a 1-dimensional array
int FindMaxSubArray(SubArrayIndices *inputArray, int *array, int arrayLength, int start, int end)
{
  /*
  printf("\n Entered FindMaxSubArray");
  printf("\n Start: %d End: %d", start, end);
  */

  /********** Asserts **********/
  assert(array != NULL);
  assert(arrayLength > 0);
  assert(start >= 0 && start <= end);
  assert(end >= start && end < arrayLength);

  /********** Handle the base case **********/
  if(start == end)
    {
      inputArray->start = start;
      inputArray->end = start;
      return array[start];
    }


  int middleIndex = (start + end) / 2;

  SubArrayIndices maxLeftSubArrayIndices;
  int maxLeftSum = FindMaxSubArray(&maxLeftSubArrayIndices, array, arrayLength, start, middleIndex);

  SubArrayIndices maxRightSubArrayIndices;
  int maxRightSum = FindMaxSubArray(&maxRightSubArrayIndices, array, arrayLength, middleIndex + 1, end);

  SubArrayIndices maxMidSubArrayIndices;
  int maxMiddleSum = FindMaxMiddleSubArray(&maxMidSubArrayIndices, array, arrayLength, start, end);

  inputArray->start = maxLeftSubArrayIndices.start;
  inputArray->end = maxLeftSubArrayIndices.end;
  int result = maxLeftSum;

  if(maxRightSum > result)
    {
      inputArray->start = maxRightSubArrayIndices.start;
      inputArray->end = maxRightSubArrayIndices.end;
      result = maxRightSum;
    }

  if(maxMiddleSum > result)
    {
      inputArray->start = maxMidSubArrayIndices.start;
      inputArray->end = maxMidSubArrayIndices.end;
      result = maxMiddleSum;
    }

  return result;
}

//Graft a smaller matrix onto a larger one. inputs: starting row, starting column to start grafting
void GraftMatrix(Matrix matrix, int row, int column, Matrix smallerMatrix)
{
  assert(row + smallerMatrix.Rows - 1 <= matrix.Rows);
  assert(column + smallerMatrix.Columns - 1 <= matrix.Columns);

  int i; 
  int j; 

  for(i = 0; i < smallerMatrix.Rows; i++)
    {
      for(j = 0; j < smallerMatrix.Columns; j++)
	{
	  ValueToMatrix(matrix, row + i, column + j, ValueFromMatrix(smallerMatrix, i + 1, j + 1));
	}
    }
}

// Add 2 matrices
void MatrixAdd(Matrix *out, Matrix A, Matrix B)
{
  assert(A.Rows == B.Rows);
  assert(A.Columns == B.Columns);

  out->Rows = A.Rows;
  out->Columns = A.Columns;

  int i, j;

  for(i = 1; i <= A.Rows; i++)
    {
      for(j = 1; j <= A.Columns; j++)
	{
	  int A_ij = ValueFromMatrix(A, i, j);
	  int B_ij = ValueFromMatrix(B, i, j);
	  int sum = A_ij + B_ij;
	  ValueToMatrix(*out, i, j, sum);
	}
    }
}

// Multiplies 2 square matrices
void MatrixMultiply(Matrix *out, Matrix A, Matrix B)
{
  CREATESUBMATRIX(subA, A, 1, 1, A.Rows, A.Columns);
  CREATESUBMATRIX(subB, B, 1, 1, B.Rows, B.Columns);

  if(out->Rows > 256)
    {
      StrassenMatrixMultiply(out, A, B);
    }else{
    RecursiveMatrixMultiply(out, subA, subB);
  }
}

// Multiply the matrix by -1
Matrix *MatrixNegate(Matrix *out, Matrix A)
{
  int matrixSize = out->Rows * out->Columns;
  int i;
  for (i = 0; i < matrixSize; i++)
    {
      *(out->Data + i) *= -1; 
    }
  return out;
}

// Subtract 2 matrices
void MatrixSubtract(Matrix *out, Matrix A, Matrix B)
{
  assert(A.Rows == B.Rows);
  assert(A.Columns == B.Columns);

  out->Rows = A.Rows;
  out->Columns = A.Columns;

  int i, j;

  for(i = 1; i <= A.Rows; i++)
    {
      for(j = 1; j <= A.Columns; j++)
	{
	  int A_ij = ValueFromMatrix(A, i, j);
	  int B_ij = ValueFromMatrix(B, i, j);
	  int dif = A_ij - B_ij;
	  ValueToMatrix(*out, i, j, dif);
	}
    }
}

// Add more than 2 matrices. numArgs is the number of matrices being added. 
void MultipleMatrixAdd(Matrix *out,int numArgs, ...)
{
  va_list argList;
  va_start(argList, numArgs);
  
  memset(out->Data, 0, out->Rows * out->Columns * sizeof(int));
  
  int i;
  for(i = 0; i < numArgs; i++)
    {
      MatrixAdd(out, *out, va_arg(argList, Matrix));
    }
}

// Pad a matrix with 0's up to a square matrix 
void PadMatrix(Matrix *out, Matrix in)
{
  int i, j;
  for(i = 1; i <= out->Rows; i++)
    {
      for(j = 1; j <= out->Columns; j++)
	{
	  if(i <= in.Rows && j <= in.Columns)
	    {
	      ValueToMatrix(*out, i, j, ValueFromMatrix(in, i, j));
	    }
	  else
	    {
	      ValueToMatrix(*out, i, j, 0);
	    }
	}
    }
}

// Multiply two square submatrices
void RecursiveMatrixMultiply(Matrix *out, SubMatrix A, SubMatrix B)
{
  int aRows = A.EndRow - A.StartRow + 1;
  int aColumns = A.EndColumn - A.StartColumn + 1;
  int bRows = B.EndRow - B.StartRow + 1;
  int bColumns = B.EndColumn - B.StartColumn + 1;

  // Make sure A, B, and out are all square matrices of the same dimension. 
  assert(aRows == aColumns);
  assert(aColumns == bRows);
  assert(bRows  == bColumns);
  assert(bColumns  == out->Rows);
  assert(out->Rows  == out->Columns);

  if(!IsPowerOfTwo(aRows))
    {
      int newDimension = 1;
      while(newDimension < aRows)
	{
	  newDimension *=2;
	}
      
      CREATEMATRIX(A_M, aRows, aColumns);
      CREATEMATRIX(B_M, bRows, bColumns);
      CREATEMATRIX(A_Padded, newDimension, newDimension);
      CREATEMATRIX(B_Padded, newDimension, newDimension);
      CREATEMATRIX(out_Padded, newDimension, newDimension);

      SubMatrixToMatrix(&A_M, A);
      SubMatrixToMatrix(&B_M, B);
      
      PadMatrix(&A_Padded, A_M);
      PadMatrix(&B_Padded, B_M);

      MatrixMultiply(&out_Padded, A_Padded, B_Padded);
      UnpadMatrix(out,out_Padded);
     
      return;
    }
  
  if(1 == aRows)
    {
      int valueA = ValueFromSubMatrix(A, 1, 1);
      int valueB = ValueFromSubMatrix(B, 1, 1);
      *(out->Data) = valueA * valueB;
      out->Rows = 1;
      out->Columns = 1;
      return;
    }

  SubMatrix ASubMatrices[4];
   
  SubMatrix *ptr = ASubMatrices;
  SplitMatrix(A, ptr, ptr + 1, ptr + 2, ptr + 3);
   
  SubMatrix BSubMatrices[4];
   
  ptr = BSubMatrices;
  SplitMatrix(B, ptr, ptr + 1, ptr + 2, ptr + 3);

  int subMatrixDimensions = aRows / 2; 

  int outputData1[subMatrixDimensions][subMatrixDimensions], outputData2[subMatrixDimensions][subMatrixDimensions], resultData[subMatrixDimensions][subMatrixDimensions];
  Matrix outputMatrix1 = {.Data = &outputData1[0][0], .Rows = subMatrixDimensions, .Columns = subMatrixDimensions};
  Matrix outputMatrix2 = {.Data = &outputData2[0][0], .Rows = subMatrixDimensions, .Columns = subMatrixDimensions};
  Matrix resultMatrix = {.Data = &resultData[0][0], .Rows = subMatrixDimensions, .Columns = subMatrixDimensions};

  RecursiveMatrixMultiply(&outputMatrix1, ASubMatrices[0], BSubMatrices[0]);
  RecursiveMatrixMultiply(&outputMatrix2, ASubMatrices[1], BSubMatrices[2]);
  MatrixAdd(&resultMatrix, outputMatrix1, outputMatrix2);
  GraftMatrix(*out, 1, 1, resultMatrix);

  RecursiveMatrixMultiply(&outputMatrix1, ASubMatrices[0], BSubMatrices[1]);
  RecursiveMatrixMultiply(&outputMatrix2, ASubMatrices[1], BSubMatrices[3]);
  MatrixAdd(&resultMatrix, outputMatrix1, outputMatrix2);
  GraftMatrix(*out, 1, subMatrixDimensions + 1, resultMatrix);

  RecursiveMatrixMultiply(&outputMatrix1, ASubMatrices[2], BSubMatrices[0]);
  RecursiveMatrixMultiply(&outputMatrix2, ASubMatrices[3], BSubMatrices[2]);
  MatrixAdd(&resultMatrix, outputMatrix1, outputMatrix2);
  GraftMatrix(*out, subMatrixDimensions+1, 1, resultMatrix);

  RecursiveMatrixMultiply(&outputMatrix1, ASubMatrices[2], BSubMatrices[1]);
  RecursiveMatrixMultiply(&outputMatrix2, ASubMatrices[3], BSubMatrices[3]);
  MatrixAdd(&resultMatrix, outputMatrix1, outputMatrix2);
  GraftMatrix(*out, subMatrixDimensions + 1, subMatrixDimensions + 1, resultMatrix);
}

// Split a matrix into 4 submatrices.
void SplitMatrix(SubMatrix matrix, SubMatrix *matrix11, SubMatrix *matrix12, SubMatrix *matrix21, SubMatrix *matrix22)
{
  assert(matrix.EndRow - matrix.StartRow == matrix.EndColumn - matrix.StartColumn);
  int middleRow = (matrix.StartRow + matrix.EndRow) / 2;
  int middleColumn = (matrix.StartColumn + matrix.EndColumn) / 2;

  matrix11->Matrix = matrix.Matrix;
  matrix11->StartRow = matrix.StartRow;
  matrix11->StartColumn = matrix.StartColumn;
  matrix11->EndRow = middleRow;
  matrix11->EndColumn = middleColumn;

  matrix12->Matrix = matrix.Matrix;
  matrix12->StartRow = matrix.StartRow;
  matrix12->StartColumn = middleColumn + 1;
  matrix12->EndRow = middleRow;
  matrix12->EndColumn = matrix.EndColumn;

  matrix21->Matrix = matrix.Matrix;
  matrix21->StartRow = middleRow +1;
  matrix21->StartColumn = matrix.StartColumn;
  matrix21->EndRow = matrix.EndRow;
  matrix21->EndColumn = middleColumn;

  matrix22->Matrix = matrix.Matrix;
  matrix22->StartRow = middleRow + 1;
  matrix22->StartColumn = middleColumn + 1;
  matrix22->EndRow = matrix.EndRow;
  matrix22->EndColumn = matrix.EndColumn;
}

void StrassenMatrixMultiply(Matrix *out, Matrix A, Matrix B)
{
  CREATESUBMATRIX(A_Sub, A, 1, 1, A.Rows, A.Columns);
  CREATESUBMATRIX(B_Sub, B, 1, 1, B.Rows, B.Columns);

  StrassenSubMultiply(ouxot, A_Sub, B_Sub);
}
// Multiply two square submatrices using Strassen's algorithm (for larger matrices)
void StrassenSubMultiply(Matrix *out, SubMatrix A, SubMatrix B)
{
  int aRows = A.EndRow - A.StartRow + 1;
  int aColumns = A.EndColumn - A.StartColumn + 1;
  int bRows = B.EndRow - B.StartRow + 1;
  int bColumns = B.EndColumn - B.StartColumn + 1;

  // Make sure A, B, and out are all square matrices of the same dimension. 
  assert(aRows == aColumns);
  assert(aColumns == bRows);
  assert(bRows  == bColumns);
  assert(bColumns  == out->Rows);
  assert(out->Rows  == out->Columns);

  if(!IsPowerOfTwo(aRows))
    {
      int newDimension = 1;
      while(newDimension < aRows)
	{
	  newDimension *=2;
	}

      CREATEMATRIX(A_M, aRows, aColumns);
      CREATEMATRIX(B_M, bRows, bColumns);
      CREATEMATRIX(A_Padded, newDimension, newDimension);
      CREATEMATRIX(B_Padded, newDimension, newDimension);
      CREATEMATRIX(out_Padded, newDimension, newDimension);

      SubMatrixToMatrix(&A_M, A);
      SubMatrixToMatrix(&B_M, B);

      PadMatrix(&A_Padded, A_M);
      PadMatrix(&B_Padded, B_M);

      MatrixMultiply(&out_Padded, A_Padded, B_Padded);
      UnpadMatrix(out, out_Padded);
      return;
    }
  
  // Base case: both input matrices are of dimension 1
  if(1 == aRows)
    {
      int valueA = ValueFromSubMatrix(A, 1, 1);
      int valueB = ValueFromSubMatrix(B, 1, 1);
      *(out->Data) = valueA * valueB;
      out->Rows = 1;
      out->Columns = 1;
      return;
    }

  // Split the 2 input matrices into 4 equal submatrices
  // 0 => A_11, 1 => A_12, 2 => A_21, 3 => A_22
  SubMatrix ASubMatrices[4];
  SubMatrix *ptr = &ASubMatrices[0];
  SplitMatrix(A, ptr, ptr + 1, ptr + 2, ptr + 3);
   
  SubMatrix BSubMatrices[4];
  ptr = &BSubMatrices[0];
  SplitMatrix(B, ptr, ptr + 1, ptr + 2, ptr + 3);

  int subMatrixDimensions = aRows / 2; 

  // The output matrices are used for intermediate matrix results. 
  CREATEMATRIX(outputMatrix1, subMatrixDimensions, subMatrixDimensions);
  CREATEMATRIX(outputMatrix2, subMatrixDimensions, subMatrixDimensions);

  // Initialize the S matrices
  Matrix sMatrices[10];
  CREATEMATRIX(sMatrix1, subMatrixDimensions, subMatrixDimensions); sMatrices[0] = sMatrix1;
  CREATEMATRIX(sMatrix2, subMatrixDimensions, subMatrixDimensions); sMatrices[1] = sMatrix2; 
  CREATEMATRIX(sMatrix3, subMatrixDimensions, subMatrixDimensions); sMatrices[2] = sMatrix3;
  CREATEMATRIX(sMatrix4, subMatrixDimensions, subMatrixDimensions); sMatrices[3] = sMatrix4;
  CREATEMATRIX(sMatrix5, subMatrixDimensions, subMatrixDimensions); sMatrices[4] = sMatrix5;
  CREATEMATRIX(sMatrix6, subMatrixDimensions, subMatrixDimensions); sMatrices[5] = sMatrix6;
  CREATEMATRIX(sMatrix7, subMatrixDimensions, subMatrixDimensions); sMatrices[6] = sMatrix7;
  CREATEMATRIX(sMatrix8, subMatrixDimensions, subMatrixDimensions); sMatrices[7] = sMatrix8
								      CREATEMATRIX(sMatrix9, subMatrixDimensions, subMatrixDimensions); sMatrices[8] = sMatrix9; 
  CREATEMATRIX(sMatrix10, subMatrixDimensions, subMatrixDimensions); sMatrices[9] = sMatrix10;
   
  // Initialize the P matrices
  CREATEMATRIX(pMatrix1, subMatrixDimensions, subMatrixDimensions);
  CREATEMATRIX(pMatrix2, subMatrixDimensions, subMatrixDimensions);
  CREATEMATRIX(pMatrix3, subMatrixDimensions, subMatrixDimensions);
  CREATEMATRIX(pMatrix4, subMatrixDimensions, subMatrixDimensions);
  CREATEMATRIX(pMatrix5, subMatrixDimensions, subMatrixDimensions);
  CREATEMATRIX(pMatrix6, subMatrixDimensions, subMatrixDimensions);
  CREATEMATRIX(pMatrix7, subMatrixDimensions, subMatrixDimensions);

  /*** End of matrix initialiations ***/ 

  /******************** Strassen algorithm computation ********************/
  // B_12 - B22 
  MatrixSubtract(&sMatrices[0], *SubMatrixToMatrix(&outputMatrix1, BSubMatrices[1]), *SubMatrixToMatrix(&outputMatrix2, BSubMatrices[3]));

  // A_11 + A+12
  MatrixAdd(&sMatrices[1], *SubMatrixToMatrix(&outputMatrix1, ASubMatrices[0]), *SubMatrixToMatrix(&outputMatrix2, ASubMatrices[1]));

  // A_21 + A_22
  MatrixAdd(&sMatrices[2], *SubMatrixToMatrix(&outputMatrix1, ASubMatrices[2]), *SubMatrixToMatrix(&outputMatrix2, ASubMatrices[3]));

  // B_21 - B_11
  MatrixSubtract(&sMatrices[3], *SubMatrixToMatrix(&outputMatrix1, BSubMatrices[2]), *SubMatrixToMatrix(&outputMatrix2, BSubMatrices[0]));

  // A_11 + A_22
  MatrixAdd(&sMatrices[4], *SubMatrixToMatrix(&outputMatrix1, ASubMatrices[0]), *SubMatrixToMatrix(&outputMatrix2, ASubMatrices[3]));

  // B_11 + B_22
  MatrixAdd(&sMatrices[5], *SubMatrixToMatrix(&outputMatrix1, BSubMatrices[0]), *SubMatrixToMatrix(&outputMatrix2, BSubMatrices[3]));

  // A_12 - A_22
  MatrixSubtract(&sMatrices[6], *SubMatrixToMatrix(&outputMatrix1, ASubMatrices[1]), *SubMatrixToMatrix(&outputMatrix2, ASubMatrices[3]));

  // B_21 + B_22
  MatrixAdd(&sMatrices[7], *SubMatrixToMatrix(&outputMatrix1, BSubMatrices[2]), *SubMatrixToMatrix(&outputMatrix2, BSubMatrices[3]));

  // A_11 - A_21
  MatrixSubtract(&sMatrices[8], *SubMatrixToMatrix(&outputMatrix1, ASubMatrices[0]), *SubMatrixToMatrix(&outputMatrix2, ASubMatrices[2]));

  // B_11 + B_12
  MatrixAdd(&sMatrices[9], *SubMatrixToMatrix(&outputMatrix1, BSubMatrices[0]), *SubMatrixToMatrix(&outputMatrix2, BSubMatrices[1]));

  // These will be used as temporary submatrices to compute the P matrices 

  CREATESUBMATRIX(tempSubMatrix1, pMatrix1, 1, 1, subMatrixDimensions, subMatrixDimensions);
  CREATESUBMATRIX(tempSubMatrix2, pMatrix1, 1, 1, subMatrixDimensions, subMatrixDimensions); 

  // P_1 = A_11 * S_1
  tempSubMatrix1.Matrix = &sMatrices[0];
  StrassenSubMultiply(&pMatrix1, ASubMatrices[0], tempSubMatrix1);

  // P_2 = S_2 * B_22
  tempSubMatrix1.Matrix = &sMatrices[1];
  StrassenSubMultiply(&pMatrix2, tempSubMatrix1,BSubMatrices[3]);

  // P_3 = S_3 * B_11
  tempSubMatrix1.Matrix = &sMatrices[2];
  StrassenSubMultiply(&pMatrix3, tempSubMatrix1, BSubMatrices[0]);

  // P_4 = A_22 * S_4
  tempSubMatrix1.Matrix = &sMatrices[3];
  StrassenSubMultiply(&pMatrix4, ASubMatrices[3], tempSubMatrix1);

  // P_5 = S_5 * S_6
  tempSubMatrix1.Matrix = &sMatrices[4];
  tempSubMatrix2.Matrix = &sMatrices[5];
  StrassenSubMultiply(&pMatrix5, tempSubMatrix1, tempSubMatrix2);

  // P_6 = S_7 * S_8
  tempSubMatrix1.Matrix = &sMatrices[6];
  tempSubMatrix2.Matrix = &sMatrices[7];
  StrassenSubMultiply(&pMatrix6, tempSubMatrix1, tempSubMatrix2);

  // P_7 = S_9 * S_10
  tempSubMatrix1.Matrix = &sMatrices[8];
  tempSubMatrix2.Matrix = &sMatrices[9];
  StrassenSubMultiply(&pMatrix7, tempSubMatrix1, tempSubMatrix2);

  // C_11 = P_5 + P_4 - P_2 + P_6

  MultipleMatrixAdd(&outputMatrix1, 3, pMatrix5, pMatrix4, pMatrix6);
  MatrixSubtract(&outputMatrix1, outputMatrix1, pMatrix2);
  GraftMatrix(*out, 1, 1, outputMatrix1);

  // C_12 = P_1 + P_2
  MatrixAdd(&outputMatrix1, pMatrix1, pMatrix2);
  GraftMatrix(*out, 1, subMatrixDimensions + 1, outputMatrix1);

  // C_21 = P_3 + P_4
  MatrixAdd(&outputMatrix1, pMatrix3, pMatrix4);
  GraftMatrix(*out, subMatrixDimensions + 1, 1, outputMatrix1);

  // C_22 = P_5 + P_1 - P_3 - P_7
  MatrixAdd(&outputMatrix1, pMatrix5, pMatrix1);
  MatrixSubtract(&outputMatrix1, outputMatrix1, pMatrix3);
  MatrixSubtract(&outputMatrix1, outputMatrix1, pMatrix7);
  GraftMatrix(*out, subMatrixDimensions + 1, subMatrixDimensions + 1, outputMatrix1);   
}
 
// Convert a SubMatrix into a Matrix (with all its memory and performance cost) 
// out matrix must have data array preallocated and Row and Column propert
// output of this function is the same as out
Matrix *SubMatrixToMatrix(Matrix *out, SubMatrix subMatrix)
{
  int subMatrixRows = subMatrix.EndRow - subMatrix.StartRow + 1;
  int subMatrixColumns = subMatrix.EndColumn - subMatrix.StartColumn + 1;

  assert(out->Rows == subMatrixRows);
  assert(out->Columns == subMatrixColumns);

  int row;
  int column;

  for( row = 1; row <= out->Rows; row++)
    {
      for(column = 1; column <= out->Columns; column++)
	{
	  ValueToMatrix(*out, row, column, ValueFromSubMatrix(subMatrix, row, column));
	}
    }

  return out;
}

// Unpad 0's from matrix
void UnpadMatrix(Matrix *out, Matrix in)
{
  assert(out->Rows <= in.Rows);
  assert(out->Columns <= in.Columns);

  int i, j;
  for(i = 1; i <= out->Rows; i++)
    {
      for(j = 1; j <= out->Columns; j++)
	{
	  ValueToMatrix(*out, i, j, ValueFromMatrix(in, i, j));
	}
    }
}

// Returns element at i, j in the matrix structure (1-based index)
int ValueFromMatrix(Matrix matrix, int row, int column)
{
  assert(row <= matrix.Rows);
  assert(column <= matrix.Columns);
  return *(matrix.Data + (row - 1) * matrix.Columns + column - 1);
}

// Get a value from a submatrix
int ValueFromSubMatrix(SubMatrix subMatrix, int row, int column)
{
  int realRow = subMatrix.StartRow + row - 1;
  int realColumn = subMatrix.StartColumn + column - 1;
  return ValueFromMatrix(*(subMatrix.Matrix), realRow, realColumn);
}

// Send value to a matrix through a submatrix
void ValueToSubMatrix(SubMatrix subMatrix, int row, int column, int value)
{
  int realRow = subMatrix.StartRow + row - 1;
  int realColumn = subMatrix.StartColumn + column - 1;
  ValueToMatrix(*subMatrix.Matrix, realRow, realColumn, value);
}

// Send value to a matrix
void ValueToMatrix(Matrix matrix, int row, int column, int value)
{
  assert(row <= matrix.Rows);
  assert(column <= matrix.Columns);
  *(matrix.Data + (row - 1) * matrix.Columns + column - 1) = value;
} 

/******************** non-exposed functions ********************/

// Used by FindMaxSubArray
int FindMaxMiddleSubArray(SubArrayIndices *inputArray, int *array, int arrayLength, int start, int end)
{
  /********** Asserts **********/
  assert(array != NULL);
  assert(arrayLength > 0);
  assert(start >= 0 && start <= end);
  assert(end >= start && end < arrayLength);

  int middleIndex = (start + end) / 2; 

  int maxLeftSum = INT_MIN;
  int currentSum = 0;

  int index; 
  for(index = middleIndex; index >= start; index--)
    {
      currentSum += array[index];
      if(currentSum > maxLeftSum)
	{
	  maxLeftSum = currentSum;
	  inputArray->start = index;
	}
    }

  int maxRightSum = INT_MIN;
  currentSum = 0;

  for(index = middleIndex + 1; index <= end; index++)
    {
      currentSum += array[index];
      if(currentSum > maxRightSum)
	{
	  maxRightSum = currentSum;
	  inputArray->end = index;
	}
    }
  return maxLeftSum + maxRightSum;
}
