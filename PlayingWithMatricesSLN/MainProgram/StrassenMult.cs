
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;


namespace MainProgram
{
	public static class StrassenMult
	{

		private class SubMatrix
		{

			public	int			Length;
			public	int			StartX;
			public	int			StartY;
			public	int			RealRowWidth;

			public	float[]		RealMatrix;


			public static void Add(SubMatrix leftSide, SubMatrix rightSide, SubMatrix result)
			{
				var offsetLH	= leftSide.StartY * leftSide.RealRowWidth + leftSide.StartX;
				var offsetRH	= rightSide.StartY * rightSide.RealRowWidth + rightSide.StartX;
				var offsetR		= result.StartY * result.RealRowWidth + result.StartX;

				for (int i = 0; i < leftSide.Length; i++)
				{
					for (int j = 0; j < leftSide.Length; j++)
						result.RealMatrix[offsetR + j] = leftSide.RealMatrix[offsetLH + j] + rightSide.RealMatrix[offsetRH + j];

					offsetLH	+= leftSide.RealRowWidth;
					offsetRH	+= rightSide.RealRowWidth;
					offsetR		+= result.RealRowWidth;
				}
			}

			public static void Subtract(SubMatrix leftSide, SubMatrix rightSide, SubMatrix result)
			{
				var offsetLH    = leftSide.StartY * leftSide.RealRowWidth + leftSide.StartX;
				var offsetRH    = rightSide.StartY * rightSide.RealRowWidth + rightSide.StartX;
				var offsetR     = result.StartY * result.RealRowWidth + result.StartX;

				for (int i = 0; i < leftSide.Length; i++)
				{
					for (int j = 0; j < leftSide.Length; j++)
						result.RealMatrix[offsetR + j] = leftSide.RealMatrix[offsetLH + j] - rightSide.RealMatrix[offsetRH + j];

					offsetLH	+= leftSide.RealRowWidth;
					offsetRH	+= rightSide.RealRowWidth;
					offsetR		+= result.RealRowWidth;
				}
			}

			public static void MatrixMult(SubMatrix leftSide, SubMatrix rightSide, SubMatrix result)
			{
				throw new NotImplementedException();
			}

		}


		public static float[] MatrixMult(float[] matrixA, float[] matrixB, Tuple<int, int> shapeA, Tuple<int, int> shapeB)
		{
			if (shapeA.Item2 != shapeB.Item1)
				throw new Exception();

			var resultRowCount      = shapeA.Item1;
			var resultColumnCount   = shapeB.Item2;
			var sharedDimension     = shapeA.Item2;
			var resultBlockSize		= resultRowCount * resultColumnCount;

			var maxD = shapeA.Item1;
			maxD = maxD < shapeA.Item2 ? shapeA.Item2 : maxD;
			maxD = maxD < shapeB.Item2 ? shapeB.Item2 : maxD;

			int n = 0;
			int t = 1;

			float[] result = new float[resultBlockSize];

			while (t < maxD)
			{
				t *= 2;
				n++;
			}

			if (t < 129)
			{
				var b_t = new float[shapeB.Item1 * shapeB.Item2];

				MatrixMult_TransposeThenDotProducts(matrixA, matrixB, result, b_t, shapeA, shapeB);
				b_t = null;

				return result;
			}

			int childLength		= t / 2;
			int childBlockSize	= childLength * childLength;
			int childN			= n - 1;

			var a11		= new SubMatrix() { RealMatrix = MonoMatrixOperations.CreateZeroMatrix(childBlockSize), StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var a12		= new SubMatrix() { RealMatrix = MonoMatrixOperations.CreateZeroMatrix(childBlockSize), StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var a21		= new SubMatrix() { RealMatrix = MonoMatrixOperations.CreateZeroMatrix(childBlockSize), StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var a22		= new SubMatrix() { RealMatrix = MonoMatrixOperations.CreateZeroMatrix(childBlockSize), StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };

			var b11		= new SubMatrix() { RealMatrix = MonoMatrixOperations.CreateZeroMatrix(childBlockSize), StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var b12		= new SubMatrix() { RealMatrix = MonoMatrixOperations.CreateZeroMatrix(childBlockSize), StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var b21		= new SubMatrix() { RealMatrix = MonoMatrixOperations.CreateZeroMatrix(childBlockSize), StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var b22		= new SubMatrix() { RealMatrix = MonoMatrixOperations.CreateZeroMatrix(childBlockSize), StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };

			var c11		= a11;
			var c12		= a12;
			var c21		= a21;
			var c22		= a22;

			var tempL	= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var tempR	= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };

			var m1		= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var m2		= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var m3		= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var m4		= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var m5		= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var m6		= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };
			var m7		= new SubMatrix() { RealMatrix = new float[childBlockSize], StartX = 0, StartY = 0, Length = childLength, RealRowWidth = childLength };


			///
			/// Break matrices { A, B } into matrices { A11, A12, A21, A22, B11, B12, B21, B22 }
			///
			FillInBlock(result: a11.RealMatrix, source: matrixA, startRow: 0,			startColumn: 0,				length: childLength, sourceHeight: shapeA.Item1, sourceWidth: shapeA.Item2);
			FillInBlock(result: a12.RealMatrix, source: matrixA, startRow: 0,			startColumn: childLength,	length: childLength, sourceHeight: shapeA.Item1, sourceWidth: shapeA.Item2);
			FillInBlock(result: a21.RealMatrix, source: matrixA, startRow: childLength, startColumn: 0,				length: childLength, sourceHeight: shapeA.Item1, sourceWidth: shapeA.Item2);
			FillInBlock(result: a22.RealMatrix, source: matrixA, startRow: childLength, startColumn: childLength,	length: childLength, sourceHeight: shapeA.Item1, sourceWidth: shapeA.Item2);

			FillInBlock(result: b11.RealMatrix, source: matrixB, startRow: 0,			startColumn: 0,				length: childLength, sourceHeight: shapeB.Item1, sourceWidth: shapeB.Item2);
			FillInBlock(result: b12.RealMatrix, source: matrixB, startRow: 0,			startColumn: childLength,	length: childLength, sourceHeight: shapeB.Item1, sourceWidth: shapeB.Item2);
			FillInBlock(result: b21.RealMatrix, source: matrixB, startRow: childLength, startColumn: 0,				length: childLength, sourceHeight: shapeB.Item1, sourceWidth: shapeB.Item2);
			FillInBlock(result: b22.RealMatrix, source: matrixB, startRow: childLength, startColumn: childLength,	length: childLength, sourceHeight: shapeB.Item1, sourceWidth: shapeB.Item2);


			///
			/// Calculte m1 to m7
			/// 
			SubMatrix.Add(a11, a22, tempL);
			SubMatrix.Add(b11, b22, tempR);
			//MatrixMult_StrassenRecursiveComp(tm1, tm2, m1, childLength);

			SubMatrix.Add(a21, a22, tempL);
			//MatrixMult_StrassenRecursiveComp(tm1, b11, m2, childLength);

			SubMatrix.Subtract(b12, b22, tempR);
			//MatrixMult_StrassenRecursiveComp(a11, tm2, m3, childLength);

			SubMatrix.Subtract(b21, b11, tempR);
			//MatrixMult_StrassenRecursiveComp(a22, tm2, m4, childLength);

			SubMatrix.Add(a11, a12, tempL);
			//MatrixMult_StrassenRecursiveComp(tm1, b22, m5, childLength);

			SubMatrix.Subtract(a21, a11, tempL);
			SubMatrix.Add(b11, b12, tempR);
			//MatrixMult_StrassenRecursiveComp(tm1, tm2, m6, childLength);

			SubMatrix.Subtract(a12, a22, tempL);
			SubMatrix.Add(b21, b22, tempR);
			//MatrixMult_StrassenRecursiveComp(tm1, tm2, m7, childLength);


			///
			/// Calculate c11, c12, c21, c22
			/// 
			SubMatrix.Add(m1, m4, tempL);
			SubMatrix.Subtract(m7, m5, tempR);
			SubMatrix.Add(tempL, tempR, c11);

			SubMatrix.Add(m3, m5, c12);

			SubMatrix.Add(m2, m4, c21);

			SubMatrix.Subtract(m1, m2, tempL);
			SubMatrix.Add(m3, m6, tempR);
			SubMatrix.Add(tempL, tempR, c22);


			///
			/// Transfer { c11, c12, c21, c22 } to result
			/// 
			for (int i = 0; i < childLength; i++)
			{
				if (resultRowCount <= i)
					break;

				int offsetR = i * resultColumnCount;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					if (resultColumnCount <= j)
						break;

					result[offsetR + j] = c11.RealMatrix[offsetC + j];
				}
			}

			for (int i = 0; i < childLength; i++)
			{
				if (i >= resultRowCount)
					break;

				int offsetR = (i * resultColumnCount) + childLength;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					if (j + childLength >= resultColumnCount)
						break;

					result[offsetR + j] = c12.RealMatrix[offsetC + j];
				}
			}

			for (int i = 0; i < childLength; i++)
			{
				if (i + childLength >= resultRowCount)
					break;

				int offsetR = (i + childLength) * resultColumnCount;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					if (j >= resultColumnCount)
						break;

					result[offsetR + j] = c21.RealMatrix[offsetC + j];
				}
			}

			for (int i = 0; i < childLength; i++)
			{
				if (i + childLength >= resultRowCount)
					break;

				int offsetR = ((i + childLength) * resultColumnCount) + childLength;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					if (j + childLength >= resultColumnCount)
						break;

					result[offsetR + j] = c22.RealMatrix[offsetC + j];
				}
			}


			///
			/// Return the result
			/// 
			return result;
		}

		private static void MatrixMult(SubMatrix leftSide, SubMatrix rightSide, SubMatrix result)
		{
			throw new NotImplementedException();
		}

		private static void MatrixMult_TransposeThenDotProducts(float[] matrixA, float[] matrixB, float[] result, float[] b_T, Tuple<int, int> shapeA, Tuple<int, int> shapeB)
		{
			var resultRowCount      = shapeA.Item1;
			var resultColumnCount   = shapeB.Item2;
			var sharedDimension     = shapeA.Item2;

			var copyRowOffset		= 0;

			for (int i = 0; i < shapeB.Item1; i++)
			{
				for (int j = 0; j < shapeB.Item2; j++)
					b_T[j * shapeB.Item1 + i] = matrixB[copyRowOffset + j];
				copyRowOffset += shapeB.Item2;
			}

			var resultSize  = resultRowCount * resultColumnCount;

			int rI          = 0;
			int aRowOffset  = 0;

			var elementCount = Vector<float>.Count;
			var simdCount = sharedDimension / elementCount;

			for (int rY = 0; rY < resultRowCount; rY++)
			{
				int bRowOffset = 0;
				for (int rX = 0; rX < resultColumnCount; rX++)
				{
					result[rI] = 0.0f;

					int i_sd = 0;
					for (int i_e = 0; i_e < simdCount; i_e++)
					{
						var vectorA = new Vector<float>(matrixA, aRowOffset + i_sd);
						var vectorB = new Vector<float>(b_T, bRowOffset + i_sd);

						result[rI] += Vector.Dot(vectorA, vectorB);

						i_sd += elementCount;
					}
					for (; i_sd < sharedDimension; i_sd++)
						result[rI] += matrixA[aRowOffset + i_sd] * b_T[bRowOffset + i_sd];
					rI++;
					bRowOffset += sharedDimension;
				}
				aRowOffset += sharedDimension;
			}
		}


		private static void FillInBlock(float[] result, float[] source, int startColumn, int startRow, int length, int sourceWidth, int sourceHeight)
		{
			for (int i = 0; i < length; i++)
			{
				var sourceRow = i + startRow;

				if (sourceHeight <= sourceRow)
					break;

				for (int j = 0; j < length; j++)
				{
					var sourceColumn = j + startColumn;
					if (sourceWidth <= sourceColumn)
						break;

					result[i * length + j] = source[sourceRow * sourceWidth + sourceColumn];
				}
			}
		}


		private static void Add(float[] leftSide, float[] rightSide, float[] result, int blockSize)
		{
			int packedCount = Vector<float>.Count;

			for (int i = 0; i < blockSize; i += packedCount)
			{
				var l = new Vector<float>(leftSide, i);
				var r = new Vector<float>(rightSide, i);

				(l + r).CopyTo(result, i);
			}
		}

		private static void Add(float[] leftSide, float[] rightSide, float[] result, int length, int rowSize, Tuple<int, int> leftStart, Tuple<int, int> rightStart)
		{
			int packedCount = Vector<float>.Count;
			int resultIndex = 0;

			int offsetL = leftStart.Item1 * rowSize + leftStart.Item2;
			int offsetR = rightStart.Item1 * rowSize + rightStart.Item2;

			for (int i = 0; i < length; i++)
			{
				for (int j = 0; j < length; j += packedCount)
				{
					var l = new Vector<float>(leftSide, offsetL + j);
					var r = new Vector<float>(rightSide, offsetR + j);

					(l + r).CopyTo(result, resultIndex);

					resultIndex += packedCount;
				}

				offsetL += rowSize;
				offsetR += rowSize;
			}
		}

		private static void Subtract(float[] leftSide, float[] rightSide, float[] result, int blockSize)
		{
			int packedCount     = Vector<float>.Count;

			for (int i = 0; i < blockSize; i += packedCount)
			{
				var l = new Vector<float>(leftSide, i);
				var r = new Vector<float>(rightSide, i);

				(l - r).CopyTo(result, i);
			}
		}

		private static void Subtract(float[] leftSide, float[] rightSide, float[] result, int length, int rowSize, Tuple<int, int> leftStart, Tuple<int, int> rightStart)
		{
			int packedCount = Vector<float>.Count;
			int resultIndex = 0;

			int offsetL = leftStart.Item1 * rowSize + leftStart.Item2;
			int offsetR = rightStart.Item1 * rowSize + rightStart.Item2;

			for (int i = 0; i < length; i++)
			{
				for (int j = 0; j < length; j += packedCount)
				{
					var l = new Vector<float>(leftSide, offsetL + j);
					var r = new Vector<float>(rightSide, offsetR + j);

					(l - r).CopyTo(result, resultIndex);

					resultIndex += packedCount;
				}

				offsetL += rowSize;
				offsetR += rowSize;
			}
		}

	}
}
