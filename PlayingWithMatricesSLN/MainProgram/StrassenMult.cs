
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

		/// <summary>
		/// This is a wrapper to part (or all) of an 'N' by 'N' matrix. It is assumed that 'N' is a power of 2 
		/// and that Vector&lt;float&gt; fits inside of 'N' nicely. We are assuming the same for Length.
		/// </summary>
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

			public static void MatrixMult(SubMatrix leftSide, SubMatrix rightSide, SubMatrix rightSide_T, SubMatrix result)
			{
				/// 
				/// Note:
				/// - leftSide.Length == rightSide.Length == rightSide_T.Length
				/// - rightSide_T.RealMatrix.Length == (leftSide.Length * leftSide.Length)
				/// 
				/// SubMatrix rightSide_T was created to be the size at which we switch to more conventional matrix 
				/// multiplication. More conventional matrix multiplication runs faster if we transpose the right hand
				/// side matrix and use that instead because it cuts down on the number of cache misses and 
				/// cache-reloads.
				/// 

				///
				/// Transpose the rightSide SubMatrix onto rightSide_T.
				/// 
				for (int i = 0; i < leftSide.Length; i++)
				{
					int sourceRow		= rightSide.StartY + i;
					int sourceOffset	= sourceRow * rightSide.RealRowWidth + rightSide.StartX;

					int offsetT = 0;
					for (int j = 0; j < leftSide.Length; j++)
					{
						//rightSide_T.RealMatrix[j * leftSide.Length + i] = rightSide.RealMatrix[sourceOffset + j];
						rightSide_T.RealMatrix[offsetT + i] = rightSide.RealMatrix[sourceOffset + j];

					}
					offsetT += leftSide.Length;
				}

				var elementCount	= Vector<float>.Count;
				var simdCount		= leftSide.Length / elementCount;

				var resultIndexOffset	= result.StartY * result.RealRowWidth + result.StartX;
				var leftIndexOffset		= leftSide.StartY * leftSide.RealRowWidth + leftSide.StartX;
				for (int rY = 0; rY < leftSide.Length; rY++)
				{
					//var resultIndexOffset = result.StartY + rY;
					//resultIndexOffset *= result.RealRowWidth;
					//resultIndexOffset += result.StartX;

					//var leftIndexOffset = leftSide.StartY + rY;
					//leftIndexOffset *= leftSide.RealRowWidth;
					//leftIndexOffset += leftSide.StartX;

					int bRowOffset = 0;

					for (int rX = 0; rX < leftSide.Length; rX++)
					{
						result.RealMatrix[resultIndexOffset + rX] = 0.0f;

						for (int rS = 0; rS < leftSide.Length; rS += elementCount)
						{
							var vectorA = new Vector<float>(leftSide.RealMatrix, leftIndexOffset + rS);
							var vectorB = new Vector<float>(rightSide_T.RealMatrix, bRowOffset + rS);

							result.RealMatrix[resultIndexOffset + rX] += Vector.Dot(vectorA, vectorB);
						}

						bRowOffset += leftSide.Length;
					}

					resultIndexOffset += result.RealRowWidth;
					leftIndexOffset += leftSide.RealRowWidth;
				}
			}

		}


		/// <summary>
		/// This needs to be equal to "2^n + 1" and "Vector&lt;float&gt;.Count &lt;= 2^n".
		/// </summary>
		private static int SwitchLength = 65;

		public static float[] MatrixMult(float[] matrixA, float[] matrixB, Tuple<int, int> shapeA, Tuple<int, int> shapeB)
		{
			if (shapeA.Item2 != shapeB.Item1)
				throw new Exception();

			var resultRowCount      = shapeA.Item1;
			var resultColumnCount   = shapeB.Item2;
			var sharedDimension     = shapeA.Item2;
			var resultBlockSize		= resultRowCount * resultColumnCount;

			// Find the largest dimension to use as our starting point for length 't'
			var maxD = shapeA.Item1;
			maxD = maxD < shapeA.Item2 ? shapeA.Item2 : maxD;
			maxD = maxD < shapeB.Item2 ? shapeB.Item2 : maxD;

			// Find length 't'
			int n = 0;
			int t = 1;

			float[] result = new float[resultBlockSize];

			while (t < maxD)
			{
				t *= 2;
				n++;
			}

			// If 't' is less than our switch length, then use a more "conventional" way to multiply matrices
			if (t < SwitchLength)
			{
				var b_t = new float[shapeB.Item1 * shapeB.Item2];

				MatrixMult_TransposeThenDotProducts(matrixA, matrixB, result, b_t, shapeA, shapeB);
				b_t = null;

				return result;
			}

			// Create the initial matrices and sub-matrices
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

			// By the time we are using 'C's sub-matrices, we no longer need 'A's sub-matrices, so we might as well re-use them.
			var c11		= a11;
			var c12		= a12;
			var c21		= a21;
			var c22		= a22;

			int limitLength = SwitchLength - 1;
			var b_T		= new SubMatrix() { RealMatrix = new float[limitLength * limitLength], StartX = 0, StartY = 0, Length = limitLength, RealRowWidth = limitLength };

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
			MatrixMult(tempL, tempR, b_T, m1);

			SubMatrix.Add(a21, a22, tempL);
			MatrixMult(tempL, b11, b_T, m2);

			SubMatrix.Subtract(b12, b22, tempR);
			MatrixMult(a11, tempR, b_T, m3);

			SubMatrix.Subtract(b21, b11, tempR);
			MatrixMult(a22, tempR, b_T, m4);

			SubMatrix.Add(a11, a12, tempL);
			MatrixMult(tempL, b22, b_T, m5);

			SubMatrix.Subtract(a21, a11, tempL);
			SubMatrix.Add(b11, b12, tempR);
			MatrixMult(tempL, tempR, b_T, m6);

			SubMatrix.Subtract(a12, a22, tempL);
			SubMatrix.Add(b21, b22, tempR);
			MatrixMult(tempL, tempR, b_T, m7);


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
			int offsetR = 0;
			for (int i = 0; i < childLength; i++)
			{
				if (resultRowCount <= i)
					break;

				//int offsetR = i * resultColumnCount;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					if (resultColumnCount <= j)
						break;

					result[offsetR + j] = c11.RealMatrix[offsetC + j];
				}

				offsetR += resultColumnCount;
			}

			offsetR = childLength;
			for (int i = 0; i < childLength; i++)
			{
				if (i >= resultRowCount)
					break;

				//int offsetR = (i * resultColumnCount) + childLength;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					if (j + childLength >= resultColumnCount)
						break;

					result[offsetR + j] = c12.RealMatrix[offsetC + j];
				}
				offsetR += childLength;
			}

			offsetR = childLength * resultColumnCount;
			for (int i = 0; i < childLength; i++)
			{
				if (i + childLength >= resultRowCount)
					break;

				//int offsetR = (i + childLength) * resultColumnCount;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					if (j >= resultColumnCount)
						break;

					result[offsetR + j] = c21.RealMatrix[offsetC + j];
				}
				offsetR += resultColumnCount;
			}

			offsetR = childLength * resultColumnCount + childLength;
			for (int i = 0; i < childLength; i++)
			{
				if (i + childLength >= resultRowCount)
					break;

				//int offsetR = ((i + childLength) * resultColumnCount) + childLength;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					if (j + childLength >= resultColumnCount)
						break;

					result[offsetR + j] = c22.RealMatrix[offsetC + j];
				}
				offsetR += resultColumnCount;
			}


			///
			/// Return the result
			/// 
			return result;
		}

		private static void MatrixMult(SubMatrix leftSide, SubMatrix rightSide, SubMatrix b_T, SubMatrix result)
		{
			if (leftSide.Length < SwitchLength)
			{
				SubMatrix.MatrixMult(leftSide, rightSide, b_T, result);
				return;
			}

			var childLength		= leftSide.Length / 2;
			var childBlockSize	= childLength * childLength;

			///
			/// Break A, B, and C down into 4 components (each).
			/// 
			var a11		= new SubMatrix() { RealMatrix = leftSide.RealMatrix,	Length = childLength, RealRowWidth = leftSide.RealRowWidth,		StartX = leftSide.StartX,					StartY = leftSide.StartY };
			var a12		= new SubMatrix() { RealMatrix = leftSide.RealMatrix,	Length = childLength, RealRowWidth = leftSide.RealRowWidth,		StartX = leftSide.StartX + childLength,		StartY = leftSide.StartY };
			var a21		= new SubMatrix() { RealMatrix = leftSide.RealMatrix,	Length = childLength, RealRowWidth = leftSide.RealRowWidth,		StartX = leftSide.StartX,					StartY = leftSide.StartY + childLength };
			var a22		= new SubMatrix() { RealMatrix = leftSide.RealMatrix,	Length = childLength, RealRowWidth = leftSide.RealRowWidth,		StartX = leftSide.StartX + childLength,		StartY = leftSide.StartY + childLength };

			var b11		= new SubMatrix() { RealMatrix = rightSide.RealMatrix,	Length = childLength, RealRowWidth = rightSide.RealRowWidth,	StartX = rightSide.StartX,					StartY = rightSide.StartY };
			var b12		= new SubMatrix() { RealMatrix = rightSide.RealMatrix,	Length = childLength, RealRowWidth = rightSide.RealRowWidth,	StartX = rightSide.StartX + childLength,	StartY = rightSide.StartY };
			var b21		= new SubMatrix() { RealMatrix = rightSide.RealMatrix,	Length = childLength, RealRowWidth = rightSide.RealRowWidth,	StartX = rightSide.StartX,					StartY = rightSide.StartY + childLength };
			var b22		= new SubMatrix() { RealMatrix = rightSide.RealMatrix,	Length = childLength, RealRowWidth = rightSide.RealRowWidth,	StartX = rightSide.StartX + childLength,	StartY = rightSide.StartY + childLength };

			var c11		= new SubMatrix() { RealMatrix = result.RealMatrix,		Length = childLength, RealRowWidth = result.RealRowWidth,		StartX = result.StartX,						StartY = result.StartY };
			var c12		= new SubMatrix() { RealMatrix = result.RealMatrix,		Length = childLength, RealRowWidth = result.RealRowWidth,		StartX = result.StartX + childLength,		StartY = result.StartY };
			var c21		= new SubMatrix() { RealMatrix = result.RealMatrix,		Length = childLength, RealRowWidth = result.RealRowWidth,		StartX = result.StartX,						StartY = result.StartY + childLength };
			var c22		= new SubMatrix() { RealMatrix = result.RealMatrix,		Length = childLength, RealRowWidth = result.RealRowWidth,		StartX = result.StartX + childLength,		StartY = result.StartY + childLength };

			///
			/// Create new matrices for our intermidiate values.
			/// 
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
			/// Calculte m1 to m7
			/// 
			SubMatrix.Add(a11, a22, tempL);
			SubMatrix.Add(b11, b22, tempR);
			MatrixMult(tempL, tempR, b_T, m1);

			SubMatrix.Add(a21, a22, tempL);
			MatrixMult(tempL, b11, b_T, m2);

			SubMatrix.Subtract(b12, b22, tempR);
			MatrixMult(a11, tempR, b_T, m3);

			SubMatrix.Subtract(b21, b11, tempR);
			MatrixMult(a22, tempR, b_T, m4);

			SubMatrix.Add(a11, a12, tempL);
			MatrixMult(tempL, b22, b_T, m5);

			SubMatrix.Subtract(a21, a11, tempL);
			SubMatrix.Add(b11, b12, tempR);
			MatrixMult(tempL, tempR, b_T, m6);

			SubMatrix.Subtract(a12, a22, tempL);
			SubMatrix.Add(b21, b22, tempR);
			MatrixMult(tempL, tempR, b_T, m7);


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
			/// We don't need to copy { c11, c12, c21, c22 } back over to the result because they are already in 
			/// the result.
			/// 
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

				int sourceOffset = sourceRow * sourceWidth;
				int resultOffset = i * length;

				for (int j = 0; j < length; j++)
				{
					var sourceColumn = j + startColumn;
					if (sourceWidth <= sourceColumn)
						break;

					result[resultOffset + j] = source[sourceOffset + sourceColumn];
				}
			}
		}

	}
}
