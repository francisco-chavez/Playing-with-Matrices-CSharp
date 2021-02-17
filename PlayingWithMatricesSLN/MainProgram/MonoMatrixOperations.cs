﻿
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace MainProgram
{

	public static class MonoMatrixOperations
	{
		public static float[] CreateRandomMonoMatrix(int seed, int width, int height, float minValue = 0.0f, float maxValue = 1.0f)
		{
			Random r = new Random(seed);
			int size = width * height;

			double range = (double) maxValue;
			range -= minValue;

			var matrix = new float[size];

			// I'm using doubles to create the initial values to ensure the greatest range in bits.
			for (int i = 0; i < size; i++)
				matrix[i] = (float) ((range * r.NextDouble()) + minValue);

			return matrix;
		}

		public static float[] MatrixMult_Conventional(float[] matrixA, float[] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			if (sizeA.Item2 != sizeB.Item1)
				throw new Exception();

			var resultRowCount = sizeA.Item1;
			var resultColumnCount = sizeB.Item2;
			var sharedDimension = sizeA.Item2;

			var resultSize = resultRowCount * resultColumnCount;
			var result = new float[resultSize];

			int rI = 0;
			int aRowOffset = 0;

			for (int rY = 0; rY < resultRowCount; rY++)
			{
				for (int rX = 0; rX < resultColumnCount; rX++)
				{
					//result[rY * resultColumnCount + rX] = 0.0f;
					result[rI] = 0.0f;
					for (int sD = 0; sD < sharedDimension; sD++)
					{
						//result[rY * resultColumnCount + rX] += matrixA[rY * sizeA.Item2 + sD] * matrixB[sD * sizeB.Item2 + rX];
						result[rI] += matrixA[aRowOffset + sD] * matrixB[sD * sizeB.Item2 + rX];
					}
					rI++;
				}
				aRowOffset += sizeA.Item2;
			}

			return result;
		}

		public static float[] MatrixMult_TransposeDotProduct(float[] matrixA, float[] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			if (sizeA.Item2 != sizeB.Item1)
				throw new Exception();

			var resultRowCount = sizeA.Item1;
			var resultColumnCount = sizeB.Item2;
			var sharedDimension = sizeA.Item2;


			var matrixB_T = new float[sizeB.Item1 * sizeB.Item2];
			for (int i = 0; i < sizeB.Item1; i++)
				for (int j = 0; j < sizeB.Item2; j++)
					matrixB_T[j * sizeB.Item1 + i] = matrixB[i * sizeB.Item2 + j];

			var resultSize = resultRowCount * resultColumnCount;
			var result = new float[resultSize];

			int rI = 0;
			int aRowOffset = 0;

			for (int rY = 0; rY < resultRowCount; rY++)
			{
				int bRowOffset = 0;
				for (int rX = 0; rX < resultColumnCount; rX++)
				{
					//result[rY * resultColumnCount + rX] = 0.0f;
					result[rI] = 0.0f;
					for (int sD = 0; sD < sharedDimension; sD++)
					{
						//result[rY * resultColumnCount + rX] += matrixA[rY * sharedDimension + sD] * matrixB_T[rX * sharedDimension + sD];
						result[rI] += matrixA[aRowOffset + sD] * matrixB_T[bRowOffset + sD];
					}
					rI++;
					bRowOffset += sharedDimension;
				}
				aRowOffset += sharedDimension;
			}

			return result;
		}

		public static float[] MatrixMult_Strassen(float[] matrixA, float[] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			if (sizeA.Item2 != sizeB.Item1)
				throw new Exception();

			var resultRowCount = sizeA.Item1;
			var resultColumnCount = sizeB.Item2;
			var sharedDimension = sizeA.Item2;

			var maxD = sizeA.Item1 > sizeA.Item2 ? sizeA.Item2 : sizeA.Item1;
			maxD = maxD < sizeB.Item2 ? sizeB.Item2 : maxD;

			int n = 0;
			int t = 1;

			while (t < maxD)
			{
				t *= 2;
				n++;
			}

			// If the matrices are below a certain size, then using the Strassen algorithm isn't worth it. Also, 
			// the dot transpose method is there for avoid cache misses, but if the matrix is small enough, then 
			// that's not really a worry either.
			if (n < 5)
			{
				return MatrixMult_Conventional(matrixA, matrixB, sizeA, sizeB);
			}

			var halfT = t / 2;

			var start_11_0 = 0;
			var start_11_1 = 0;
			//var end_11_0 = (halfT) - 1;
			//var end_11_1 = (halfT) - 1;

			var start_12_0 = 0;
			var start_12_1 = halfT;
			//var end_12_0 = (halfT) - 1;
			//var end_12_1 = t - 1;

			var start_21_0 = halfT;
			var start_21_1 = 0;
			//var end_21_0 = t - 1;
			//var end_21_1 = (halfT) - 1;

			var start_22_0 = halfT;
			var start_22_1 = halfT;
			//var end_22_0 = t - 1;
			//var end_22_1 = t - 1;


			int blockSize = halfT * halfT;

			var a11 = new float[blockSize];
			var a12 = new float[blockSize];
			var a21 = new float[blockSize];
			var a22 = new float[blockSize];

			var b11 = new float[blockSize];
			var b12 = new float[blockSize];
			var b21 = new float[blockSize];
			var b22 = new float[blockSize];

			FillInBlock(a11, matrixA, start_11_1, start_11_0, halfT, sizeA.Item2, sizeA.Item1);
			FillInBlock(a12, matrixA, start_12_1, start_12_0, halfT, sizeA.Item2, sizeA.Item1);
			FillInBlock(a21, matrixA, start_21_1, start_21_0, halfT, sizeA.Item2, sizeA.Item1);
			FillInBlock(a22, matrixA, start_22_1, start_22_0, halfT, sizeA.Item2, sizeA.Item1);

			FillInBlock(b11, matrixB, start_11_1, start_11_0, halfT, sizeB.Item2, sizeB.Item1);
			FillInBlock(b12, matrixB, start_12_1, start_12_0, halfT, sizeB.Item2, sizeB.Item1);
			FillInBlock(b21, matrixB, start_21_1, start_21_0, halfT, sizeB.Item2, sizeB.Item1);
			FillInBlock(b22, matrixB, start_22_1, start_22_0, halfT, sizeB.Item2, sizeB.Item1);


			var m1 = new float[blockSize];
			var m2 = new float[blockSize];
			var m3 = new float[blockSize];
			var m4 = new float[blockSize];
			var m5 = new float[blockSize];
			var m6 = new float[blockSize];
			var m7 = new float[blockSize];

			var tm1 = new float[blockSize];
			var tm2 = new float[blockSize];

			///
			/// Calculte m1 to m7
			/// 
			Add(a11, a22, tm1, blockSize);
			Add(b11, b22, tm2, blockSize);
			MatrixMult_StrassenRecursiveComp(tm1, tm2, m1, halfT);


			Add(a21, a22, tm1, blockSize);
			MatrixMult_StrassenRecursiveComp(tm1, b11, m2, halfT);


			Subtract(b12, b22, tm2, blockSize);
			MatrixMult_StrassenRecursiveComp(a11, tm2, m3, halfT);

			Subtract(b21, b11, tm2, blockSize);
			MatrixMult_StrassenRecursiveComp(a22, tm2, m4, halfT);

			Add(a11, a12, tm1, blockSize);
			MatrixMult_StrassenRecursiveComp(tm1, b22, m5, halfT);

			Subtract(a21, a11, tm1, blockSize);
			Add(b11, b12, tm2, blockSize);
			MatrixMult_StrassenRecursiveComp(tm1, tm2, m6, halfT);

			Subtract(a12, a22, tm1, blockSize);
			Add(b21, b22, tm2, blockSize);
			MatrixMult_StrassenRecursiveComp(tm1, tm2, m7, halfT);


			float[] c11 = new float[blockSize];
			float[] c12 = new float[blockSize];
			float[] c21 = new float[blockSize];
			float[] c22 = new float[blockSize];

			///
			/// Calculate c11, c12, c21, c22
			/// 
			Add(m1, m4, tm1, blockSize);
			Subtract(m7, m5, tm2, blockSize);
			Add(tm1, tm2, c11, blockSize);

			Add(m3, m5, c12, blockSize);

			Add(m2, m4, c21, blockSize);

			Subtract(m1, m2, tm1, blockSize);
			Add(m3, m6, tm2, blockSize);
			Add(tm1, tm2, c22, blockSize);

			float[] result = new float[resultRowCount * resultColumnCount];

			///
			/// Transfer { c11, c12, c21, c22 } to result
			/// 
			for (int i = 0; i < halfT; i++)
			{
				if (i >= resultRowCount)
					break;

				int offsetR = i * resultColumnCount;
				int offsetC = i * halfT;
				for (int j = 0; j < halfT; j++)
				{
					if (j >= resultColumnCount)
						break;

					result[offsetR + j] = c11[offsetC + j];
				}
			}

			for (int i = 0; i < halfT; i++)
			{
				if (i >= resultRowCount)
					break;

				int offsetR = (i * resultColumnCount) + halfT;
				int offsetC = i * halfT;

				for (int j = 0; j < halfT; j++)
				{
					if (j + halfT >= resultColumnCount)
						break;

					result[offsetR + j] = c12[offsetC + j];
				}
			}

			for (int i = 0; i < halfT; i++)
			{
				if (i + halfT >= resultRowCount)
					break;

				int offsetR = (i + halfT) * resultColumnCount;
				int offsetC = i * halfT;

				for (int j = 0; j < halfT; j++)
				{
					if (j >= resultColumnCount)
						break;

					result[offsetR + j] = c21[offsetC + j];
				}
			}

			for (int i = 0; i < halfT; i++)
			{
				if (i + halfT >= resultRowCount)
					break;

				int offsetR = ((i + halfT) * resultColumnCount) + halfT;
				int offsetC = i * halfT;

				for (int j = 0; j < halfT; j++)
				{
					if (j + halfT >= resultColumnCount)
						break;

					result[offsetR + j] = c22[offsetC + j];
				}
			}

			return result;
		}


		private static void MatrixMult_StrassenRecursiveComp(float[] matrixA, float[] matrixB, float[] result, int length)
		{
			var blockSize = length * length;

			// So, once an matrix falls below a certain size, we no longer need to worry so much about cache misses 
			// and such. At that point, we can just use the conventional algorithm.
			if (length < 3)
			{
				var tempResult = MatrixMult_Conventional(matrixA, matrixB, new Tuple<int, int>(length, length), new Tuple<int, int>(length, length));
				Array.Copy(tempResult, result, blockSize);
				return;
			}

			var childLength = length / 2;
			var childBlockSize = childLength * childLength;

			var a11 = new float[childBlockSize];
			var a12 = new float[childBlockSize];
			var a21 = new float[childBlockSize];
			var a22 = new float[childBlockSize];

			var b11 = new float[childBlockSize];
			var b12 = new float[childBlockSize];
			var b21 = new float[childBlockSize];
			var b22 = new float[childBlockSize];

			//			dest,	source,		start column	start row,		length along row/column
			FillInBlock(a11, matrixA, 0, 0, childLength);
			FillInBlock(a12, matrixA, childLength, 0, childLength);
			FillInBlock(a21, matrixA, 0, childLength, childLength);
			FillInBlock(a22, matrixA, childLength, childLength, childLength);

			FillInBlock(b11, matrixB, 0, 0, childLength);
			FillInBlock(b12, matrixB, childLength, 0, childLength);
			FillInBlock(b21, matrixB, 0, childLength, childLength);
			FillInBlock(b22, matrixB, childLength, childLength, childLength);


			var m1 = new float[childBlockSize];
			var m2 = new float[childBlockSize];
			var m3 = new float[childBlockSize];
			var m4 = new float[childBlockSize];
			var m5 = new float[childBlockSize];
			var m6 = new float[childBlockSize];
			var m7 = new float[childBlockSize];

			var tm1 = new float[childBlockSize];
			var tm2 = new float[childBlockSize];


			///
			/// Calculte m1 to m7
			/// 
			Add(a11, a22, tm1, childBlockSize);
			Add(b11, b22, tm2, childBlockSize);
			MatrixMult_StrassenRecursiveComp(tm1, tm2, m1, childLength);


			Add(a21, a22, tm1, childBlockSize);
			MatrixMult_StrassenRecursiveComp(tm1, b11, m2, childLength);


			Subtract(b12, b22, tm2, childBlockSize);
			MatrixMult_StrassenRecursiveComp(a11, tm2, m3, childLength);

			Subtract(b21, b11, tm2, childBlockSize);
			MatrixMult_StrassenRecursiveComp(a22, tm2, m4, childLength);

			Add(a11, a12, tm1, childBlockSize);
			MatrixMult_StrassenRecursiveComp(tm1, b22, m5, childLength);

			Subtract(a21, a11, tm1, childBlockSize);
			Add(b11, b12, tm2, childBlockSize);
			MatrixMult_StrassenRecursiveComp(tm1, tm2, m6, childLength);

			Subtract(a12, a22, tm1, childBlockSize);
			Add(b21, b22, tm2, childBlockSize);
			MatrixMult_StrassenRecursiveComp(tm1, tm2, m7, childLength);


			float[] c11 = new float[childBlockSize];
			float[] c12 = new float[childBlockSize];
			float[] c21 = new float[childBlockSize];
			float[] c22 = new float[childBlockSize];

			///
			/// Calculate c11, c12, c21, c22
			/// 
			Add(m1, m4, tm1, childBlockSize);
			Subtract(m7, m5, tm2, childBlockSize);
			Add(tm1, tm2, c11, childBlockSize);

			Add(m3, m5, c12, childBlockSize);

			Add(m2, m4, c21, childBlockSize);

			Subtract(m1, m2, tm1, childBlockSize);
			Add(m3, m6, tm2, childBlockSize);
			Add(tm1, tm2, c22, childBlockSize);

			///
			/// Transfer { c11, c12, c21, c22 } to result
			/// 
			for (int i = 0; i < childLength; i++)
			{
				int offsetR = i * length;
				int offsetC = i * childLength;
				for (int j = 0; j < childLength; j++)
				{
					result[offsetR + j] = c11[offsetC + j];
				}
			}

			for (int i = 0; i < childLength; i++)
			{
				int offsetR = (i * length) + childLength;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					result[offsetR + j] = c12[offsetC + j];
				}
			}

			for (int i = 0; i < childLength; i++)
			{

				int offsetR = (i + childLength) * length;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{

					result[offsetR + j] = c21[offsetC + j];
				}
			}

			for (int i = 0; i < childLength; i++)
			{
				int offsetR = ((i + childLength) * length) + childLength;
				int offsetC = i * childLength;

				for (int j = 0; j < childLength; j++)
				{
					result[offsetR + j] = c22[offsetC + j];
				}
			}

			//int foo = 0;
			//foo++;
		}

		private static void FillInBlock(float[] result, float[] source, int startColumn, int startRow, int length, int sourceWidth, int sourceHeight)
		{
			for (int i = 0; i < length; i++)
			{
				int currentRowSource = i + startRow;

				for (int j = 0; j < length; j++)
				{
					float value;

					var indexR = (i * length) + j;

					int currentColumnSource = j + startColumn;

					if (currentColumnSource < sourceWidth && currentRowSource < sourceHeight)
					{
						var indexS = (currentRowSource * sourceWidth) + currentColumnSource;
						value = source[indexS];
					}
					else
					{
						value = 0.0f;
					}

					result[indexR] = value;
				}
			}
		}

		/// <summary>
		/// This will copy a section of <paramref name="source"/> into <paramref name="result"/>.
		/// </summary>
		/// <param name="result">This is a matrix that's <paramref name="length"/> by <paramref name="length"/> in size. We will copy our data into this matrix.</param>
		/// <param name="source">This is a matrix that's 2 * <paramref name="length"/> by 2 * <paramref name="length"/> in size. This is the matrix we will by copying from.</param>
		/// <param name="startColumn">The column in the source matrix that we start coping data from.</param>
		/// <param name="startRow">The row in the source matrix that we start coping data from.</param>
		/// <param name="length">The dimensional size reference used in our matrices.</param>
		public static void FillInBlock(float[] result, float[] source, int startColumn, int startRow, int resultLength)
		{
			for (int i = 0; i < resultLength; i++)
			{
				int currentRowSource = i + startRow;

				int offsetSource = (i + startRow) * (resultLength * 2) + startColumn;
				int offsetResult = i * resultLength;

				for (int j = 0; j < resultLength; j++)
				{
					var indexR = offsetResult + j;
					var indexS = offsetSource + j;

					result[indexR] = source[indexS];
				}
			}
		}

		public static void Add(float[] left, float[] right, float[] result, int length)
		{
			for (int i = 0; i < length; i++)
				result[i] = left[i] + right[i];
		}

		public static void Subtract(float[] left, float[] right, float[] result, int length)
		{
			for (int i = 0; i < length; i++)
				result[i] = left[i] - right[i];
		}



	}

}