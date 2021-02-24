
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

			var a11 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			var a12 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			var a21 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			var a22 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);

			var b11 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			var b12 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			var b21 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			var b22 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);

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
