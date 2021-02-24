
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

		private static float[] A11;
		private static float[] A12;
		private static float[] A21;
		private static float[] A22;

		private static float[] B11;
		private static float[] B12;
		private static float[] B21;
		private static float[] B22;


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

			A11 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			A12 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			A21 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			A22 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);

			B11 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			B12 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			B21 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);
			B22 = MonoMatrixOperations.CreateZeroMatrix(childBlockSize);

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

	}
}
