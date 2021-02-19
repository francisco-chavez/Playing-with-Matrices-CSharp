
using System;
using System.Numerics;


namespace MainProgram
{
	public static class DoubleCheck
	{

		public static double[] MatrixMult(float[] matrixA, float[] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			if (sizeA.Item2 != sizeB.Item1)
				throw new Exception();

			var resultRowCount      = sizeA.Item1;
			var resultColumnCount   = sizeB.Item2;
			var sharedDimension     = sizeA.Item2;

			var blockSizeA  = sizeA.Item1 * sizeA.Item2;
			var aPrime      = new double[blockSizeA];
			for (int i = 0; i < blockSizeA; i++)
				aPrime[i] = matrixA[i];


			var matrixB_T = new double[sizeB.Item1 * sizeB.Item2];
			for (int i = 0; i < sizeB.Item1; i++)
			{
				var sourceOffset = i * sizeB.Item2;

				for (int j = 0; j < sizeB.Item2; j++)
					matrixB_T[j * sizeB.Item1 + i] = matrixB[sourceOffset + j];
			}

			var resultSize = resultRowCount * resultColumnCount;
			var result      = new double[resultSize];

			int rI          = 0;
			int aRowOffset  = 0;

			var elementCount = Vector<double>.Count;
			var simdCount = sharedDimension / elementCount;

			for (int rY = 0; rY < resultRowCount; rY++)
			{
				int bRowOffset = 0;
				for (int rX = 0; rX < resultColumnCount; rX++)
				{
					result[rI] = 0.0;

					int i_sd = 0;
					for (int i_e = 0; i_e < simdCount; i_e++)
					{
						var vectorA = new Vector<double>(aPrime, aRowOffset + i_sd);
						var vectorB = new Vector<double>(matrixB_T, bRowOffset + i_sd);

						result[rI] += Vector.Dot(vectorA, vectorB);

						i_sd += elementCount;
					}
					for (; i_sd < sharedDimension; i_sd++)
						result[rI] += matrixA[aRowOffset + i_sd] * matrixB_T[bRowOffset + i_sd];
					rI++;
					bRowOffset += sharedDimension;
				}
				aRowOffset += sharedDimension;
			}

			return result;
		}
	}
}
