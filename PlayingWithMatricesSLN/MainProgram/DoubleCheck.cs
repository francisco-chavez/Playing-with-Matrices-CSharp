
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace MainProgram
{
	public static class DoubleCheck
	{
		public static double[] MatrixMult_TransposeDotProduct(float[] matrixA, float[] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			if (sizeA.Item2 != sizeB.Item1)
				throw new Exception( );

			var resultRowCount		= sizeA.Item1;
			var resultColumnCount	= sizeB.Item2;
			var sharedDimension		= sizeA.Item2;

			var blockSizeA	= sizeA.Item1 * sizeA.Item2;
			var aPrime		= new double[blockSizeA];
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
			var result		= new double[resultSize];

			int rI			= 0;
			int aRowOffset	= 0;

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
						result[rI] += aPrime[aRowOffset + sD] * matrixB_T[bRowOffset + sD];
					}
					rI++;
					bRowOffset += sharedDimension;
				}
				aRowOffset += sharedDimension;
			}

			return result;
		}

	}
}
