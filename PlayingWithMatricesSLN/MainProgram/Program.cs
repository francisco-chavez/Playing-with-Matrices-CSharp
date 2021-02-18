
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace MainProgram
{
	class Program
	{
		static void Main(string[] args)
		{
			// int rWidth		= 1500;
			// int rHeight		= 1500;
			// int dotLength	= 1500;
			// 
			// var matrixA		= CreateRandomMonoMatrix(seed: 0, height: rHeight, width: dotLength);
			// var matrixB		= CreateRandomMonoMatrix(seed: 1, height: dotLength, width: rWidth);

			//int rWidth		= 1;
			//int rHeight		= 1;
			//int dotLength	= 3;
			//
			//var matrixA		= new float[] { 1.0f, 2.0f, 3.0f };
			//var matrixB		= new float[] { 4.0f, 5.0f, 6.0f };


			int rWidth		= 40;
			int rHeight		= 1155;
			int dotLength	= 800;

			var matrixA = MonoMatrixOperations.CreateRandomMatrix(0, 
																  height: rHeight, 
																  width: dotLength, 
																  minValue: -1.0f, 
																  maxValue: +1.0f);

			var matrixB	= MonoMatrixOperations.CreateRandomMatrix(5, 
																  height: dotLength, 
																  width: rWidth, 
																  minValue: -1.0f, 
																  maxValue: +2.0f);

			var sizeA		= new Tuple<int, int>(rHeight, dotLength);
			var sizeB		= new Tuple<int, int>(dotLength, rWidth);


			var valueCheck = DoubleCheck.MatrixMult_TransposeDotProduct(matrixA, matrixB, sizeA, sizeB);

			var conventionalResult	= MonoMatrixOperations.MatrixMult_Conventional(matrixA, matrixB, sizeA, sizeB);
			//var conventionalResult = MatrixMult_MonoArray_Conventional(matrixB, matrixA, sizeB, sizeA);
			var dotTResult			= MonoMatrixOperations.MatrixMult_TransposeDotProduct(matrixA, matrixB, sizeA, sizeB);
			//var dotTResult		= MatrixMult_MonoArray_TransposeDotProduct(matrixB, matrixA, sizeB, sizeA);
			var strassenResult		= MonoMatrixOperations.MatrixMult_Strassen(matrixA, matrixB, sizeA, sizeB);
			//var strassenResult	= MatrixMult_MonoArray_Strassen(matrixB, matrixA, sizeB, sizeA);

			var resultSize = rWidth * rHeight;

			double conventionalError = 0.0;
			double dotTError		= 0.0;
			double strassenError	= 0.0;

			for (int i = 0; i < resultSize; i++)
			{
				conventionalError += Math.Abs(valueCheck[i] - conventionalResult[i]);
				dotTError += Math.Abs(valueCheck[i] - dotTResult[i]);
				strassenError += Math.Abs(valueCheck[i] - strassenResult[i]);

				var error_c_i = conventionalResult[i] - valueCheck[i];
				var error_s_i = strassenResult[i] - valueCheck[i];

				var breakPointB = 0;
				breakPointB++;
			}

			int breakPoint = 0;
			breakPoint++;


			//DateTime t0;
			//DateTime t1;
			//TimeSpan monoConventional;
			//TimeSpan monoDotT;
			//
			//t0 = DateTime.Now;
			////var conventionalResult = MatrixMult_MonoArray_Conventional(matrixA, matrixB, sizeA, sizeB);
			//var conventionalResult = MatrixMult_MonoArray_Conventional(matrixB, matrixA, sizeB, sizeA);
			//t1 = DateTime.Now;
			//
			//monoConventional = t1 - t0;
			//
			//t0 = DateTime.Now;
			////var dotTResult = MatrixMult_MonoArray_TransposeDotProduct(matrixA, matrixB, sizeA, sizeB);
			//var dotTResult = MatrixMult_MonoArray_TransposeDotProduct(matrixB, matrixA, sizeB, sizeA);
			//t1 = DateTime.Now;
			//
			//monoDotT = t1 - t0;
			//
			//Console.WriteLine("[A] * [B] Conventional: {0}s", monoConventional.TotalSeconds);
			//Console.WriteLine("[A] * [B] Dot Product: {0}s", monoDotT.TotalSeconds);

			Console.WriteLine();
			Console.Write("Press any key to Exit.");
			var userInput = Console.ReadKey(true);
		}





		public static float[,] CreateRandomNeatMatrix(int seed, int width, int height, float minValue = 0.0f, float maxValue = 1.0f)
		{
			Random r = new Random(seed);

			double range = (double) maxValue;
			range -= minValue;

			var matrix = new float[height, width];

			// I'm using doubles to create the initial values to ensure the greatest range in bits.
			for (int h = 0; h < height; h++)
				for (int w = 0; w < width; w++)
					matrix[h, w] = (float) ((range * r.NextDouble()) + minValue);

			return matrix;
		}

		public static float[,] MatrixMult_NeatArray_Conventional(float[,] matrixA, float[,] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			if (sizeA.Item2 != sizeB.Item1)
				throw new Exception();

			var resultRowCount		= sizeA.Item1;
			var resultColumnCount	= sizeB.Item2;
			var sharedDimension		= sizeA.Item2;

			var result = new float[resultRowCount, resultColumnCount];

			for (int rY = 0; rY < resultRowCount; rY++)
				for (int rX = 0; rX < resultColumnCount; rX++)
				{
					result[rY, rX] = 0.0f;
					for (int sD = 0; sD < sharedDimension; sD++)
						result[rY, rX] += matrixA[rY, sD] * matrixB[sD, rX];
				}

			return result;
		}

	}
}
