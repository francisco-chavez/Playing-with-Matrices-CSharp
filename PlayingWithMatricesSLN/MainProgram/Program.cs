
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

			//int rWidth		= 1;
			//int rHeight		= 1;
			//int dotLength	= 3;
			//
			//var matrixA		= new float[] { 1.0f, 2.0f, 3.0f };
			//var matrixB		= new float[] { 4.0f, 5.0f, 6.0f };

			int rHeight		=  9000;
			int dotLength	= 10000;
			int rWidth		=  8000;

			var matrixA = MonoMatrixOperations.CreateRandomMatrix(0, 
																  height: rHeight, 
																  width: dotLength, 
																  minValue: -1.0f, 
																  maxValue: +4.0f);

			var matrixB	= MonoMatrixOperations.CreateRandomMatrix(5, 
																  height: dotLength, 
																  width: rWidth, 
																  minValue: -1.0f, 
																  maxValue: +3.0f);

			var sizeA		= new Tuple<int, int>(rHeight, dotLength);
			var sizeB		= new Tuple<int, int>(dotLength, rWidth);



			DateTime t0;
			DateTime t1;
			TimeSpan monoConventional;
			TimeSpan monoDotTTime;
			TimeSpan monoDotTSIMDTime;
			TimeSpan monoStrassenTime;
			TimeSpan newStrassenTime;


			GC.Collect();

			//t0 = DateTime.Now;
			//var oldSchool			= MonoMatrixOperations.MatrixMult_Conventional(matrixA, matrixB, sizeA, sizeB);
			//GC.Collect();
			//t1 = DateTime.Now;
			//monoConventional		= t1 - t0;
			
			t0 = DateTime.Now;
			var dotTResult			= MonoMatrixOperations.MatrixMult_TransposeDotProduct(matrixA, matrixB, sizeA, sizeB);
			GC.Collect();
			t1 = DateTime.Now;
			monoDotTTime			= t1 - t0;
			
			t0 = DateTime.Now;
			var dotTsimdResult		= MonoMatrixOperations.MatrixMult_TransposeDotProductVector(matrixA, matrixB, sizeA, sizeB);
			GC.Collect();
			t1 = DateTime.Now;
			monoDotTSIMDTime		= t1 - t0;
			
			t0 = DateTime.Now;
			var strassenResult		= MonoMatrixOperations.MatrixMult_Strassen(matrixA, matrixB, sizeA, sizeB);
			GC.Collect();
			t1 = DateTime.Now;
			monoStrassenTime		= t1 - t0;

			t0 = DateTime.Now;
			var strassenV2Result	= StrassenMult.MatrixMult(matrixA, matrixB, sizeA, sizeB);
			GC.Collect();
			t1 = DateTime.Now;
			newStrassenTime			= t1 - t0;



			// var valueCheck = DoubleCheck.MatrixMult(matrixA, matrixB, sizeA, sizeB);
			// var resultSize = rWidth * rHeight;
			// 
			// var dotTError			= 0.0;
			// var strassenError		= 0.0;
			// var strassen2Error		= 0.0;
			// 
			// var maxDotTError		= 0.0;
			// var maxStrassenError	= 0.0;
			// var maxStrassen2Error	= 0.0;
			// 
			// for (int i = 0; i < resultSize; i++)
			// {
			// 	var error = Math.Abs(valueCheck[i] - dotTsimdResult[i]);
			// 	dotTError += error;
			// 	maxDotTError = Math.Max(maxDotTError, error);
			// 
			// 	error = Math.Abs(valueCheck[i] - strassenResult[i]);
			// 	strassenError += error;
			// 	maxStrassenError = Math.Max(maxStrassenError, error);
			// 
			// 	error = Math.Abs(valueCheck[i] - strassenV2Result[i]);
			// 	strassen2Error += error;
			// 	maxStrassen2Error = Math.Max(maxStrassen2Error, error);
			// }
			// 
			// int breakPoint = 0;
			// breakPoint++;


			//Console.WriteLine("[A] * [B] Conventional:          {0}s", monoConventional.TotalSeconds);
			Console.WriteLine("[A] * [B] Dot Product:           {0}s", monoDotTTime.TotalSeconds);
			Console.WriteLine("[A] * [B] Dot Product with SIMD: {0}s", monoDotTSIMDTime.TotalSeconds);
			Console.WriteLine("[A] * [B] Strassen:              {0}s", monoStrassenTime.TotalSeconds);
			Console.WriteLine("[A] * [B] Strassen 2:            {0}s", newStrassenTime.TotalSeconds);

			string userInput;
			Console.WriteLine();
			Console.Write("Press Enter twice to Exit.");
			Console.ReadLine();
			userInput = Console.ReadLine();
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
