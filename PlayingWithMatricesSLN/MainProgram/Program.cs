
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
			int rWidth		= 1500;
			int rHeight		= 1500;
			int dotLength	= 1500;

			var matrixA		= CreateRandomMonoMatrix(seed: 0, height: rHeight, width: dotLength);
			var matrixB		= CreateRandomMonoMatrix(seed: 1, height: dotLength, width: rWidth);

			var sizeA		= new Tuple<int, int>(rHeight, dotLength);
			var sizeB		= new Tuple<int, int>(dotLength, rWidth);

			// int rWidth		= 1;
			// int rHeight		= 1;
			// int dotLength	= 3;
			   
			// var matrixA		= new float[] { 1.0f, 2.0f, 3.0f };
			// var matrixB		= new float[] { 4.0f, 5.0f, 6.0f };
			   
			// var sizeA		= new Tuple<int, int>(rHeight, dotLength);
			// var sizeB		= new Tuple<int, int>(dotLength, rWidth);

			DateTime t0;
			DateTime t1;
			TimeSpan monoConventional;
			TimeSpan monoDotT;

			t0 = DateTime.Now;
			var contentionalResult = MatrixMult_MonoArray_Conventional(matrixA, matrixB, sizeA, sizeB);
			t1 = DateTime.Now;

			monoConventional = t1 - t0;

			t0 = DateTime.Now;
			var dotTResult = MatrixMult_MonoArray_TransposeDotProduct(matrixA, matrixB, sizeA, sizeB);
			t1 = DateTime.Now;

			monoDotT = t1 - t0;

			Console.WriteLine("[A] * [B] Conventional: {0}s", monoConventional.TotalSeconds);
			Console.WriteLine("[A] * [B] Dot Product: {0}s", monoDotT.TotalSeconds);

			Console.WriteLine();
			Console.Write("Press any key to Exit.");
			var userInput = Console.ReadKey(true);
		}

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

		public static float[] MatrixMult_MonoArray_Conventional(float[] matrixA, float[] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			if (sizeA.Item2 != sizeB.Item1)
				throw new Exception();

			var resultRowCount		= sizeA.Item1;
			var resultColumnCount	= sizeB.Item2;
			var sharedDimension		= sizeA.Item2;

			var resultSize = resultRowCount * resultColumnCount;
			var result = new float[resultSize];
			for (int i = 0; i < resultSize; i++)
				result[i] = 0.0f;

			for (int rY = 0; rY < resultRowCount; rY++)
				for (int rX = 0; rX < resultColumnCount; rX++)
					for (int sD = 0; sD < sharedDimension; sD++)
						result[rY * resultColumnCount + rX] += matrixA[rY * sizeA.Item2 + sD] * matrixB[sD * sizeB.Item2 + rX];

			return result;
		}

		public static float[] MatrixMult_MonoArray_TransposeDotProduct(float[] matrixA, float[] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			if (sizeA.Item2 != sizeB.Item1)
				throw new Exception();

			var resultRowCount		= sizeA.Item1;
			var resultColumnCount	= sizeB.Item2;
			var sharedDimension		= sizeA.Item2;


			var matrixB_T = new float[sizeB.Item1 * sizeB.Item2];
			for (int i = 0; i < sizeB.Item1; i++)
				for (int j = 0; j < sizeB.Item2; j++)
					matrixB_T[j * sizeB.Item1 + i] = matrixB[i * sizeB.Item2 + j];

			var resultSize = resultRowCount * resultColumnCount;
			var result = new float[resultSize];
			for (int i = 0; i < resultSize; i++)
				result[i] = 0.0f;

			for (int rY = 0; rY < resultRowCount; rY++)
				for (int rX = 0; rX < resultColumnCount; rX++)
					for (int sD = 0; sD < sharedDimension; sD++)
						result[rY * resultColumnCount + rX] += matrixA[rY * sharedDimension + sD] * matrixB_T[rX * sharedDimension + sD];

			return result;
		}

		public static float[] MatrixMult_MonoArray_Strassen(float[] matrixA, float[] matrixB, Tuple<int, int> sizeA, Tuple<int, int> sizeB)
		{
			throw new NotImplementedException();
		}

	}
}
