
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
			throw new NotImplementedException();
		}
	}
}
