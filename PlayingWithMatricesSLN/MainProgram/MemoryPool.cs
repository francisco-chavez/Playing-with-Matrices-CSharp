
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace MainProgram
{
	public static class MemoryPool
	{

		private static Dictionary<int, List<float[]>> _freeArrays = new Dictionary<int, List<float[]>>();


		public static float[] GetArray(int blockSize)
		{
			if (!_freeArrays.ContainsKey(blockSize))
			{
				_freeArrays.Add(blockSize, new List<float[]>(21));
				return new float[blockSize];
			}

			if (_freeArrays[blockSize].Count == 0)
				return new float[blockSize];

			var arrayCollection = _freeArrays[blockSize];
			var array = arrayCollection[arrayCollection.Count - 1];
			arrayCollection.RemoveAt(arrayCollection.Count - 1);

			return array;
		}

		public static float[] GetZerosArray(int blockSize)
		{
			var array = GetArray(blockSize);

			for (int i = 0; i < blockSize; i++)
				array[i] = 0.0f;

			return array;
		}

		public static float[] GetRandomArray(int blockSize, int randomSeed, float min = 0.0f, float max = 1.0f)
		{
			var array = GetArray(blockSize);

			Random randomGen = new Random(randomSeed);
			var range = ((double) max) - ((double) min);

			for (int i = 0; i < blockSize; i++)
				array[i] = (float) (range * randomGen.NextDouble() + min);

			return array;
		}

		public static void ReturnArray(int blockSize, float[] array)
		{
			//if (!_freeArrays.ContainsKey(blockSize))
			//	_freeArrays.Add(blockSize, new List<float[]>(21));

			_freeArrays[blockSize].Add(array);
		}

		public static void ClearMemoryPool()
		{
			foreach (var blockSize in new List<int>(_freeArrays.Keys))
			{
				_freeArrays[blockSize].Clear();
				_freeArrays[blockSize] = null;
			}

			_freeArrays.Clear();
		}

	}
}
