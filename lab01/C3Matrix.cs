using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab01
{
    public class C3Matrix
    {
        public C3Matrix(int aSize)
        {
            mSize = aSize;

            mDiag = new double[aSize];
            unDiag = new double[aSize];
            upDiag = new double[aSize];
            for(int i=0; i < aSize; i++)
            {
                Array.Clear(mDiag, 0, mSize);
                Array.Clear(unDiag, 0, mSize);
                Array.Clear(upDiag, 0, mSize);
            }
        }

        public double[] upD
        {
            get
            {
                return upDiag;
            }

        }

        public double[] unD
        {
            get
            {
                return unDiag;
            }

        }

        public double[] mD
        {
            get
            {
                return mDiag;
            }

        }

        public void Clear()
        {
            Array.Clear(mDiag, 0, mSize);
            Array.Clear(unDiag, 0, mSize);
            Array.Clear(upDiag, 0, mSize);
        }

        // Размерность матрицы
        private int mSize;
        private double [] mDiag;  // массив главной диагонали
        private double [] unDiag; // массив поддиагонали
        private double [] upDiag; // массив наддиагонали
    }
}
