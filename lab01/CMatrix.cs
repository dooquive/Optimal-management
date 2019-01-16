using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab01
{
    public class CMatrix
    {
        private int n; // количество строк
        private int m; // количество столбцов
        private double[,] M;

        public CMatrix(int a_n, int a_m)
        {
            n = a_n;
            m = a_m;
            M = new double[a_n,a_m];
            for (int i = 0; i < a_n; i++)
            {
                for (int j=0; j<a_m; j++)
                {
                    M[i,j] = 0;
                }
            }
        }

        public CMatrix()
        {
            n = 1;
            m = 1;
            M = new double[1,1];
            M[0,0] = 0;
        }

        public int numRows
        {
            get { return n; }
            set { n = value; }
        }

        public int numCols
        {
            get { return m; }
            set { m = value; }
        }

        public double this[int i, int j]
        {
            //Метод доступа для чтения
            get
            {
                return M[i, j];
            }
            //Метод доступа для установки
            set
            {
                M[i, j] = value;
            }
        }

        public static CMatrix operator +(CMatrix A, CMatrix B)
        {
            if (A.n == B.n && A.m == B.m)
            {
                CMatrix C = new CMatrix(A.n, A.m);

                for (int i = 0; i < A.n; i++)
                {
                    for (int j = 0; j < A.m; j++)
                    {
                        C[i, j] = A[i, j] + B[i, j];
                    }
                }
                return C;
            }
            return new CMatrix();
        }

        public static CMatrix operator -(CMatrix A, CMatrix B)
        {
            if (A.n == B.n && A.m == B.m)
            {
                CMatrix C = new CMatrix(A.n, A.m);

                for (int i = 0; i < A.n; i++)
                {
                    for (int j = 0; j < A.m; j++)
                    {
                        C[i, j] = A[i, j] - B[i, j];
                    }
                }
                return C;
            }
            return new CMatrix();
        }

        public static CMatrix operator *(double a_val, CMatrix A)
        {
            for (int i = 0; i < A.n; i++)
            {
                for(int j = 0; j < A.m; j++)
                {
                    A[i, j] *= a_val;
                }
            }
            return A;
        }

        public double Norm()
        {
            double S = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    S += Math.Pow(M[i, j], 2);
                }
            }
            return Math.Pow(S, 0.5);
        }

    }
}
