using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace lab01
{
    public partial class Form1 : Form
    {
        // Матрица прогонки для Х
        C3Matrix TridiagMatrX;
        // Матрица прогонки для Psi
        C3Matrix TridiagMatrPsi;

        // Матрица для решения первой задачи методом прогонки
        CMatrix mMatrX;
        // Матрица для решения второй (вспомогательной) задачи методом прогонки
        CMatrix mMatrPsi;

        // Количество узлов сетки по длине стержня
        int mlGridSize = 100;
        // Количество узлов сетки по времени
        int mtGridSize = 100;
        
        double a2 = 1; // Постоянная a^2
        double l = 1;  // Длина стержня, м
        double n = 5;  // Конcтанта ню
        double T = 20; // Время наблюдения, сек
        double EPS = 0.001;
        double R = 7; // 

        // Параметры для выбора шага спуска по методу 3 (с использованием константы Липшица)
        double eps0 = 0.1;
        double eps1 = 0.1;
        double alphaLip;

        // Функция p(t) - температура внешней среды, град. Цельсия
        const int p = 1;

        // Шаги сетки: h - вдоль стержня, tau - по времени
        double h, tau;

        // матрица управления по f
        CMatrix F;
        public double Func_f(double s, double t)
        {
            return Math.Pow(s, 2) + t;
            //return 2d;
        }

        // новое управление, полученное на итерации
        CMatrix F_1;

        // Массив значений исходной функции y(s)
        double[] X;
        double[] Y;

        CPlot pl;

        // Функция phi(s)
        public double Func_phi(double s)
        {
            return 0.5 * s + 1;
        }

        // Функция y(s)
        public double Func_y(double s)
        {
            //return Math.Cos(s);
            return 0.5d * Math.Exp(s*s/2d);
            // return 0.5 * Math.Cos(s*7d);
            //return s * s;
        }

        // Подинтегральная функция
        public double Func_I(double f_k, double psi_k, double alpha_k)
        {
            return Math.Pow((f_k - alpha_k * psi_k), 2);
        }

        // Функция численного интегрирования методом трапеций (двойной интеграл)
        public double DoubleIntegral(CMatrix a_F, CMatrix a_matrPsi, double a_alpha)
        {
            // Значение двойного интеграла
            double S = 0;

            // Массив значений одномерных интегралов по переменной s для каждого временного слоя
            double[] S_str = new double[a_F.numRows];

            for(int i = 0; i < a_F.numRows; i++)
            {
                for (int j = 1; j < a_F.numCols - 1; j++)
                {
                    S_str[i] += Func_I(a_F[i,j], a_matrPsi[i,j], a_alpha);
                }
                S_str[i] += (Func_I(a_F[i,0], a_matrPsi[i,0], a_alpha) + 
                             Func_I(a_F[i, a_F.numCols - 1], a_matrPsi[i, a_F.numCols - 1], a_alpha)
                            ) / 2d;
                S_str[i] *= h;
            }

            for(int k = 1; k < a_F.numRows - 1; k++)
            {
                S += S_str[k];
            }
            S += (S_str[0] + S_str[a_F.numRows - 1])/2d;
            S *= tau;

            return S;
        }

        // Функция численного интегрирования методом трапеций (двойной интеграл)
        public double NormF(CMatrix a_F1, CMatrix a_F2)
        {
            // Значение двойного интеграла
            double S = 0;

            // Массив значений одномерных интегралов по переменной s для каждого временного слоя
            double[] S_str = new double[a_F1.numRows];

            for (int i = 0; i < a_F1.numRows; i++)
            {
                for (int j = 1; j < a_F1.numCols - 1; j++)
                {
                    S_str[i] += Math.Pow(a_F1[i, j] - a_F2[i,j],2);
                }
                S_str[i] += (Math.Pow(a_F1[i, 0] - a_F2[i, 0],2) +
                             Math.Pow(a_F1[i, a_F1.numCols-1] - a_F2[i, a_F1.numCols - 1],2)
                            ) / 2d;
                S_str[i] *= h;
            }

            for (int k = 1; k < a_F1.numRows - 1; k++)
            {
                S += S_str[k];
            }
            S += (S_str[0] + S_str[a_F1.numRows - 1]) / 2d;
            S *= tau;

            return Math.Sqrt(S);
        }

        // Функция метода прогонки для получения решения на одном временном слое
        public void TridiagonalMatrixAlgorithm(int a_n, double []a_x, C3Matrix a_matr, double []a_b)
        {
            double m;
            double[] gamma = new double[a_n];
            for(int i = 0; i < a_n; i++)
            {
                gamma[i] = a_matr.mD[i];
            }

            for (int i = 1; i < a_n; i++)
            {
                m = a_matr.unD[i] / gamma[i-1];
                gamma[i] = gamma[i] - m * a_matr.upD[i - 1];
                a_b[i] = a_b[i] - m * a_b[i - 1];
            }

            a_x[a_n - 1] = a_b[a_n - 1] / gamma[a_n - 1];

            for (int i = a_n - 2; i >= 0; i--)
            {
                a_x[i] = (a_b[i] - a_matr.upD[i] * a_x[i + 1]) / gamma[i];
            }

            return;
        }

        public double Functional(CMatrix a_Matr)
        {
            double s = 0;

            for(int i = 0; i < a_Matr.numCols; i++)
            {
                s += Math.Pow(a_Matr[a_Matr.numRows-1,i] - Func_y(i*h), 2);
            }

            return s;
        }

        // Выбор шага спуска способ 3 (с использованием константы Липшица)
        public double getAlphaLip()
        {
            // Получение c0
            double val1 = (a2 * a2 * n * n + 2d * l) / (a2 * n);
            double val2 = (2d * l * l) / (a2);
            double c0 = Math.Max(val1, val2);

            // получение c1
            val1 = (a2 * n) / eps0;
            val2 = 1d / (a2 * eps0);
            double c1 = Math.Max(val1, val2);

            // Вычисление константы Липшица
            double lip = Math.Sqrt(2d * c0 * c1);

            // Вычисление alpha
            Random rnd = new Random();
            double val = rnd.NextDouble();

            return val * eps0 + (1d - val) * eps1;
        }

        // Вычисление проекции точки на множество
        public CMatrix CalcProjection(CMatrix a_F)
        {
            CMatrix Pu = new CMatrix(mtGridSize, mlGridSize);
            CMatrix f_1  = new CMatrix(mtGridSize, mlGridSize);

            double norm = NormF(a_F, f_1);

            if(norm > R)
            {
                for(int i = 0; i < Pu.numRows; i++)
                {
                    for(int j = 0; j < Pu.numCols; j++)
                    {
                        Pu[i,j] = (R / norm) * a_F[i,j];
                    }
                }
                return Pu;
            }
            else
            {
                return a_F;
            }
        }

        // Вычисление шага спуска через проекцию
        public double GetAlphaByProjection(CMatrix a_X, CMatrix a_F)
        {
            // Матрица для новой функции x(s,t)
            CMatrix x_new = new CMatrix(mtGridSize, mlGridSize);
            // Матрица для хранения нового управления f(s,t)
            CMatrix f_new = new CMatrix(mtGridSize, mlGridSize);

            // Начальное значение шага спуска
            double alpha = 5;

            double I;            // Двойной интеграл для вычисления управления
            CMatrix Projection;  // Проекция f на множество U
            double norm;         // Норма ||U_k - U_k1||
            double Ju, Ju1;      // Функционалы от Uk и Uk_1

            do
            {
                I = DoubleIntegral(a_F, mMatrPsi, alpha);

                if (I <= R * R)
                {
                    for (int i = 0; i < f_new.numRows; i++)
                    {
                        for (int j = 0; j < f_new.numCols; j++)
                        {
                            f_new[i,j] = a_F[i, j] - alpha * mMatrPsi[i, j];
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < F_1.numRows; i++)
                    {
                        for (int j = 0; j < F_1.numCols; j++)
                        {
                            f_new[i, j] = (R * (a_F[i, j] - alpha * mMatrPsi[i, j])) / Math.Sqrt(I);
                        }
                    }
                }

                Projection = CalcProjection(f_new);
                norm = NormF(a_F, f_new);

                calcMatrX(x_new, f_new, TridiagMatrX);
                Ju = Functional(a_X);
                Ju1 = Functional(x_new);

                alpha /= 2.0;

            } while ((Ju - Ju1) < EPS * norm * norm);

            return alpha * 2.0;
        }

        // Решение краевой задачи для получения матрицы иксов
        public void calcMatrX(CMatrix a_matrX, CMatrix a_F, C3Matrix a_3M)
        {
            for (int i = 0; i < a_matrX.numCols; i++)
            {
                a_matrX[0, i] = Func_phi(i * h);
            }

            double[] b = new double[a_matrX.numCols];
            double[] x = new double[a_matrX.numCols];

            for (int j = 1; j < a_matrX.numRows; j++)
            {
                for (int k = 0; k < a_matrX.numCols; k++)
                {
                    if(k == 0)
                    {
                        b[k] = (a_F[j, k] + a_matrX[j-1, k] / tau) * (h * h) / 2 * a2;
                    }
                    else if (k == a_matrX.numCols - 1)
                    {
                        b[k] = (a_F[j, k] + a_matrX[j - 1, k] / tau) * (h * h) / (2d * a2) + p * n * h;
                    }
                    else if (k > 0 && k < a_matrX.numCols - 1)
                    {
                        b[k] = a_matrX[j - 1, k] + a_F[j, k] * tau;
                    }
                }

                TridiagonalMatrixAlgorithm(a_matrX.numCols, x, a_3M, b);

                for(int q = 0; q < a_matrX.numCols; q++)
                {
                    a_matrX[j, q] = x[q];
                }
            }
        }

        // Решение краевой задачи для получения матрицы кси
        public void calcMatrPsi(CMatrix a_matrPsi, C3Matrix a_3M)
        {
            for (int i = 0; i < a_matrPsi.numCols; i++)
            {
                a_matrPsi[mtGridSize-1, i] = 2*(mMatrX[mtGridSize - 1,i] - Func_y(i*h));
            }

            double[] b = new double[a_matrPsi.numCols];
            double[] x = new double[a_matrPsi.numCols];

            for (int j = a_matrPsi.numRows-2; j >= 0; j--)
            {
                for(int k = 0; k < a_matrPsi.numCols; k++)
                {
                    if(k == 0)
                    {
                        b[k] = a_matrPsi[j + 1, k] * (h * h / (2 * a2 * tau));
                    }
                    else if (k == a_matrPsi.numCols - 1)
                    {
                        b[k] = a_matrPsi[j + 1, k] * (h * h) / (2 * a2 * tau);
                    }
                    else if (k > 0 && k < a_matrPsi.numCols - 1)
                    {
                        b[k] = a_matrPsi[j + 1, k];
                    }
                }
                
                TridiagonalMatrixAlgorithm(a_matrPsi.numCols, x, a_3M, b);

                for (int q = 0; q < a_matrPsi.numCols; q++)
                {
                    a_matrPsi[j, q] = x[q];
                }
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            // Трёхдиагональные матрицы
            mMatrX = new CMatrix(mtGridSize, mlGridSize);
            mMatrPsi = new CMatrix(mtGridSize, mlGridSize);

            // Вычисление шагов сетки
            h = l / mlGridSize;
            tau = T / mtGridSize;

            // Задание матрицы начального управления
            F = new CMatrix(mtGridSize, mlGridSize);
            for (int i = 0; i < F.numRows; i++)
            {
                for(int j = 0; j < F.numCols; j++)
                {
                    F[i, j] = Func_f(i*tau,j*h);
                }
            }

            // новое управление, полученное на итерации
            F_1 = new CMatrix(mtGridSize, mlGridSize);

            // Заполнения TridiagMatrX
            TridiagMatrX = new C3Matrix(mlGridSize);
            TridiagMatrX.mD[0] = 1 + (h*h) / (a2*2*tau);
            TridiagMatrX.upD[0] = -1;

            TridiagMatrX.mD[mlGridSize - 1] = 1 + (h * h) / (a2 * 2 * tau) + h*n;
            TridiagMatrX.unD[mlGridSize - 1] = -1;

            for (int i = 1; i < mlGridSize - 1; i++)
            {
                TridiagMatrX.unD[i] = -a2 * tau / (h * h);
                TridiagMatrX.mD[i] = 1 + 2 * a2 * tau / (h * h);
                TridiagMatrX.upD[i] = -a2 * tau / (h * h);
            }

            // Заполнения TridiagMatrPsi
            TridiagMatrPsi = new C3Matrix(mlGridSize);
            TridiagMatrPsi.mD[0] = -(1 - (h * h) / (a2 * 2 * tau));
            TridiagMatrPsi.upD[0] = 1;

            TridiagMatrPsi.mD[mlGridSize - 1] = -1 - n * h + (h * h) / (2 * a2 * tau);
            TridiagMatrPsi.unD[mlGridSize - 1] = 1;

            for (int i = 1; i < mlGridSize - 1; i++)
            {
                TridiagMatrPsi.unD[i] = -a2 * tau / (h*h);
                TridiagMatrPsi.mD[i] = 1 + 2 * a2 * tau / (h*h);
                TridiagMatrPsi.upD[i] = -a2 * tau / (h*h);
            }

            // Шаг спуска по Липшицу
            alphaLip = getAlphaLip();

            // Инициализация графика
            pl = new CPlot(chart1);
            pl.Init();

            // Вывод исходной функции
            X = new double[mlGridSize];
            Y = new double[mlGridSize];
            for (int i = 0; i < mlGridSize; i++)
            {
                X[i] = i * h;
                Y[i] = Func_y(X[i]);
            }
            pl.Show(X, Y, 0);
        }

        private void gradientProjectionMethodToolStripMenuItem_Click(object sender, EventArgs e)
        {
            lbl_numIter.Visible = true;

            double numIter = 0;   // Количество итераций
            double I;             // Величина двойного интеграла
            double norm;          //   || u_k - u_k1 ||

            do
            {
                if(numIter != 0)
                {
                    for (int i = 0; i < F_1.numRows; i++)
                    {
                        for (int j = 0; j < F_1.numCols; j++)
                        {
                            F[i, j] = F_1[i, j];
                        }
                    }
                }

                // Решение первой краевой задачи
                calcMatrX(mMatrX, F, TridiagMatrX);
                // Решение второй краевой задачи
                calcMatrPsi(mMatrPsi, TridiagMatrPsi);

                I = DoubleIntegral(F, mMatrPsi, alphaLip);

                if (I <= R * R)
                {
                    for (int i = 0; i < F_1.numRows; i++)
                    {
                        for (int j = 0; j < F_1.numCols; j++)
                        {
                            F_1[i, j] = F[i, j] - alphaLip * mMatrPsi[i, j];
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < F_1.numRows; i++)
                    {
                        for (int j = 0; j < F_1.numCols; j++)
                        {
                            F_1[i, j] = (R * (F[i, j] - alphaLip * mMatrPsi[i, j])) / Math.Sqrt(I);
                        }
                    }
                }

                numIter++;

                norm = NormF(F, F_1);

                pl.Show(X, Y, 0);
                pl.Show(mMatrX, X, 1);
                
                lbl_numIter.Text = "Number of iterations: " + numIter;
                this.Refresh();

            } while (norm >= EPS);

            // Запись матрицы управления в файл
            string filePath = ".\\dataLip.txt";
            using (StreamWriter sw = new StreamWriter(filePath, true, System.Text.Encoding.Default))
            {
                for (int i = 0; i < F_1.numCols; i++)
                {
                    for (int j = 0; j < F_1.numRows; j++)
                    {
                        sw.Write(Math.Round(F[i, j], 5));
                        sw.Write(" ");
                    }
                    sw.WriteLine("");
                }
            }
        }

        private void gradientProjectionMethodProjectionToolStripMenuItem_Click(object sender, EventArgs e)
        {
            lbl_numIter.Visible = true;

            double numIter = 0;   // Количество итераций
            double I;             // Величина двойного интеграла
            double norm;          //   || u_k - u_k1 ||
            double alpha;

            do
            {
                // Начало итерации метода проекции градиента
                if (numIter != 0)
                {
                    for (int i = 0; i < F_1.numRows; i++)
                    {
                        for (int j = 0; j < F_1.numCols; j++)
                        {
                            F[i, j] = F_1[i, j];
                        }
                    }
                }

                // Решение первой краевой задачи
                calcMatrX(mMatrX, F, TridiagMatrX);
                // Решение второй краевой задачи
                calcMatrPsi(mMatrPsi, TridiagMatrPsi);

                // Вычисление шага спуска через проекцию
                alpha = GetAlphaByProjection(mMatrX, F);

                I = DoubleIntegral(F, mMatrPsi, alpha);

                if (I <= R * R)
                {
                    for (int i = 0; i < F_1.numRows; i++)
                    {
                        for (int j = 0; j < F_1.numCols; j++)
                        {
                            F_1[i, j] = F[i, j] - alpha * mMatrPsi[i, j];
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < F_1.numRows; i++)
                    {
                        for (int j = 0; j < F_1.numCols; j++)
                        {
                            F_1[i, j] = (R * (F[i, j] - alpha * mMatrPsi[i, j])) / Math.Sqrt(I);
                        }
                    }
                }

                numIter++;

                norm = NormF(F, F_1);

                pl.Show(X, Y, 0);
                pl.Show(mMatrX, X, 1);

                lbl_numIter.Text = "Number of iterations: " + numIter;
                this.Refresh();

            } while (norm >= EPS);

            // Запись матрицы управления в файл
            string filePath = ".\\dataProj.txt";
            using (StreamWriter sw = new StreamWriter(filePath, true, System.Text.Encoding.Default))
            {
                for (int i = 0; i < F_1.numCols; i++)
                {
                    for (int j = 0; j < F_1.numRows; j++)
                    {
                        sw.Write(Math.Round(F[i, j], 5));
                        sw.Write(" ");
                    }
                    sw.WriteLine("");
                }
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            // Заполнение матрицы начального управления
            for (int i = 0; i < F.numRows; i++)
            {
                for (int j = 0; j < F.numCols; j++)
                {
                    F[i, j] = Func_f(i * tau, j * h);
                }
            }

            // Обнуление нового управление, полученное на итерации
            F_1 = new CMatrix(mtGridSize, mlGridSize);

            // Инициализация графика
            pl = new CPlot(chart1);
            pl.Init();
            pl.Show(X, Y, 0);

            lbl_numIter.Text = "Number of iterations: ";
        }

        private void exitToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }

        public Form1()
        {
            InitializeComponent();
            this.Load += Form1_Load;  // подписка на событие Form1_Load
        }
    }
}
