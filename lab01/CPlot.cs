using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.DataVisualization.Charting;

namespace lab01
{
    public class CPlot
    {
        Chart c;

        public CPlot(Chart chart)
        {
            c = chart;
        }

        public void Init()
        {
            c.Series.Clear();
            c.Series.Add(new Series("Source function"));
            c.Series.Add(new Series("Function calculated"));

            for (int i = 0; i < 2; i++)
            {
                c.Series[i].ChartType = SeriesChartType.Line;
                c.Series[i].BorderWidth = 4;
            }
        }

        public void Show(CMatrix a_M, double[] x, int SeriesNumer)
        {
            c.Series[SeriesNumer].Points.Clear();
            for (int i = 0; i < a_M.numCols; i++)
                c.Series[SeriesNumer].Points.AddXY(x[i], a_M[a_M.numRows - 1, i]);
        }

        public void Show(double[] x, double []y, int SeriesNumer)
        {
            c.Series[SeriesNumer].Points.Clear();
            for (int i = 0; i < y.Length; i++)
                c.Series[SeriesNumer].Points.AddXY(x[i], y[i]);
        }

    }
}
