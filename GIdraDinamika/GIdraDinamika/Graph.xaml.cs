using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;

namespace GIdraDinamika
{
    /// <summary>
    /// Логика взаимодействия для Graph.xaml
    /// </summary>
    public partial class Graph : Window
    {
        double[,] u;
        double[,] v;
        int N, M;
        double dx, dy,u0;

        public Graph(int N, int M, double[,] u, double[,] v, double dx, double dy,double u0)
        {
            this.N = N;
            this.M = M;
            this.u = u;
            this.v = v;
            this.dx = dx;
            this.dy = dy;
            this.u0 = u0;
            InitializeComponent();
            this.plotView.Model = ZeroCrossing();
        }
        public PlotModel ZeroCrossing()
        {
            var plotModel = new PlotModel();
            plotModel.PlotAreaBorderThickness = new OxyThickness(0.0);
            plotModel.PlotMargins = new OxyThickness(10);

            var LinearAxis = new LinearAxis();
            LinearAxis.Maximum = (M + 1) * dy;
            LinearAxis.Minimum = 0;
            LinearAxis.PositionAtZeroCrossing = true;
            LinearAxis.TickStyle = TickStyle.Crossing;
            LinearAxis.MajorGridlineStyle = LineStyle.Solid;
            plotModel.Axes.Add(LinearAxis);
            var secondLinearAxis = new LinearAxis();
            secondLinearAxis.Maximum = (N + 1) * dx;
            secondLinearAxis.Minimum = 0;
            secondLinearAxis.PositionAtZeroCrossing = true;
            secondLinearAxis.TickStyle = TickStyle.Crossing;
            secondLinearAxis.MajorGridlineStyle = LineStyle.Solid;
            secondLinearAxis.Position = AxisPosition.Bottom;
            plotModel.Axes.Add(secondLinearAxis);
            // plotModel.Series.Add(new FunctionSeries(fun, 0, 3, 0.01, nameof(function)));

            //u = new double[N + 2, M + 2];
            //v = new double[N + 2, M + 2];
            //for(int i = 0; i < N + 2; i++)
            //{
            //    for(int j = 0; j < M + 2; j++)
            //    {
            //        u[i, j] = 0;
            //        v[i, j] = 0;
            //    }
            //}
            //for (int i = 1; i < N + 1; i++)
            //{
            //    for (int j = 1; j < M + 1; j++)
            //    {
            //        u[i, j] = 1;
            //        v[i, j] = 1;
            //    }
            //}
            //Gidra gid = new Gidra();
            //gid.culculate();
            //u = gid.get_u();
            //v = gid.get_v();
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < M + 2; j++)
                {
                    double d;
                    if (dx >= dy)
                    {
                        d = dy;
                    }
                    else
                    {
                        d = dx;
                    }
                    //d = d * 2 / 3;
                    double modul = Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
                    if (modul >= u0)
                    {
                        d =d* 3 / 4;
                    }
                    else
                    {
                        if (modul > u0 / 2)
                        {
                            d =d* 1 / 2;
                        }
                        else
                        {
                            d =d* 1 / 4;
                        }
                    }
                    double x_m = (i * dx + 5 * (i * dx + u[i, j] * d / modul)) / (1 + 5);
                    double y_m = (j * dy + 5 * (j * dy + v[i, j] * d / modul)) / (1 + 5);
                    double a = modul / u[i, j] / d, b = -modul / v[i, j] / d, c = -i * dx * modul / u[i, j] / d + j * dy * modul / v[i, j] / d;

                    var line = new LineSeries();
                    line.Color = OxyColors.Black;
                    line.Points.Add(new DataPoint(i * dx, j * dy));
                    line.Points.Add(new DataPoint(i * dx + u[i, j] * d / modul, j * dy + v[i, j] * d / modul));
                    line.Points.Add(new DataPoint((-c - 0.3 + b * b * x_m / a - b * y_m) / (b * b / a + a),
                        (-c - 0.3 - a * x_m + a * a * y_m / b) / (b + a * a / b)));
                    line.Points.Add(new DataPoint(i * dx + d * u[i, j] / modul, j * dy + d * v[i, j] / modul));
                    line.Points.Add(new DataPoint((-c + 0.3 + b * b * x_m / a - b * y_m) / (b * b / a + a),
                        (-c + 0.3 - a * x_m + a * a * y_m / b) / (b + a * a / b)));

                    plotModel.Series.Add(line);

                }
            }

            //var series1 = new LineSeries { Title = $"Series {m}" };
            //for (int i = 0; i <= n; i++)
            //{
            //    z = a + i * h;
            //    series1.Points.Add(new DataPoint(z, y[i, m]));
            //}
            //plotModel.Series.Add(series1);

            return plotModel;
        }
    }
}
