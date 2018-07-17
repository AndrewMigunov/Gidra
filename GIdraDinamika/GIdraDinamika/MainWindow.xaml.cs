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
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.IO;

namespace GIdraDinamika
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }
        private void OK_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                int N = Convert.ToInt32(NTextbox.Text), M = Convert.ToInt32(MTextbox.Text);
                double Re = Convert.ToDouble(ReTextbox.Text), L = Convert.ToDouble(LTextbox.Text);
                double H = Convert.ToDouble(HTextbox.Text), l = Convert.ToDouble(lTextbox.Text);
                double epsPhi = Convert.ToDouble(epsPhiTextbox.Text), epsOmega = Convert.ToDouble(epsOmegaTextbox.Text);
                double U0 = Convert.ToDouble(U0Textbox.Text);
                double[,] u, v,phi,omega;
                u = new double[N + 2, M + 2];
                v = new double[N + 2, M + 2];
                phi= new double[N + 2, M + 2];
                omega= new double[N + 2, M + 2];
                double dx = L / (N +1);
                double dy = H / (M +1);
                Gidra gidra = new Gidra(N, M, Re, epsOmega, epsPhi, U0, dx, dy, l);
                gidra.culculate();
                u = gidra.get_u();
                v = gidra.get_v();
                phi = gidra.get_phi();
                omega = gidra.get_omega();
                Graph g = new Graph(N, M, u, v, dx, dy,U0);
                StreamWriter phiFile, omegaFile;
                phiFile = new StreamWriter("C:\\с#/GIdraDinamika/GIdraDinamika/phi.txt");
                omegaFile = new StreamWriter("C:\\с#/GIdraDinamika/GIdraDinamika/omega.txt");
                
                //for (int j = 0; j < M + 2; j++)
                //{
                //    for (int i = 0; i < N + 2; i++)
                //    {
                //        phiFile.Write($"{Math.Round(phi[i, j],2).ToString()} ");
                //        omegaFile.Write($"{Math.Round(omega[i, j],2).ToString()} ");
                //    }
                //    phiFile.Write("\n");
                //    omegaFile.Write("\n");
                //}
                
                for(int i = 0; i < N + 2; i++)
                {
                    for(int j = 0; j < M + 2; j++)
                    {
                        string phiLine = $"{Math.Round(i * dx, 2)} {Math.Round(j * dy, 2)} {Math.Round(phi[i, j], 2).ToString()} \n";
                        for(int k = 0; k < phiLine.Length; k++)
                        {
                            if (phiLine[k] == ',')
                            {
                                phiLine=phiLine.Remove(k, 1).Insert(k, ".");
                            }
                        }
                        string omegaLine = $"{Math.Round(i * dx, 2)} {Math.Round(j * dy, 2)} {Math.Round(omega[i, j], 2).ToString()} \n";
                        for (int k = 0; k < omegaLine.Length; k++)
                        {
                            if (omegaLine[k] == ',')
                            {                                
                                omegaLine=omegaLine.Remove(k, 1).Insert(k, ".");
                            }
                        }
                        phiFile.Write(phiLine);
                        omegaFile.Write(omegaLine);
                    }
                }
                
                phiFile.Close();
                omegaFile.Close();
                g.ShowDialog();

            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
            }

}
    }
}
