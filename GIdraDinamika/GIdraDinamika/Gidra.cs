using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GIdraDinamika
{
    class Gidra
    {
        int N, M,MaxIterPhi=100,MaxIterOmega=10000,MinIterPhi=30,MinIterOmega=30;
        double[,] u, v, phi, omega, NewPhi, NewOmega, NewPhi_t;
        double Re, u0, dt = 0.00001, epsPhi, epsOmega, T = 1;
        double dx, dy;
        double l;
        bool[,] triangle;

        public Gidra(int N, int M, double Re, double epsOmega, double epsPhi, double u0, double dx, double dy, double l)
        {
            this.N = N;
            this.M = M;
            this.Re = Re;
            this.epsOmega = epsOmega;
            this.epsPhi = epsPhi;
            this.u0 = u0;
            this.dx = dx;
            this.dy = dy;
            this.l = l;

            phi = new double[N + 2, M + 2];
            omega = new double[N + 2, M + 2];
            NewPhi = new double[N + 2, M + 2];
            NewPhi_t = new double[N + 2, M + 2];
            NewOmega = new double[N + 2, M + 2];
           
            triangle = new bool[N + 2, M + 2];
            u = new double[N + 2, M + 2];
            v = new double[N + 2, M + 2];
            double H = dy * (M + 1);
            int delta = (int)(l / dx);
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < M + 2; j++)
                {
                    phi[i, j] = 0;
                    omega[i, j] = 0;
                    u[i, j] = 0;
                    v[i, j] = 0;
                    NewPhi[i, j] = 0;
                    NewOmega[i, j] = 0;
                    NewPhi_t[i, j] = 0;
                }
            }
            //индексы треугольника
            for (int i = 0; i <= N + 1; i++)
            {
                for (int j = 0; j <= M + 1; j++)
                {
                    triangle[i, j] = false;
                }

            }

            int k = 0;
            for (int i = N / 2; i <= N / 2 + delta; i++)
            {
                for (int j = (M+1) / 2 - delta + k; j <= (M+1) / 2 + delta - k; j++)
                {
                    triangle[i, j] = true;
                }
                k++;
            }

            //начальные условия на границах

            //на левой границе
            for (int j = 0; j < M + 2; j++)
            {
                double y = j * dy;
                u[0, j] = u0 * y * (H - y);
                phi[0, j] = u0 * (y * y * H / 2 - y * y * y / 3);
                NewPhi[0, j] = u0 * (y * y * H / 2 - y * y * y / 3);
            }
            
            //на верхней границе
            for (int i = 0; i < N + 1; i++)
            {
                phi[i, M + 1] = u0 * H * H * H / 6;
                NewPhi[i, M + 1] = u0 * H * H * H / 6;
            }
            //на нижней границе нули

            // правой границе
            for (int j = 0; j < M + 2; j++)
            {
                phi[N + 1, j] = phi[N, j];
                NewPhi[N + 1, j] = phi[N, j];
            }
            //на треугольнике
            for(int i = 0; i < N + 2; i++)
            {
                for(int j = 0; j < M + 2; j++)
                {
                    //боковая сторона
                    if(triangle[i,j]&&!triangle[i-1, j])
                    {
                        phi[i, j] = u0 * H * H * H / 12;
                        NewPhi[i, j] = u0 * H * H * H / 12;
                        continue;
                    }
                    //верхняя сторона
                    if (triangle[i, j] && !triangle[i, j + 1])
                    {
                        phi[i, j] = u0 * H * H * H / 12;
                        NewPhi[i, j] = u0 * H * H * H / 12;
                        continue;
                    }
                    //нижняя сторона
                    if (triangle[i, j] && !triangle[i, j - 1])
                    {
                        phi[i, j] = u0 * H * H * H / 12;
                        NewPhi[i, j] = u0 * H * H * H / 12;
                        continue;
                    }
                }

            }
            //for (int j = -delta + M / 2; j <= M / 2 + delta; j++)
            //{
            //    phi[N / 2, j] = u0 * H * H * H / 12;
            //    NewPhi[N / 2, j] = u0 * H * H * H / 12;
            //}
            //k = 0;
            //for (int j = M / 2 + delta; j >= M / 2; j--)
            //{
            //    phi[N / 2 + k, j] = u0 * H * H * H / 12;
            //    NewPhi[N / 2 + k, j] = u0 * H * H * H / 12;
            //    k++;
            //}
            //k = 0;
            //for (int j = M / 2 - delta; j <= M / 2; j++)
            //{
            //    phi[N / 2 + k, j] = u0 * H * H * H / 12;
            //    NewPhi[N / 2 + k, j] = u0 * H * H * H / 12;
            //    k++;
            //}
            NewVortex();
            NewBorder();
            NewSpeeds();
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < M + 2; j++)
                {
                    phi[i, j] = NewPhi[i, j];
                    omega[i, j] = NewOmega[i, j];
                }
            }

        }
        
        private void NewVortex()
        {
            for (int i = 1; i < N + 1; i++)
            {
                for (int j = 1; j < M + 1; j++)
                {
                    if (triangle[i, j])
                    {
                        continue;
                    }
                    double c_x = dt * u[i, j] / dx;
                    double c_y = dt * v[i, j] / dy;
                    double d_x = dt / (Re * dx * dx);
                    double d_y = dt / (Re * dy * dy);
                    NewOmega[i, j] = omega[i, j] - c_x * (omega[i + 1, j] - omega[i - 1, j]) / 2 -
                        c_y * (omega[i, j + 1] - omega[i, j - 1]) / 2 +
                        d_x * (omega[i + 1, j] + omega[i - 1, j] - 2 * omega[i, j]) +
                        d_y * (omega[i, j + 1] + omega[i, j - 1] - 2 * omega[i, j]);
                }
            }
            //Parallel.For(1, N + 1, i => VortexParallel(i));
        }
        //private void VortexParallel(int i)
        //{
        //    for (int j = 1; j < M + 1; j++)
        //    {
        //        if (triangle[i, j])
        //        {
        //            continue;
        //        }
        //        double c_x = dt * u[i, j] / dx;
        //        double c_y = dt * v[i, j] / dy;
        //        double d_x = dt / (Re * dx * dx);
        //        double d_y = dt / (Re * dy * dy);
        //        NewOmega[i, j] = omega[i, j] - c_x * (omega[i + 1, j] - omega[i - 1, j]) / 2 -
        //            c_y * (omega[i, j + 1] - omega[i, j - 1]) / 2 +
        //            d_x * (omega[i + 1, j] + omega[i - 1, j] - 2 * omega[i, j]) +
        //            d_y * (omega[i, j + 1] + omega[i, j - 1] - 2 * omega[i, j]);
        //    }

        //}
        private void NewCurrent()
        {
            double beta = dx / dy;
            for (int i = 1; i < N + 1; i++)                
            {
                for (int j = 1; j < M + 1; j++)
                {
                    if (triangle[i, j])
                    {
                        continue;
                    }
                    NewPhi_t[i, j] = (NewPhi[i + 1, j] + NewPhi_t[i - 1, j] +
                        beta * beta * (NewPhi[i, j + 1] + NewPhi_t[i, j - 1]) -
                        dx * dx * NewOmega[i, j]) / (2 * (1 + beta * beta));
                    //NewPhi_t[i, j] = NewPhi[i, j] + 1.5 * (NewPhi[i + 1, j] + NewPhi_t[i - 1, j]
                    //    + beta * beta * NewPhi[i, j + 1] + beta * beta * NewPhi_t[i, j - 1]
                    //    - dx * dx * NewOmega[i, j] - 2 * (1 + beta * beta) * NewPhi[i, j])
                    //    / (2 * (1 + beta * beta));
                }
            }           

        }
        private void NewBorder()
        {
            //пересчет граничных условий

            //на левой границе
            NewOmega[0, 0] = NewOmega[1, 1] / (dx * dx + dy * dy);
            NewOmega[0, M+1] = NewOmega[1, M] / (dx * dx + dy * dy);
            for (int j = 1; j < M + 1; j++)
            {
                double H = dy * (M + 1);
                 //NewOmega[0, j] = (NewPhi[1, j] - 2 * NewPhi[0, j]) / (dx * dx) +
                  //   (NewPhi[0, j - 1] + NewPhi[0, j + 1] - 2 * NewPhi[0, j]) / (dy * dy);
                NewOmega[0, j] = u0 * H - u0 * 2 * j * dy + (NewPhi[0, j] + NewPhi[2, j] - 2 * NewPhi[1, j]) / (dx * dx);
            }
            //на верхней границе
            for(int i = 1; i < N + 1; i++)
            {
                //NewOmega[i, M + 1] = 3 * (NewPhi[i, M] - NewPhi[i, M + 1]) / (dy * dy)-NewOmega[i,M]/2;
                NewOmega[i, M + 1] = 2 * (NewPhi[i, M] - NewPhi[i, M + 1]) / (dy * dy);
            }
            //на нижней границе
            for (int i = 1; i < N + 1; i++)
            {
                NewOmega[i, 0] = 2 * (NewPhi[i, 1] - NewPhi[i, 0]) / (dy * dy);
            }
            //на правой границе
            for(int j = 0; j < M + 2; j++)
            {
                NewOmega[N + 1, j] = NewOmega[N, j];
            }
            //на треугольнике
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < M + 2; j++)
                {
                    //верхняя вершина
                    if (triangle[i, j] && !triangle[i - 1, j]
                        && !triangle[i + 1, j] && !triangle[i, j + 1])
                    {
                        NewOmega[i, j] = (NewPhi[i - 1, j] + NewPhi[i + 1, j] -
                            2 * NewPhi[i, j]) / (dx * dx) + (2*NewPhi[i, j + 1] -
                            2 * NewPhi[i, j]) / (dy * dy);
                        continue;
                    }
                    //нижняя вершина
                    if (triangle[i, j] && !triangle[i - 1, j]
                        && !triangle[i + 1, j] && !triangle[i, j - 1])
                    {
                        NewOmega[i, j] = (NewPhi[i - 1, j] +NewPhi[i + 1, j] -
                            2 * NewPhi[i, j]) / (dx * dx) + (2*NewPhi[i, j - 1] -
                            2 * NewPhi[i, j]) / (dy * dy);
                        continue;
                    }
                    //правая вершина
                    if (triangle[i, j] && !triangle[i, j - 1]
                        && !triangle[i + 1, j] && !triangle[i, j + 1])
                    {
                        //NewOmega[i, j] = (NewPhi[i + 1, j] -
                        //    2 * NewPhi[i, j]) / (dx * dx) + (NewPhi[i, j + 1] +
                        //    NewPhi[i, j - 1] -
                        //    2 * NewPhi[i, j]) / (dy * dy);
                        NewOmega[i, j] = (/*NewPhi[i - 1, j] +*/ 2 * NewPhi[i + 1, j] -
                            2 * NewPhi[i, j]) / (dx * dx) + ( NewPhi[i, j - 1] + NewPhi[i, j - 1] -
                            2 * NewPhi[i, j]) / (dy * dy);
                        continue;
                    }
                    //боковая сторона
                    if (triangle[i, j] && !triangle[i - 1, j])
                    {
                        NewOmega[i, j] = 2 * (NewPhi[i - 1, j] - NewPhi[i, j]) / (dx * dx);
                        continue;
                    }
                    //верхняя сторона
                    if (triangle[i, j] && !triangle[i, j + 1])
                    {
                        NewOmega[i, j] =  2*(NewPhi[i, j + 1] + NewPhi[i + 1, j]
                          - 2 * NewPhi[i, j]) / (dx * dx);
                        continue;
                    }
                    //нижняя сторона
                    if (triangle[i, j] && !triangle[i, j - 1])
                    {
                        NewOmega[i, j] =  2*(NewPhi[i, j - 1] + NewPhi[i + 1, j]
                          - 2 * NewPhi[i, j]) / (dx * dx);
                        continue;
                    }
                }
            }
                //int delta = (int)(l / dx);
                //NewOmega[N / 2, M / 2 + delta] = (NewPhi[N / 2 - 1, M / 2 + delta] + NewPhi[N / 2 + 1, M / 2 + delta] -
                //    2 * NewPhi[N / 2, M / 2 + delta]) / (dx * dx) + (2*NewPhi[N / 2, M / 2 + delta + 1] -
                //    2 * NewPhi[N / 2, M / 2 + delta]) / (dy * dy);
                //NewOmega[N / 2, M / 2 - delta] = (NewPhi[N / 2 - 1, M / 2 - delta] + NewPhi[N / 2 + 1, M / 2 - delta] -
                //    2 * NewPhi[N / 2, M / 2 - delta]) / (dx * dx) + (2*NewPhi[N / 2, M / 2 - delta - 1] -
                //    2 * NewPhi[N / 2, M / 2 - delta]) / (dy * dy);
                //NewOmega[N / 2 + delta, M / 2] = (2*NewPhi[N / 2 + delta + 1, M / 2] -
                //    2 * NewPhi[N / 2 + delta, M / 2]) / (dx * dx) + (NewPhi[N / 2 + delta, M / 2 + 1] -
                //    2 * NewPhi[N / 2 + delta, M / 2] + NewPhi[N / 2 + delta, M / 2 - 1]) / (dy * dy);
                //for (int j = -delta + M / 2 + 1; j <= M / 2 + delta - 1; j++)
                //{

                //    NewOmega[N / 2, j] = 2 * (NewPhi[N / 2 - 1, j] - NewPhi[N / 2, j]) / (dx * dx);
                //}
                //int k = 1;
                //for (int j = M / 2 + delta - 1; j > M / 2; j--)
                //{

                //    NewOmega[N / 2 + k, j] = 2 * (NewPhi[N / 2 + k, j + 1] + NewPhi[N / 2 + k + 1, j]
                //        - 2 * NewPhi[N / 2 + k, j]) / (dx * dx);
                //    //NewOmega[N / 2 + k, j] = NewPhi[N / 2 + k + 1, j + 1] / (dx * dx);
                //    k++;
                //}
                //k = 1;
                //for (int j = M / 2 - delta + 1; j < M / 2; j++)
                //{
                //    NewOmega[N / 2 + k, j] = 2 * (NewPhi[N / 2 + k, j - 1] + NewPhi[N / 2 + k + 1, j]
                //        - 2 * NewPhi[N / 2 + k, j]) / (dx * dx);
                //    //NewOmega[N / 2 + k, j] = NewPhi[N / 2 + k + 1, j - 1] / (dx * dx);
                //    k++;
                //}
            }
        private void NewSpeeds()
        {
            for (int i = 1; i < N + 1; i++)
            {
                for (int j = 1; j < M + 1; j++)
                {
                    if (triangle[i, j])
                    {
                        continue;
                    }
                    u[i, j] = (NewPhi[i, j + 1] - NewPhi[i, j]) / dy;
                    v[i, j] = -(NewPhi[i+1 , j] - NewPhi[i, j]) / dx;
                }
            }


        }
        public void culculate()
        {
            int IterOmega = 0;
            while (IterOmega <= MaxIterOmega)
            {
                NewVortex();
                if (!IsStacOmega()&&IterOmega>=MinIterOmega)
                {
                    break;
                }
                for (int i = 0; i < N + 2; i++)
                {
                    for(int j = 0; j < M + 2; j++)
                    {                        
                        NewPhi_t[i, j] = NewPhi[i, j];
                    }
                }
                int IterPhi = 0;
                while (IterPhi <= MaxIterPhi)
                {
                    NewCurrent();
                    if (!IsStacPhi()&&IterPhi>=MinIterPhi)
                    {
                        break;
                    }
                    //заполняем правую стенку
                    for (int j = 0; j < M + 2; j++)
                    {
                        NewPhi_t[N + 1, j] = NewPhi_t[N, j];
                    }

                    for (int i = 0; i < N + 2; i++)
                    {
                        for (int j = 0; j < M + 2; j++)
                        {
                            NewPhi[i, j] = NewPhi_t[i, j];
                        }
                    }
                    IterPhi++;

                }
                // Parallel.Invoke(NewBorder, NewSpeeds);
                NewBorder();
                NewSpeeds();
                for (int i = 0; i < N + 2; i++)
                {
                    for (int j = 0; j < M + 2; j++)
                    {
                        phi[i, j] = NewPhi[i, j];
                        omega[i, j] = NewOmega[i, j];
                    }
                }
                IterOmega++;

            }

            
        }

        private bool IsStacPhi()
        {
            for (int i = 1; i < N + 1; i++)
            {
                for(int j = 1; j < M + 1; j++)
                {
                    if (triangle[i, j])
                    {
                        continue;
                    }
                    if (Math.Abs(NewPhi[i, j] - NewPhi_t[i, j])>=epsPhi)
                    {
                        return true;
                    }
                }
            }
            return false;

        }
        private bool IsStacOmega()
        {
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < M + 2; j++)
                {
                    if (Math.Abs(NewOmega[i, j] - omega[i, j]) >= epsOmega)
                    {
                        return true;
                    }
                }
            }
            return false;
        }
       

       
        public double[,] get_u()
        {
            return u;
        }
        public double[,] get_v()
        {
            return v;
        }
        public double[,] get_phi()
        {
            return phi;
        }
        public double[,] get_omega()
        {
            return omega;
        }


    }
}
