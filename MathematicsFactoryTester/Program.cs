using System;
using MathematicsFactory;
using System.Collections.Generic;
using System.Diagnostics;


namespace MathematicsFactoryTester
{
    class Program
    {
        
        static void Main(string[] args)
        {
            Mth.DecimalPlace = 4;
            Console.TreatControlCAsInput = false;

            //Matrix<Doub> A = Matrix<Doub>.Rand(7);
            //Console.WriteLine(A);
            //Console.WriteLine(A["", Stp(0, 2, 6)]);

            double v = SpecialFunctions.Erf(0.6);
            Console.WriteLine(v);

            //Matrix<Complex> A = new Complex[,] { { new Complex(99, -21), new Complex(23, 13), new Complex(-18, 17) },
            //                                     { new Complex(11, -12), new Complex(92, 51), new Complex(13, -27) },
            //                                     { new Complex(31, -14), new Complex(57, 20), new Complex(99, -89) } };

            //A.Choldc();
            //Console.WriteLine(A);
            //Console.WriteLine(A.L_chol);
            //Console.WriteLine(A.U_chol);

            //Matrix<Doub> B = new Doub[,] { { 9, 1, 1 }, { 1, 9, 1 }, { 1, 1, 9 } };
            //B.Choldc();
            //Console.WriteLine(B);
            //Console.WriteLine(B.L_chol);
            //Console.WriteLine(B.U_chol);

            //Matrix<Complex> C = new Complex[,] { { new Complex(99, -21), new Complex(23, 13), new Complex(-18, 17) },
            //                                     { new Complex(11, -12), new Complex(92, 51), new Complex(13, -27) },
            //                                     { new Complex(31, -14), new Complex(57, 20), new Complex(99, -89) } };
            //Mth.DecimalPlace = 4;
            //C.Choldc();
            //Console.WriteLine(C);
            //Console.WriteLine(C.L_chol);
            //Console.WriteLine(C.U_chol);
            //Console.WriteLine(C.Det());

            //Matrix<Complex> D = new Complex[,] { { new Complex(99, -21), new Complex(23, 13), new Complex(-18, 1) },
            //                                     { new Complex(23, 13), new Complex(92, 51), new Complex(13, -27) },
            //                                     { new Complex(-18, 1), new Complex(13, -27), new Complex(99, -89) } };
            //Mth.DecimalPlace = 4;
            //D.Choldc();
            //Console.WriteLine(D);
            //Console.WriteLine(D.L_chol);
            //Console.WriteLine(D.U_chol);
            //Console.WriteLine(D.Det());
            //Console.WriteLine(Matrix<Complex>.Mult(D.L_chol, D.U_chol.Conj()));

            //Matrix<Doub> E = new Doub[,] { { 99, 23, -18 },
            //                               { 23, 92,  13 },
            //                               { -18, 13, 99 } };
            //Mth.DecimalPlace = 4;
            //E.Choldc();
            //Console.WriteLine(E);
            //Console.WriteLine(E.L_chol);
            //Console.WriteLine(E.U_chol);
            //Console.WriteLine(E.Det());
            //Console.WriteLine(E.L_chol*E.U_chol);

            //Matrix<Complex> E = new Complex[,] { { new Complex(99, -21), new Complex(23, 13), new Complex(-18, 1) },
            //                                     { new Complex(23, 13), new Complex(92, 51), new Complex(13, -27) },
            //                                     { new Complex(-18, 1), new Complex(13, -27), new Complex(99, -89) } };
            //Mth.DecimalPlace = 4;
            //E.MakeLU();
            //Console.WriteLine(E);
            //Console.WriteLine(E.L_lu);
            //Console.WriteLine(E.U_lu);


            //Matrix<Doub> F = new Doub[,] { { 10, -1,  0,  0,  0,  0,  0},
            //                               { -1, 10, -1,  0,  0,  0,  0},
            //                               {  0, -1, 10, -1,  0,  0,  0},
            //                               {  0,  0, -1, 10, -1,  0,  0},
            //                               {  0,  0,  0, -1, 10, -1,  0},
            //                               {  0,  0,  0,  0, -1, 10, -1},
            //                               {  0,  0,  0,  0,  0, -1, 10}};
            //Matrix<Doub> v = new Doub[,] { { 1, 2, 3},
            //                               { 1, 2, 3},
            //                               { 1, 2, 3},
            //                               { 1, 9, 3},
            //                               { 1, 2, 3},
            //                               { 1, 2, 3},
            //                               { 1, 2, 3}};
            //Mth.DecimalPlace = 4;
            //Console.WriteLine(Matrix<Doub>.TriSolve(F, v));


            //Matrix<Doub> G = new Doub[,] { { 10, -1,  0,  0,  0,  0,  0},
            //                               { -1, 10, -1,  0,  0,  0,  0},
            //                               {  0, -1, 10, -1,  0,  0,  0},
            //                               {  0,  0, -1, 10, -1,  0,  0},
            //                               {  0,  0,  0, -1, 10, -1,  0},
            //                               {  0,  0,  0,  0, -1, 10, -1},
            //                               {  0,  0,  0,  0,  0, -1, 10}};
            //Mth.DecimalPlace = 4;
            //Console.WriteLine(Matrix<Doub>.TriSolve(F, Matrix<Doub>.Eye(7)));


            //Matrix<Complex> G = new Complex[,] { { new Complex(99, -21), new Complex(23, 13), 0 },
            //                                     { new Complex(23, 13), new Complex(92, 51), new Complex(13, -27) },
            //                                     { 0, new Complex(13, -27), new Complex(99, -89) } };
            //Mth.DecimalPlace = 4;
            //Console.WriteLine(G);
            //Matrix<Complex> Ginv = Matrix<Complex>.TriSolve(G, Matrix<Complex>.Eye(3));
            //Console.WriteLine(Matrix<Complex>.Mult(Ginv, G));


            //Doub[,] Hmt = new Doub[10,10];
            //Doub[] Hb = new Doub[10];
            //Random rand = new Random();
            //for (int i = 0; i < 10; i++)
            //{
            //    if (i > 0) Hmt[i, i - 1] = rand.NextDouble();
            //    if (i < 9) Hmt[i, i + 1] = rand.NextDouble();
            //    Hmt[i, i] = rand.NextDouble();
            //    Hb[i] = rand.NextDouble();
            //}
            //Matrix<Doub> H = Hmt, Hsol;
            //Mth.DecimalPlace = 4;
            //Stopwatch Timer = new Stopwatch();
            //Timer.Start();
            //for (int k = 0; k < 1000000; k++) Hsol = Matrix<Doub>.TriSolve(H, Hb);
            //Timer.Stop();
            //Console.WriteLine(H);
            //Console.WriteLine(Hsol);

            //Complex[,] Kmt = new Complex[3000, 3000];
            //Complex[] Kb = new Complex[3000];
            //Random rand = new Random();
            //for (int i = 0; i < 3000; i++)
            //{
            //    if (i > 0) Kmt[i, i - 1] = Complex.Random(rand);
            //    if (i < 2999) Kmt[i, i + 1] = Complex.Random(rand);
            //    Kmt[i, i] = Complex.Random(rand);
            //    Kb[i] = Complex.Random(rand);
            //}
            //Matrix<Complex> K = Kmt;
            //Mth.DecimalPlace = 4;
            //Matrix<Complex> Ksol = Matrix<Complex>.TriSolve(H, Hb);
            //Console.WriteLine(Hsol);


            //Matrix<Complex> K = new Complex[,] { { new Complex(99, -21), new Complex(23, 13), new Complex(-18, 1) },
            //                                     { new Complex(23, 13), new Complex(92, 51), new Complex(13, -27) },
            //                                     { new Complex(-18, 1), new Complex(13, -27), new Complex(99, -89) } };
            //Mth.DecimalPlace = 4;
            //K.MakeLU();
            //Console.WriteLine(K);
            //Console.WriteLine(K.L_lu);
            //Console.WriteLine(K.U_lu);
            //Console.WriteLine(Matrix<Complex>.Solve(K, Matrix<Complex>.Eye(3)));

            //Mth.DecimalPlace = 4;
            //Matrix<Complex> L = Matrix<Complex>.Round(Matrix<Complex>.Rand(5, 5),4);
            //Console.WriteLine(L);
            //Console.WriteLine(Matrix<Complex>.Solve(L, Matrix<Complex>.Eye(5)));

            //Matrix<Doub> M= Matrix<Doub>.Add(Matrix<Doub>.Rand(5, 5), Matrix<Doub>.Rand(5, 5));
            //var Ms = Matrix <Doub>.Asin(M);
            //Console.WriteLine(M);
            //Console.WriteLine(Ms);

            //Matrix<Doub> N = Matrix<Doub>.Rand(5, 5);
            //var Ns = Matrix<Doub>.Asin(N);
            //Console.WriteLine(N);
            //Console.WriteLine(Ns);

            //Matrix<Doub> O = Matrix<Doub>.Rand(3);
            //Matrix<Doub> P = Matrix<Doub>.Rand(3);
            //Console.WriteLine(O);
            //Console.WriteLine(P);
            //Console.WriteLine(P.DotMult(O));

            //Matrix<Int> x = Stp(1, 6), y = Stp(1, 8);
            //Matrix<Int>[] XY = Matrix<Int>.Meshgrid(x, y);
            //Console.WriteLine(XY[0]);
            //Console.WriteLine(XY[1]);


            //Matrix<Doub> A = Matrix<Doub>.Rand(10);
            //Matrix <Int> I = new Int[,] { { 1, 2 }, { 3, 4 } };
            //Console.WriteLine(A);
            //Console.WriteLine(A[I,I]);

            //Matrix<Doub> A = Matrix<Doub>.Rand(6);
            //Matrix<Int> I = new Int[,] { { 1, 2 }, { 3, 4 } };
            //Console.WriteLine(A);
            //Console.WriteLine(A[I]);
            //Console.WriteLine(A[I+1]);


            //Matrix<Doub> Points = new Doub[,] { { 0.44759,         0.27961,         0.13765 },
            //                                    { 0.55758,         0.71711,         0.81953 },
            //                                    { 0.44433,         0.32083,         0.81114 },
            //                                    { 0.60475,         0.60659,         0.69274 },
            //                                    { 0.79989,         0.94045,         0.20218 },
            //                                    { 0.37050,         0.38219,         0.29579 },
            //                                    { 0.04182,         0.26368,         0.18067 },
            //                                    { 0.20866,         0.57174,         0.58434 },
            //                                    { 0.47552,         0.67745,         0.48055 },
            //                                    { 0.92824,         0.97647,         0.27657 },
            //                                    { 0.89711,         0.57505,         0.80801 },
            //                                    { 0.17432,         0.19620,         0.31602 },
            //                                    { 0.34450,         0.11481,         0.00413 },
            //                                    { 0.64081,         0.56888,         0.83586 },
            //                                    { 0.70656,         0.90002,         0.92896 },
            //                                    { 0.07741,         0.49129,         0.41020 },
            //                                    { 0.06571,         0.36607,         0.45184 },
            //                                    { 0.55905,         0.25269,         0.71434 },
            //                                    { 0.54724,         0.23893,         0.70786 },
            //                                    { 0.55331,         0.08205,         0.03967 }};
            //Delaunay DT = new Delaunay(Points);
        }


        static Doub[] Stp(Doub start, Doub end)
        {
            Doub step = 1; int N =(int)((end - start) / step + 1);
            Doub[] ans = new Doub[N]; ans[0] = start;
            for (int i = 1; i < N; i++) ans[i] = ans[i - 1] + step;
            return ans;
        }

        static Int[] Stp(Int start, Int end)
        {
            Int step = 1; int N = (int)((end - start) / step + 1);
            Int[] ans = new Int[N]; ans[0] = start;
            for (int i = 1; i < N; i++) ans[i] = ans[i - 1] + step;
            return ans;
        }

        static int[] Stp(int start, int step, int end)
        {
            List<int> ans = new List<int>(); ans.Add(start);
            while (start < end)
            { start += step; ans.Add(start); }
            return ans.ToArray();
        }

    }
}
