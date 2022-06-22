using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    public static class Mth
    {
        public static int DecimalPlace = 4;
        public static Doub PI = 3.1415926535897932384626433832795029,
                            E = 2.7182818284590452353602874713526625,
                            y = 0.5772156649015328606051209008240243,
                          Phi = 1.6180339887498948420458683436568311;
        public static Complex I = new Complex(0, 1);
        public static Random rand = new Random();

        #region Transcedentals
        #region Absolute
        public static Int Abs(Int i)
        { return (int)i.Abs(); }
        public static Doub Abs(Doub d)
        { return d.Abs(); }
        public static Doub Abs(Complex c)
        { return c.Abs(); }
        #endregion

        #region Square root
        public static dynamic Sqrt(Doub d)
        {
            Complex c = d.Sqrt();
            if (c.IsReal())return c.Real;
            else return c;
        }
        public static Complex Sqrt(Complex c)
        { return c.Sqrt(); }
        #endregion

        #region Sine
        public static Doub Sin(Doub d)
        { return d.Sin(); }
        public static Complex Sin(Complex c)
        { return c.Sin(); }
        #endregion

        #region Cosine
        public static Doub Cos(Doub d)
        { return d.Cos(); }
        public static Complex Cos(Complex c)
        { return c.Cos(); }
        #endregion

        #region Tangent
        public static Doub Tan(Doub d)
        { return d.Tan(); }
        public static Complex Tan(Complex c)
        { return c.Tan(); }
        #endregion

        #region ArcSine
        public static dynamic Asin(Doub d)
        {
            Complex c = d.Asin();
            if (c.IsReal()) return c.Real;
            else return c;
        }
        public static Complex Asin(Complex c)
        { return c.Asin(); }
        #endregion

        #region ArcCosine
        public static dynamic Acos(Doub d)
        {
            Complex c = d.Acos();
            if (c.IsReal()) return c.Real;
            else return c;
        }
        public static Complex Acos(Complex c)
        { return c.Acos(); }
        #endregion

        #region ArcTangent
        public static Doub Atan(Doub d)
        { return d.Atan(); }
        public static Complex Atan(Complex c)
        { return c.Atan(); }
        public static Doub Atan2(Doub y, Doub x)
        { return y.Atan2(x); }
        public static Complex Atan2(Complex y, Complex x)
        { throw new Exception("Argument must be real"); }
        #endregion

        #region Exponential
        public static Doub Exp(Doub d)
        { return d.Exp(); }
        public static Complex Exp(Complex c)
        { return c.Exp(); }
        #endregion

        #region Logarithm
        public static dynamic Log(Doub d)
        {
            if (d == 0)
                return double.NegativeInfinity;
            else
            {
                Complex c = d.Log();
                if (c.IsReal()) return c.Real;
                else return c;
            }
        }
        public static dynamic Log(Complex c)
        {
            if (Abs(c) == 0)
                return double.NegativeInfinity;
            else
                return c.Log();
        }
        #endregion

        #region SineH
        public static Doub Sinh(Doub d)
        { return d.Sinh(); }
        public static Complex Sinh(Complex c)
        { return c.Sinh(); }
        #endregion

        #region CosineH
        public static Doub Cosh(Doub d)
        { return d.Cosh(); }
        public static Complex Cosh(Complex c)
        { return c.Cosh(); }
        #endregion

        #region TangentH
        public static Doub Tanh(Doub d)
        { return d.Tanh(); }
        public static Complex Tanh(Complex c)
        { return c.Tanh(); }
        #endregion

        #region ArcSineH
        public static dynamic Asinh(Doub d)
        {
            Complex c = d.Asinh();
            if (c.IsReal()) return c.Real;
            else return c;
        }
        public static Complex Asinh(Complex c)
        { return c.Asinh(); }
        #endregion

        #region ArcCosineH
        public static dynamic Acosh(Doub d)
        {
            Complex c = d.Acosh();
            if (c.IsReal()) return c.Real;
            else return c;
        }
        public static Complex Acosh(Complex c)
        { return c.Acosh(); }
        #endregion

        #region ArcTangentH
        public static dynamic Atanh(Doub d)
        {
            Complex c = d.Atanh();
            if (c.IsReal()) return c.Real;
            else return c;
        }
        public static Complex Atanh(Complex c)
        { return c.Atanh(); }
        #endregion


        #region ArrayOperators
        public static int[] Stp(int start, int end)
        {
            int step = 1; int N = (end - start) / step + 1;
            int[] ans = new int[N]; ans[0] = start;
            for (int i = 1; i < N; i++) ans[i] = ans[i - 1] + step;
            return ans;
        }

        public static int[] Stp(int start, int step, int end)
        {
            int N = (end - start) / step + 1;
            int[] ans = new int[N]; ans[0] = start;
            for (int i = 1; i < N; i++) ans[i] = ans[i - 1] + step;
            return ans;
        }

        public static double[] Stp(double start, double end)
        {
            double step = 1; int N = (int)((end - start) / step + 1);
            double[] ans = new double[N]; ans[0] = start;
            for (int i = 1; i < N; i++) ans[i] = ans[i - 1] + step;
            return ans;
        }

        public static double[] Stp(double start, double step, double end)
        {
            int N = (int)((end - start) / step + 1);
            double[] ans = new double[N]; ans[0] = start;
            for (int i = 1; i < N; i++) ans[i] = ans[i - 1] + step;
            return ans;
        }

        public static dynamic Find(bool[] A, int k = 1)
        {
            List<int> ans = new List<int>();
            for (int i = 0; i < A.Length; i++)
            {
                if (A[i]) ans.Add(i);
            }
            if (ans.Count == 1) return ans[0];
            else return ans.ToArray();
        }

        public static dynamic Find(bool[,] A, int k = 1)
        {
            List<int> ans = new List<int>();
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    if (A[i, j]) ans.Add(i);
                }
            }
            if (ans.Count == 1) return ans[0];
            else return ans.ToArray();
        }
        #endregion

        #endregion

    }
}
