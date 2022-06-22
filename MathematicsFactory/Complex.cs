using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    public struct Complex : IComparable<Complex>, ITranscedentals<Complex>, IArithmetics<Complex>
    {
        public Doub Real { get { return x; } }
        public Doub Imaginary { get { return y; } }
        double x, y;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="real"></param>
        /// <param name="imaginary"></param>
        public Complex(Doub real, Doub imaginary)
        {
            x = (double)real; 
            y = (double)imaginary;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="real"></param>
        /// <param name="imaginary"></param>
        public Complex(double real, double imaginary)
        {
            x = real;
            y = imaginary;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="real"></param>
        /// <param name="imaginary"></param>
        public Complex(double real)
        {
            x = real;
            y = 0;
        }

        /// <summary>
        /// Implicit operator for conversion from double to complex
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static implicit operator Complex(Doub value)
        {
            return new Complex(value, 0);
        }

        /// <summary>
        /// Implicit operator for conversion from double to complex
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static implicit operator Complex(Int value)
        {
            return new Complex(value, 0);
        }

        /// <summary>
        /// Implicit operator for conversion from double to complex
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static implicit operator Complex(double value)
        {
            return new Complex(value, 0);
        }

        /// <summary>
        /// Converts a Complex number to String for printing
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            if (IsReal())
            {
                return ((Doub)x).ToString();
            }
            else
            {
                double ax, ay, xy = Math.Max(ax = Math.Abs(x), ay = Math.Abs(y));
                int Dp = Mth.DecimalPlace, decplace = 0;
                if (xy != 0) decplace = Math.Abs((int)Math.Floor(Math.Log10(xy)));
                string pad = "000000000000000000000000000000000", fmtx, fmty, ans;
                if (decplace >= 2)
                {
                    fmtx = "{0," + (4 + Dp) + ":0." + pad.Substring(0, Dp) + "e+00}";
                    fmty = "{0," + (-(4 + Dp)) + ":0." + pad.Substring(0, Dp) + "e+00i}";
                }
                else
                {
                    fmtx = "{0," + (4 + Dp) + ":0." + pad.Substring(0, Dp) + "}";
                    fmty = "{0," + (-(4 + Dp)) + ":0." + pad.Substring(0, Dp) + "i}";
                }
                ans = String.Format(fmtx, x) + (y >= 0.0 ? " + " : " - ") + String.Format(fmty, ay);
                return ans;
            }
        }

        public string ToString4Array()
        {
            double ay = Math.Abs(y); int Dp = Mth.DecimalPlace, F = 4 + Dp;
            string pad = "0000000000000000000000", zeros = pad.Substring(0, Dp), fmtx, fmty, ans;
            fmtx = "{0," + (4 + Dp) + ":0." + zeros + "}"; fmty = "{0," + (-(4 + Dp)) + ":0." + zeros + "i}";
            ans = String.Format(fmtx, x) + (y >= 0.0 ? " + " : " - ") + String.Format(fmty, ay);
            return ans;
        }

        public int CompareTo(Complex other)
        {
            if (Abs() == other.Abs())
                return Arg().CompareTo(other.Arg());
            return Abs().CompareTo(other.Abs());
        }

        public override bool Equals(object obj)
        {
            Complex item = (Complex)obj;
            if (ReferenceEquals(item, null))
                return false; 
            if (Abs() == item.Abs())
                return Arg().Equals(item.Arg());
            else
                return Abs().Equals(this.Abs());
        }

        public override int GetHashCode()
        {
            return Abs().GetHashCode();
        }
        
        

        /// <summary>
        /// 
        /// </summary>
        /// <param name="abs"></param>
        /// <param name="angle"></param>
        /// <returns></returns>
        public static Complex Cart(double abs, double angle)
        {
            return new Complex(abs * Math.Cos(angle), abs * Math.Sin(angle));
        }

        public static double[] Pol(Complex c)
        {
            double[] ans = new double[2];
            ans[0] = c.Abs();
            ans[1] = c.Arg();
            return ans;
        }

        public static Complex[] Root(Complex c, int N)
        {
            double pi, ang, angN, abs, absN, del; 
            pi = 3.1415926535897932384626433832795;
            ang = c.Arg(); abs = c.Abs(); angN = ang / N;
            absN = Math.Pow(abs, (1.0 / N)); del = 2 * pi / N;
            Complex[] ans = new Complex[N];
            for (int n = 0; n < N; n++)
            {
                angN += del;
                ans[n] = Cart(absN, angN);
            }
            return ans;
        }

        #region ITranscedental functions

        public double Abs()
        {
            double ax, ay, axy, ayx;
            if ((ax = Math.Abs(x)) >= (ay = Math.Abs(y)))
            {
                ayx = ay / ax;
                return ax * Math.Sqrt(1 + ayx * ayx);
            }
            else
            {
                axy = ax / ay;
                return ay * Math.Sqrt(1 + axy * axy);
            }
        }

        public double Arg()
        {
            return Math.Atan2(y, x);
        }

        public Complex Sign()
        {
            double abs = Abs();
            return new Complex( x/abs, y/abs);
        }

        public Complex Sin()
        {
            Complex i2, cp, cm;
            i2 = new Complex(0, 2);
            cp = new Complex(-y, x);
            cm = new Complex(y, -x);
            return (cp.Exp() - cm.Exp()) / i2;
        }

        public Complex Cos()
        {
            Complex cp, cm;
            cp = new Complex(-y, x);
            cm = new Complex(y, -x);
            return (cp.Exp() + cm.Exp()) / 2;
        }

        public Complex Tan()
        {
            return Sin() / Cos();
        }

        public Complex Sinh()
        {
            Complex cp, cm;
            cp = new Complex(x, y);
            cm = new Complex(-x, -y);
            return (cp.Exp() - cm.Exp()) / 2;
        }

        public Complex Cosh()
        {
            Complex cp, cm;
            cp = new Complex(x, y);
            cm = new Complex(-x, -y);
            return (cp.Exp() + cm.Exp()) / 2;
        }

        public Complex Tanh()
        {
            return Sinh() / Cosh();
        }

        public Complex Asin()
        {
            Complex ci, ci2p1, d, e;
            ci = new Complex(-y, x);
            ci2p1 = ci * ci + 1;
            d = ci2p1.Sqrt() + ci;
            e = d.Log();
            return new Complex(e.y, -e.x);
        }

        public Complex Acos()
        {
            Complex c = this, ci, ci2p1, d, di, dipc, e;
            ci = new Complex(-y, x);
            ci2p1 = ci * ci + 1;
            d = ci2p1.Sqrt();
            di = new Complex(-d.y, d.x);
            dipc = di + c;
            e = dipc.Log();
            return new Complex(e.y, -e.x);
        }

        public Complex Atan()
        {
            Complex OnepC, OnemC, q, e;
            OnepC = new Complex(1 - y, x);
            OnemC = new Complex(1 + y, -x);
            q = OnemC / OnepC; e = q.Log();
            return new Complex(-0.5*e.y, 0.5*e.x);
        }

        public Complex Asinh()
        {
            Complex c = this, c2p1 = c * c + 1; 
            return (c + c2p1.Sqrt()).Log();
        }

        public Complex Acosh()
        {
            Complex c = this, c2m1 = c * c - 1;
            return (c + c2m1.Sqrt()).Log();
        }

        public Complex Atanh()
        {
            Complex OnemC, OnepC, q;
            OnemC = new Complex(1 - x, -y);
            OnepC = new Complex(1 + x, y);
            q = OnepC / OnemC;
            return 0.5 * q.Log();
        }

        public Complex Sqrt()
        {
            double ax, ay = 0, ax2, ay2, w = 0;
            if (x == y)
            { w = 0; }
            else if ((ax = Math.Abs(x)) >= (ay = Math.Abs(y)))
            {
                ax2 = ax * ax; ay2 = ay * ay;
                w = Math.Sqrt(ax * 0.5 * (1 + Math.Sqrt(1 + (ay2 / ax2))));
            }
            else
            {
                ax2 = ax * ax; ay2 = ay * ay;
                w = Math.Sqrt(ay * 0.5 * ((ax / ay) + Math.Sqrt(1 + (ax2 / ay2))));
            }

            if (w == 0) return 0;
            else
            {
                if (x >= 0) return new Complex(w, y / (2 * w));
                else
                {
                    if (y >= 0) return new Complex(ay / (2 * w), w);
                    else return new Complex(ay / (2 * w), -w);
                }
            }
        }

        public Complex Sqrt(Complex Other) => Other.Sqrt();

        public Complex Log()
        {
            return new Complex(Math.Log(Abs()), Arg());
        }

        public Complex Exp()
        {
            return Cart(Math.Exp(x), y);
        }

        public Complex Sqr()
        {
            return this*this;
        }
        
        public bool IsReal()
        {
            if (x == 0 && y == 0)
                return true;
            return Math.Abs(y / Abs()) < 1e-14;
        }

        public bool IsComplex()
        {
            return !IsReal();
        }

        public bool IsImaginary()
        {
            if (x == 0 && y == 0)
                return false;
            return Math.Abs(x / Abs()) < 1e-14;
        }

        public dynamic REAL() => x;

        public dynamic IMAG() => y;

        public Complex Conj() => new Complex(x, -y);

        public Complex One() => 1;

        public Complex Zero() => 0;

        public Complex BIGNO() => 1e30;

        public Complex TINYNO() => 1e-30;

        public Complex Rand() => 
            new Complex(Mth.rand.NextDouble(), Mth.rand.NextDouble());

        public Complex Round(int decpts) => 
            new Complex(Math.Round(x, decpts), Math.Round(y, decpts));


        #endregion

        #region IArithmetics functions
        public dynamic Add(Complex other) => this + other;
        public dynamic Add(Doub other) => this + other;
        public dynamic Add(Int other) => this + other;

        public dynamic Mult(Complex other) => this * other;
        public dynamic Mult(Doub other) => this * other;
        public dynamic Mult(Int other) => this * other;

        public dynamic Sub(Complex other) => this - other;
        public dynamic Sub(Doub other) => this - other;
        public dynamic Sub(Int other) => this - other;

        public dynamic Div(Complex other) => this / other;
        public dynamic Div(Doub other) => this / other;
        public dynamic Div(Int other) => this / other;

        public dynamic Pow(Complex other) => this ^ other;
        public dynamic Pow(Doub other) => this ^ other;
        public dynamic Pow(Int other) => this ^ other;
        #endregion

        #region Operators workhorse
        private static Complex Conjugate(Complex c)
        {
            return new Complex(c.Real, -c.Imaginary);
        }

        private static Complex Add(Complex c1, Complex c2)
        {
            return new Complex(c1.x + c2.x, c1.y + c2.y);
        }

        private static Complex Sub(Complex c1, Complex c2)
        {
            return new Complex(c1.x - c2.x, c1.y - c2.y);
        }

        private static Complex Mult(Complex c1, Complex c2)
        {
            double a = c1.x, b = c1.y, c = c2.x, d = c2.y;
            double ac = a * c, bd = b * d, apb = a + b, cpd = c + d;
            return new Complex(ac - bd, apb * cpd - ac - bd);
        }
        
        private static Complex Div(Complex c1, Complex c2)
        {
            double a = c1.x, b = c1.y, c = c2.x, d = c2.y, numr, numi, denr, dc, cd;
            if (Math.Abs(c) >= Math.Abs(d))
            { dc = d / c; numr = a + b * dc; numi = b - a * dc; denr = c + d * dc; }
            else
            { cd = c / d; numr = a * cd + b; numi = b * cd - a; denr = c * cd + d; }
            return new Complex(numr / denr, numi / denr);
        }

        private static Complex Pow(Complex c, double N)
        {
            double ang, angN, abs, absN;
            ang = c.Arg(); abs = c.Abs(); angN = ang * N;
            absN = Math.Pow(abs, N);
            Complex ans = Cart(absN, angN);
            return ans;
        }

        private static Complex Pow(Complex c1, Complex c2)
        {
            return (c2 * c1.Log()).Exp();
        }
        #endregion

        #region operators

        public static dynamic operator +(Complex c1, Complex c2)
        {
            Complex c = Add(c1, c2);
            if (c.IsReal())
                return c.Real;
            else
                return c;
        }

        public static dynamic operator -(Complex c1)
        {
            Complex c = new Complex(-c1.x, -c1.y);
            if (c.IsReal())
                return c.Real;
            else
                return c;
        }
        
        public static dynamic operator -(Complex c1, Complex c2)
        {
            Complex c = Sub(c1, c2);
            if (c.IsReal())
                return c.Real;
            else
                return c;
        }

        public static dynamic operator *(Complex c1, Complex c2)
        {
            Complex c = Mult(c1, c2);
            if (c.IsReal())
                return c.Real;
            else
                return c;
        }
        
        public static dynamic operator /(Complex c1, Complex c2)
        {
            Complex c = Div(c1, c2);
            if (c.IsReal())
                return c.Real;
            else
                return c;
        }

        public static dynamic operator ~(Complex c1)
        {
            Complex c = Conjugate(c1);
            if (c.IsReal())
                return c.Real;
            else
                return c;
        }
        
        public static dynamic operator ^(Complex c1, Complex c2)
        {
            return Pow(c1, c2);
        }
        
        public static bool operator ==(Complex left, Complex right)
        {
            return left.Equals(right);
        }
        
        public static bool operator !=(Complex left, Complex right)
        {
            return !left.Equals(right);
        }

        #endregion
    }

    /// <summary>
    /// Exception from the Complex Class
    /// </summary>
    public class CException : Exception
    {
        /// <summary>
        /// Exception from the Complex Class
        /// </summary>
        /// <param name="Message">Message to be displayed</param>
        public CException(string Message)
            : base(Message)
        { }
    }
}
