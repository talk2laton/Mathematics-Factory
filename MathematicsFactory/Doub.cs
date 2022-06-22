using System;

namespace MathematicsFactory
{
    public struct Doub : IComparable<Doub>, ITranscedentals<Doub>, IArithmetics<Doub>
    {
        double d;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        public Doub (double v = 0)
        {
            d = v;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        public static implicit operator Doub(Int v)=> new Doub((int)v);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        public static implicit operator Doub(double v)=> new Doub(v);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        public static explicit operator double(Doub v)=> v.d;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(Doub other)
        {
            return d.CompareTo(other.d);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj)
        {
            Doub item = (Doub)obj;
            if (ReferenceEquals(item, null))
                return false;
            return d.Equals(item.d);
        }

        public override int GetHashCode()
        {
            return d.GetHashCode();
        }

        /// <summary>
        /// Converts a doub number to String for printing
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            int Dp = Mth.DecimalPlace, decplace = 0;
            if (d != 0) decplace = (int)Math.Floor(Math.Log(Math.Abs(d))); 
            string pad = "0000000000000000000000"; string myformat;
            if (Math.Abs(decplace) >=2)
                myformat = "{0," + (4 + Dp) + ":0." + pad.Substring(0, Dp) + "e+00}";
            else
                myformat = "{0," + (4 + Dp) + ":0." + pad.Substring(0, Dp)+"}";
            string ans = String.Format(myformat, d);
            return ans;
        }

        public string ToString4Array()
        {
            int Dp = Mth.DecimalPlace; string pad = "0000000000000000000000", myformat;
            myformat = "{0," + (4 + Dp) + ":0." + pad.Substring(0, Dp) + "}";
            string ans = String.Format(myformat, d);
            return ans;
        }

        #region ITranscedental functions

        public double Abs()
        {
            return Math.Abs(d);
        }

        public double Arg()
        {
            return d >= 0 ? 0 : 3.1415926535897932384626433832795029;
        }

        public Doub Sign()
        {
            return Math.Sign(d);
        }

        public Doub Sin()
        {
            return Math.Sin(d);
        }

        public Doub Cos()
        {
            return Math.Cos(d);
        }

        public Doub Tan()
        {
            return Math.Tan(d);
        }

        public Doub Sinh()
        {
            return Math.Sinh(d);
        }

        public Doub Cosh()
        {
            return Math.Cosh(d);
        }

        public Doub Tanh()
        {
            return Math.Tanh(d);
        }

        public Complex Asin()
        {
            Complex c = new Complex(d);
            return c.Asin();
        }

        public Complex Acos()
        {
            Complex c = new Complex(d);
            return c.Acos();
        }

        public Doub Atan()
        {
            return Math.Atan(d);
        }

        public Doub Asinh()
        {
            Complex c = d + Math.Sqrt(d * d + 1);
            return c.Log().Real;
        }

        public Complex Acosh()
        {
            Complex c = d * d - 1;
            return (d + c.Sqrt()).Log();
        }

        public Complex Atanh()
        {
            Complex q = (1 + d) / (1 - d);
            return 0.5 * q.Log();
        }

        public Complex Sqrt()
        {
            if (d >= 0) return Math.Sqrt(d);
            else
                return new Complex(0, Math.Sqrt(Math.Abs(d)));
        }

        public Doub Sqrt(Doub Other)
        {
            if (d >= 0) return Math.Sqrt(Other.d);
            else
                return 0;
        }

        public Complex Log()
        {
            return new Complex(Math.Log(Abs()), Arg());
        }

        public Doub Exp()
        {
            return Math.Exp(d);
        }

        public Doub Sqr()
        {
            return d * d;
        }

        public Doub Atan2(Doub x)
        {
            return Math.Atan2(d, x.d);
        }

        public bool IsReal() => true;

        public dynamic REAL() => d;

        public dynamic IMAG() => 0;

        public Doub Conj() => d;

        public Doub One() => 1;

        public Doub Zero() => 0;

        public Doub BIGNO() => 1e30;

        public Doub TINYNO() => 1e-30;

        public Doub Rand() => Mth.rand.NextDouble();

        public Doub Round(int decpts) => Math.Round(d, decpts);

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

        #region Operators
        public static Doub operator -(Doub A)
        { return -A.d; }

        public static Doub operator -(Doub A, Doub B)
        { return A.d - B.d; }

        public static Doub operator +(Doub A, Doub B)
        { return A.d + B.d; }

        public static Doub operator *(Doub A, Doub B)
        { return A.d*B.d; }

        public static Doub operator /(Doub A, Doub B)
        { return A.d / B.d; }

        public static Doub operator *(Doub A, bool b)
        { return b ? A : 0; }

        public static dynamic operator ^(Doub A, Doub b)
        {
            Complex c = A; Complex d = b;
            Complex c2d = c ^ d;
            if (c2d.IsReal())
                return c2d.Real;
            else
                return c2d;
        }

        public static bool operator ==(Doub left, Doub right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(Doub left, Doub right)
        {
            return !left.Equals(right);
        }

        public static bool operator >(Doub A, Doub B)
        { return A.CompareTo(B) > 0; }

        public static bool operator <(Doub A, Doub B)
        { return A.CompareTo(B) < 0; }

        public static bool operator >=(Doub A, Doub B)
        { return A.CompareTo(B) >= 0; }

        public static bool operator <=(Doub A, Doub B)
        { return A.CompareTo(B) <= 0; }

        #endregion
    }
}
