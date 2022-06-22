using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    public struct Int : IComparable<Int>, ITranscedentals<Int>, IArithmetics<Int>
    {
        int i;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        public Int(int v = 0)
        {
            i = v;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        public static implicit operator Int(int v)=> new Int(v);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        public static explicit operator int(Int v) => v.i;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(Int other)
        {
            return i.CompareTo(other.i);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj)
        {
            Int item = (Int)obj;
            if (ReferenceEquals(item, null))
                return false;
            return i.Equals(item.i);
        }

        public override int GetHashCode()
        {
            return i.GetHashCode();
        }

        /// <summary>
        /// Converts a doub number to String for printing
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return i.ToString();
        }

        public string ToString4Array()
        {
            string myformat = "{0,9}";
            string ans = String.Format(myformat, i);
            return ans;
        }

        #region ITranscedental functions

        public double Abs()
        {
            return Math.Abs(i);
        }

        public double Arg()
        {
            return i >= 0 ? 0 : 3.1415926535897932384626433832795029;
        }

        public Int Sign()
        {
            return Math.Sign(i);
        }

        public Int Sin()
        {
            return (int)Math.Sin(i);
        }

        public Int Cos()
        {
            return (int)Math.Sin(i);
        }

        public Int Tan()
        {
            return (int)Math.Tan(i);
        }

        public Int Sinh()
        {
            return (int)Math.Sinh(i);
        }

        public Int Cosh()
        {
            return (int)Math.Cosh(i);
        }

        public Int Tanh()
        {
            return (int)Math.Tanh(i);
        }

        public Complex Asin()
        {
            Complex c = new Complex(i);
            return c.Asin();
        }

        public Complex Acos()
        {
            Complex c = new Complex(i);
            return c.Acos();
        }

        public Int Atan()
        {
            return (int)Math.Atan(i);
        }

        public Int Asinh()
        {
            Complex c = i + Math.Sqrt(i * i + 1);
            return (int)c.Log().Real;
        }

        public Complex Acosh()
        {
            Complex c = i * i - 1;
            return (i + c.Sqrt()).Log();
        }

        public Complex Atanh()
        {
            Complex q = (1 + i) / (1 - i);
            return 0.5 * q.Log();
        }

        public Complex Sqrt()
        {
            if (i >= 0) return Math.Sqrt(i);
            else
                return new Complex(0, Math.Sqrt(Math.Abs(i)));
        }

        public Int Sqrt(Int Other)
        {
            if (i >= 0) return (int)Math.Sqrt(Other.i);
            else
                return 0;
        }

        public Complex Log()
        {
            return new Complex(Math.Log(Abs()), Arg());
        }

        public Int Exp()
        {
            return (int)Math.Exp(i);
        }

        public Int Sqr()
        {
            return i*i;
        }

        public Int Atan2(Int x)
        {
            return (int)Math.Atan2(i, x.i);
        }

        public bool IsReal() => true;

        public dynamic REAL() => i;

        public dynamic IMAG() => 0;

        public Int Conj() => i;

        public Int One() => 1;

        public Int Zero() => 0;

        public Int BIGNO() => int.MaxValue;

        public Int TINYNO() => 0;

        public Int Rand() => Mth.rand.Next();

        public Int Round(int decpts) => i;


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
        public static Int operator -(Int A)
        { return -A.i; }

        public static Int operator -(Int A, Int B)
        { return A.i - B.i; }

        public static Int operator +(Int A, Int B)
        { return A.i + B.i; }

        public static Int operator *(Int A, Int B)
        { return A.i * B.i; }

        public static Int operator /(Int A, Int B)
        { return A.i / B.i; }

        public static Int operator *(Int A, bool b)
        { return b ? A : 0; }

        public static dynamic operator ^(Int A, Int b)
        {
            Complex c = A.i; Complex d = b.i;
            Complex c2d = c ^ d;
            if (c2d.IsReal())
                return c2d.Real;
            else
                return c2d;
        }

        public static bool operator ==(Int left, Int right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(Int left, Int right)
        {
            return !left.Equals(right);
        }

        public static bool operator >(Int A, Int B)
        { return A.CompareTo(B) > 0; }

        public static bool operator <(Int A, Int B)
        { return A.CompareTo(B) < 0; }

        public static bool operator >=(Int A, Int B)
        { return A.CompareTo(B) >= 0; }

        public static bool operator <=(Int A, Int B)
        { return A.CompareTo(B) <= 0; }

        #endregion

    }
}
