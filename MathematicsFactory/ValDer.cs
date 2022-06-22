using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    public class ValDer
    {
        public dynamic Val { get { return val; } set { val = value; } }
        public dynamic Der { get { return der; } set { der = value; } }
        dynamic val, der;
        private int depth = 1;
        int Depth { get { return MakeDepth(); } }
        public bool IsConstant { get; set; }

        public ValDer(dynamic Vl, dynamic Dr)
        {
            val = Vl;
            der = Dr;
            IsConstant = ReferenceEquals(der, 0);
        }

        public int MakeDepth()
        {
            dynamic test = this.val;
            Type unknown = test.GetType();
            bool Isprimitive = unknown == typeof(Complex) ||
                               unknown == typeof(double) ||
                               unknown == typeof(int);
            int dp = depth;
            while (!Isprimitive)
            {
                dp++;
                test = test.val;
                unknown = test.GetType();
                Isprimitive = unknown == typeof(Complex) ||
                              unknown == typeof(double) ||
                              unknown == typeof(int);
            }
            return dp;
        }

        public dynamic Derivative(int order)
        {
            int N = Depth;
            dynamic ans = this;
            while (N > 0)
            { ans = N > (order) ? ans.val : ans.der; N--; }
            return ans;
        }

        public ValDer(dynamic Val)
        {
            val = Val;
            der = 0;
            IsConstant = der == 0;
        }

        /// <summary>
        /// Implicit conversion from type double to ValDer
        /// </summary>
        /// <param name="value"></param>
        static public implicit operator ValDer(double value)
        {
            return new ValDer(value);
        }

        /// <summary>
        /// Explicit conversion from type double to ValDer
        /// </summary>
        /// <param name="value"></param>
        static public explicit operator double(ValDer value)
        {
            return value.val;
        }

        /// <summary>
        /// This method prepares the information for diaply by Console.writeline or write 
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string s1 = "f", s2 = "'";
            int N = this.Depth;
            string s = s1 + " = " + Derivative(0) + "\n";
            for (int i = 1; i < Depth; i++)
            {
                s1 = s1 + s2;
                s = s + s1 + " = " + Derivative(i) + "\n";
            }
            s1 = s1 + s2;
            s = s + s1 + " = " + Derivative(Depth);
            return s;
        }

        /// <summary>
        /// Method performs addition with automatic differentiation
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        public static ValDer vdAdd(ValDer m1, ValDer m2)
        {
            dynamic Val = m1.val + m2.val;
            dynamic Der = m2.der + m1.der;
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method performs multiplication with automatic differentiation
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        public static ValDer vdTimes(ValDer m1, ValDer m2)
        {
            dynamic Val = m1.val * m2.val;
            dynamic Der = m1.val * m2.der + m1.der * m2.val;
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method performs division with automatic differentiation
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        public static ValDer vdDivide(ValDer m1, ValDer m2)
        {
            dynamic Val = m1.val / m2.val;
            dynamic Der = (m2.val * m1.der - m1.val * m2.der) / (m2.val * m2.val);
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method raise a number to a given Power with automatic differentiation
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        public static ValDer vdPow(ValDer m1, ValDer m2)
        {
            ValDer ans = new ValDer(0);
            if (m2.IsConstant)
            {
                ans.val = pow(m1.val, m2.val);
                ans.der = m1.der * m2.val * pow(m1.val, m2.val - 1);
            }
            else
            {
                ans = exp(m2 * Log(m1));
            }
            return ans;
        }

        /// <summary>
        /// Method computes the squareroot of a number with automatic differentiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Sqrt(ValDer m)
        {
            dynamic Val = sqrt(m.val);
            dynamic Der = 0.5 * m.der / sqrt(m.val);
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the Sine of a number with automatic differentiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Sin(ValDer m)
        {
            dynamic Val = sin(m.val);
            dynamic Der = cos(m.val) * m.der;
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the Cosine of a number with automatic differentiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Cos(ValDer m)
        {
            dynamic Val = cos(m.val);
            dynamic Der = -sin(m.val) * m.der;
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the tangent of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Tan(ValDer m)
        {
            dynamic Val = tan(m.val);
            dynamic Der = m.der / (cos(m.val) * cos(m.val));
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the hyperbolic sine of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Sinh(ValDer m)
        {
            dynamic Val = sinh(m.val);
            dynamic Der = cosh(m.val) * m.der;
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the hyperbolic cosine of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Cosh(ValDer m)
        {
            dynamic Val = cosh(m.val);
            dynamic Der = sinh(m.val) * m.der;
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the hyperbolic tangent of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Tanh(ValDer m)
        {
            dynamic Val = tanh(m.val);
            dynamic Der = m.der / (cosh(m.val) * cosh(m.val));
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the arcsine of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Asin(ValDer m)
        {
            dynamic Val = asin(m.val);
            dynamic Der = m.der / Sqrt(1 - m.val * m.val);
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the arccosine of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Acos(ValDer m)
        {
            dynamic Val = acos(m.val);
            dynamic Der = -m.der / sqrt(1 - m.val * m.val);
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the arctan of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Atan(ValDer m)
        {
            dynamic Val = atan(m.val);
            dynamic Der = m.der / (1 + m.val * m.val);
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the natural logarithm of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Ln(ValDer m)
        {
            dynamic Val = Log(m.val);
            dynamic Der = m.der / (m.val);
            return new ValDer(Val, Der);
        }

        /// <summary>
        /// Method computes the logarithm of a number to a base with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Log(ValDer vd, double n)
        {
            return Ln(vd) / Math.Log(n);
        }

        /// <summary>
        /// Method computes the exponent of a number with automatic differntiation
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static ValDer Exp(ValDer m)
        {
            dynamic Val = exp(m.val);
            dynamic Der = m.der * Val;
            return new ValDer(Val, Der);
        }

        public static ValDer Asinh(ValDer m)
        {
            dynamic Val = asinh(m.val);
            dynamic Der = m.der / sqrt(m.val * m.val + 1);
            return new ValDer(Val, Der);
        }

        public static ValDer Acosh(ValDer m)
        {
            dynamic Val = acosh(m.val);
            dynamic Der = m.der / sqrt(m.val * m.val - 1);
            return new ValDer(Val, Der);
        }

        public static ValDer Atanh(ValDer m)
        {
            dynamic Val = atanh(m.val);
            dynamic Der = m.der / (m.val * m.val - 1);
            return new ValDer(Val, Der);
        }


        #region private transcedental function

        static dynamic sqrt(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Sqrt();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Sqrt(c);
            }
            else
            {
                if (c > 0)
                { return Math.Sqrt(c); }
                else
                { return c.Sqrt(); }
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic Log(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Log();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Log(c);
            }
            else
            {
                if (c > 0)
                { return Math.Log(c); }
                else
                { return c.Log(); }
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic pow(dynamic c, dynamic n)
        {
            Type unknown1 = c.GetType();
            Type unknown2 = n.GetType();
            if (unknown2 == typeof(Complex))
            {
                if (unknown1 == typeof(Complex))
                {
                    return c ^ n;
                }
                else if (unknown1 == typeof(ValDer))
                {
                    return vdPow(c, n);
                }
                else
                {
                    return new Complex(c) ^ n;
                }
            }
            else if (unknown2 == typeof(ValDer))
            {
                if (unknown1 == typeof(Complex))
                {
                    return Exp(c.Log() * n);
                }
                else if (unknown1 == typeof(ValDer))
                {
                    return vdPow(c, n);
                }
                else
                {
                    return Exp(Math.Log(c) * n);
                }
            }
            else
            {
                if (unknown1 == typeof(Complex))
                {
                    return c ^ new Complex(n);
                }
                else if (unknown1 == typeof(ValDer))
                {
                    return vdPow(c, n);
                }
                else
                {
                    return Math.Pow(c, n);
                };
            }
        }

        static dynamic exp(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Exp();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Exp(c);
            }
            else
            {
                return Math.Exp(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic sin(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Sin();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Sin(c);
            }
            else
            {
                return Math.Sin(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic cos(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Cos();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Cos(c);
            }
            else
            {
                return Math.Cos(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic tan(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Tan();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Tan(c);
            }
            else
            {
                return Math.Tan(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic asin(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Asin();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Asin(c);
            }
            else
            {
                return Math.Asin(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic acos(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Acos();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Acos(c);
            }
            else
            {
                return Math.Acos(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic atan(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Atan();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Atan(c);
            }
            else
            {
                return Math.Atan(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic sinh(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Sinh();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Sinh(c);
            }
            else
            {
                return Math.Sinh(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic cosh(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Cosh();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Cosh(c);
            }
            else
            {
                return Math.Cosh(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic tanh(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Tanh();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Tanh(c);
            }
            else
            {
                return Math.Tanh(c);
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic asinh(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Asinh();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Asinh(c);
            }
            else
            {
                return c.Asinh().Real;
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic acosh(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Acosh(c);
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Acosh(c);
            }
            else
            {
                return c.Acosh().Real;
            }
            throw new Exception("the argument is not a number type");
        }

        static dynamic atanh(dynamic c)
        {
            Type unknown = c.GetType();
            if (unknown == typeof(Complex))
            {
                return c.Atanh();
            }
            else if (unknown == typeof(ValDer))
            {
                return ValDer.Atanh(c);
            }
            else
            {
                return c.Atanh().Real;
            }
            throw new Exception("the argument is not a number type");
        }

        #endregion

        #region O P E R A T O R S

        public static ValDer operator -(ValDer m)
        { return vdTimes(new ValDer(-1), m); }

        public static ValDer operator +(ValDer m1, ValDer m2)
        { return vdAdd(m1, m2); }

        public static ValDer operator +(double n, ValDer m)
        { return vdAdd(new ValDer(n), m); }

        public static ValDer operator +(ValDer m, double n)
        { return vdAdd(new ValDer(n), m); }

        public static ValDer operator -(ValDer m1, ValDer m2)
        { return vdAdd(m1, -m2); }

        public static ValDer operator -(double n, ValDer m)
        { return vdAdd(new ValDer(n), -m); }

        public static ValDer operator -(ValDer m, double n)
        { return vdAdd(new ValDer(-n), m); }

        public static ValDer operator *(ValDer m1, ValDer m2)
        { return vdTimes(m1, m2); }

        public static ValDer operator *(double n, ValDer m)
        { return vdTimes(new ValDer(n), m); }

        public static ValDer operator *(ValDer m, double n)
        { return vdTimes(new ValDer(n), m); }

        public static ValDer operator /(ValDer m1, ValDer m2)
        { return vdDivide(m1, m2); }

        public static ValDer operator /(double n, ValDer m)
        { return vdDivide(new ValDer(n), m); }

        public static ValDer operator /(ValDer m, double n)
        { return vdDivide(m, new ValDer(n)); }

        public static ValDer operator ^(ValDer n, ValDer m)
        { return vdPow(n, m); }

        #endregion

    }
}
