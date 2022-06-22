using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    public interface ITranscedentals<T>
    {
        double Abs();
        double Arg();
        T Sign();
        T Sin();
        T Cos();
        T Tan();
        T Sinh();
        T Cosh();
        T Tanh();
        Complex Asin();
        Complex Acos();
        T Atan();
        T Asinh();
        Complex Acosh();
        Complex Atanh();
        Complex Sqrt();
        T Sqrt(T other);
        Complex Log();
        T Conj();
        T Exp();
        T Sqr();
        T One();
        T Zero();
        T TINYNO();
        T BIGNO();
        T Rand();
        T Round(int decpts);
        bool IsReal();
        string ToString4Array();
        dynamic REAL();
        dynamic IMAG();
    }
}
