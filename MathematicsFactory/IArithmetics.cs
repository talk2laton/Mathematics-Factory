using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    public interface IArithmetics<T>
    {
        dynamic Add(Complex other);
        dynamic Add(Doub other);
        dynamic Add(Int other);
        dynamic Add(T other);

        dynamic Mult(Complex other);
        dynamic Mult(Doub other);
        dynamic Mult(Int other);
        dynamic Mult(T other);

        dynamic Sub(Complex other);
        dynamic Sub(Doub other);
        dynamic Sub(Int other);
        dynamic Sub(T other);

        dynamic Div(Complex other);
        dynamic Div(Doub other);
        dynamic Div(Int other);
        dynamic Div(T other);

        dynamic Pow(Complex other);
        dynamic Pow(Doub other);
        dynamic Pow(Int other);
        dynamic Pow(T other);
    }
}
