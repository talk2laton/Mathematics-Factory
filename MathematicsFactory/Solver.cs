using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    class Solver
    {
        /// <summary>
        /// Returns the gradient of a scaler multivariate function with respect to its variable
        /// </summary>
        /// <param name="fun">Multivariate Scaler function</param>
        /// <param name="x">Column Matrix of varaible valuesat which the gradiate is to be computed</param>
        /// <returns>Column Matrix of partial derivatives of the function with respect to the varaibles</returns>
        public static Matrix<Doub> Grad(Func<Matrix<Doub>, Doub> fun, Matrix<Doub> x)
        {
            int I = x.Numel;
            Doub del = 1e-6, deli, fm = 0, fp = 0;
            Matrix<Doub> ans = new Doub[I];
            for (int i = 0; i < I; i++)
            {
                deli = del * (1 + x[i].Abs());
                Matrix<Doub> Xm = x.Duplicate();
                Matrix<Doub> Xp = x.Duplicate();
                Xm[i] = Xm[i] - deli; fm = fun(Xm);
                Xp[i] = Xp[i] + deli; fp = fun(Xp);
                ans[i] = (fp - fm) / (2 * deli);
            }
            return ans;
        }

        /// <summary>
        /// Returns the gradient of a scaler multivariate function with respect to its variable
        /// </summary>
        /// <param name="fun">Multivariate Scaler function</param>
        /// <param name="x">Column Matrix of varaible valuesat which the gradiate is to be computed</param>
        /// <returns>Column Matrix of partial derivatives of the function with respect to the varaibles</returns>
        public static Matrix<Doub> Jacobian(Func<Matrix<Doub>, Matrix<Doub>> fun, Matrix<Doub> x, Matrix<Doub> f)
        {
            int I = f.Numel, J = x.Numel;
            Doub del = 1e-6, delj;
            Matrix<Doub> ans = new Doub[I, J], fm = new Doub[I], fp = new Doub[I];
            for (int j = 0; j < J; j++)
            {
                delj = del * (1 + x[j].Abs());
                Matrix<Doub> Xp = x.Duplicate();
                Xp[j] = Xp[j] + delj; fp = fun(Xp);
                //ans[j, "All"] = (fp - f) / (delj);
            }
            return ans;
        }
    }
}
