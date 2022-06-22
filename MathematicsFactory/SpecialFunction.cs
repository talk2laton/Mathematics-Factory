using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    /// <summary>
    /// Handles all computaton invloving special functions
    /// </summary>
    public class SpecialFunctions
    {
        static double[] a = new double[101];
        static double[] Fvals = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
                                 362880, 3628800, 39916800, 479001600,
                                 6227020800, 87178291200, 1307674368000,
                                 20922789888000, 355687428096000, 6402373705728000};

        static double Fs = 0, Fc = 0, Fx = 0;
        /// <summary>
        /// Returns the value ln[Γ(xx)] for xx>0.
        /// </summary>
        /// <param name="xx">value of argument for which LnGamma is to be computed</param>
        /// <returns>ln[Γ(xx)]</returns>
        public static double LnGamma(double xx)
        {
            double x, y, tmp, ser;
            double[] cof = {76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
            int j;
            y = xx;
            x = xx;
            tmp = x + 5.5;
            tmp -= (x + 0.5) * Math.Log(tmp);
            ser = 1.000000000190015;
            for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
            return -tmp + Math.Log(2.5066282746310005 * ser / x);
        }

        /// <summary>
        /// BesselJ computes the bessel function of the first kind for order n and argument x;
        /// </summary>
        /// <param name="n">Order of Bessel Function</param>
        /// <param name="x">value of argument for which Bessel Function is to be computed</param>
        /// <returns>J_n(x)</returns>
        public static double BesselJ(int n, double x)
        {
            double ACC = 40.0, BIGNO = 1.0e10, BIGNI = 1.0e-10;
            int j, m;
            bool jsum;
            double ax, bj, bjm, bjp, sum, tox, ans;
            if (n == 0)
            {
                ans = BesselJ0(x);
            }
            else if (n == 1)
            {
                ans = BesselJ1(x);
            }
            else
            {

                ax = Math.Abs(x);
                if (ax == 0.0)
                {
                    return 0.0;
                }
                else if (ax > (double)n) // Upwards recurrence fromJ0andJ1.
                {
                    tox = 2.0 / ax; bjm = BesselJ0(ax); bj = BesselJ1(ax);
                    for (j = 1; j < n; j++)
                    {
                        bjp = j * tox * bj - bjm;
                        bjm = bj;
                        bj = bjp;
                    }
                    ans = bj;
                }
                else //Downwards recurrence from an evenmhere com-puted. 
                {
                    tox = 2.0 / ax;
                    m = 2 * ((n + (int)Math.Sqrt(ACC * n)) / 2);
                    jsum = false; // jsum will alternate between 0 and 1;
                    bjp = ans = sum = 0.0;
                    bj = 1.0;
                    for (j = m; j > 0; j--) // The downward recurrence.
                    {
                        bjm = j * tox * bj - bjp;
                        bjp = bj;
                        bj = bjm;
                        if (Math.Abs(bj) > BIGNO) // Renormalize to prevent overflows.
                        {
                            bj *= BIGNI;
                            bjp *= BIGNI;
                            ans *= BIGNI;
                            sum *= BIGNI;
                        }
                        if (jsum) //Accumulate the sum.
                        {
                            sum += bj;
                        }
                        jsum = !jsum; //Change 0 to 1 or vice versa.
                        if (j == n)  //Save the unnormalized answer.
                        {
                            ans = bjp;
                        }
                    }
                    sum = 2.0 * sum - bj; // Compute (5.5.16)
                    ans /= sum; // and use it to normalize the answer.
                }
                if ((x < 0.0) && (n > 0))
                {
                    ans = -ans;
                }
            }
            return ans;
        }

        /// <summary>
        /// BesselY computes the bessel function of the second kind for order n and argument x;
        /// </summary>
        /// <param name="n">Order of Bessel Function</param>
        /// <param name="x">value of argument for which Bessel Function is to be computed</param>
        /// <returns>Y_n(x)</returns>
        public static double BesselY(int n, double x)
        {
            int j;
            double by, bym, byp, tox, ans;
            if (n == 0)
            {
                ans = BesselY0(x);
            }
            else if (n == 1)
            {
                ans = BesselY1(x);
            }
            else
            {

                tox = 2.0 / x;
                by = BesselY1(x); //Starting values for the recurrence.
                bym = BesselY0(x);
                for (j = 1; j < n; j++) // Recurrence (6.5.7).
                {
                    byp = j * tox * by - bym;
                    bym = by;
                    by = byp;
                }
                ans = by;
            }
            return ans;
        }
        

        /// <summary>
        /// BesselJ computes the bessel function of the first kind for order n and argument x;
        /// </summary>
        /// <param name="n">Order of Bessel Function</param>
        /// <param name="x">value of argument for which Bessel Function is to be computed</param>
        /// <returns>J_n(x)</returns>
        public static Complex BesselJ(int n, Complex x)
        {
            Complex ACC = 40.0, BIGNO = 1.0e10, BIGNI = 1.0e-10;
            int j, m;
            bool jsum;
            Complex ax, bj, bjm, bjp, sum, tox, ans;
            if (n == 0)
            {
                ans = BesselJ0(x);
            }
            else if (n == 1)
            {
                ans = BesselJ1(x);
            }
            else
            {

                ax = x.Abs();
                if (x.Abs() == 0.0)
                {
                    return 0.0;
                }
                else if (x.Abs() > n) // Upwards recurrence fromJ0andJ1.
                {
                    tox = 2.0 / ax; bjm = BesselJ0(ax); bj = BesselJ1(ax);
                    for (j = 1; j < n; j++)
                    {
                        bjp = j * tox * bj - bjm;
                        bjm = bj;
                        bj = bjp;
                    }
                    ans = bj;
                }
                else //Downwards recurrence from an evenmhere com-puted. 
                {
                    tox = 2.0 / ax;
                    m = 2 * ((n + Sqrt(ACC * n)) / 2);
                    jsum = false; // jsum will alternate between 0 and 1;
                    bjp = ans = sum = 0.0;
                    bj = 1.0;
                    for (j = m; j > 0; j--) // The downward recurrence.
                    {
                        bjm = j * tox * bj - bjp;
                        bjp = bj;
                        bj = bjm;
                        if (bj.Abs() > BIGNO.Abs()) // Renormalize to prevent overflows.
                        {
                            bj *= BIGNI;
                            bjp *= BIGNI;
                            ans *= BIGNI;
                            sum *= BIGNI;
                        }
                        if (jsum) //Accumulate the sum.
                        {
                            sum += bj;
                        }
                        jsum = !jsum; //Change 0 to 1 or vice versa.
                        if (j == n)  //Save the unnormalized answer.
                        {
                            ans = bjp;
                        }
                    }
                    sum = 2.0 * sum - bj; // Compute (5.5.16)
                    ans /= sum; // and use it to normalize the answer.
                }
                if ((x.Real < 0.0) && (n > 0))
                {
                    ans = -ans;
                }
            }
            return ans;
        }

        /// <summary>
        /// BesselY computes the bessel function of the second kind for order n and argument x;
        /// </summary>
        /// <param name="n">Order of Bessel Function</param>
        /// <param name="x">value of argument for which Bessel Function is to be computed</param>
        /// <returns>Y_n(x)</returns>
        public static Complex BesselY(int n, Complex x)
        {
            int j;
            Complex by, bym, byp, tox, ans;
            if (n == 0)
            {
                ans = BesselY0(x);
            }
            else if (n == 1)
            {
                ans = BesselY1(x);
            }
            else
            {

                tox = 2.0 / x;
                by = BesselY1(x); //Starting values for the recurrence.
                bym = BesselY0(x);
                for (j = 1; j < n; j++) // Recurrence (6.5.7).
                {
                    byp = j * tox * by - bym;
                    bym = by;
                    by = byp;
                }
                ans = by;
            }
            return ans;
        }

        /// <summary>
        /// BesselI computes the modified bessel function of the first kind for order n and argument x;
        /// </summary>
        /// <param name="n">Order of Bessel Function</param>
        /// <param name="x">value of argument for which Bessel Function is to be computed</param>
        /// <returns>I_n(x)</returns>
        public static Complex BesselI(int n, Complex x)
        {
            double ACC = 40.0, BIGNO = 1.0e10, BIGNI = 1.0e-10;
            int j;
            Complex bi, bim, bip, tox, ans;
            if (n == 0)
            {
                ans = BesselI0(x);
            }
            else if (n == 1)
            {
                ans = BesselI1(x);
            }
            else
            {
                if (Abs(x) == 0.0)
                {
                    ans = 0.0;
                }
                else //Downwards recurrence from an evenmhere com-puted. 
                {
                    tox = 2.0 / Abs(x);
                    bip = ans = 0.0;
                    bi = 1.0;
                    for (j = 2 * (n + (int)Math.Sqrt(ACC * n)); j > 0; j--) // Downward recurrence from even
                    {
                        bim = bip + j * tox * bi;
                        bip = bi;
                        bi = bim;
                        if (Abs(bi) > BIGNO) // Renormalize to prevent overflows.
                        {
                            ans *= BIGNI;
                            bi *= BIGNI;
                            bip *= BIGNI;
                        }
                        if (j == n) ans = bip;
                    }
                    ans *= BesselI0(x) / bi;
                    if ((x.Real < 0.0) && (n > 0))
                    {
                        ans = -ans;
                    }
                }
            }
            return ans;
        }

        /// <summary>
        /// BesselK computes the modified bessel function of the second kind for order n and argument x;
        /// </summary>
        /// <param name="n">Order of Bessel Function</param>
        /// <param name="x">value of argument for which Bessel Function is to be computed</param>
        /// <returns>K_n(x)</returns>
        public static Complex BesselK(int n, Complex x)
        {
            int k;
            Complex bk, bkm, bkp, tox, ans;
            if (n == 0)
            {
                ans = BesselK0(x);
            }
            else if (n == 1)
            {
                ans = BesselK1(x);
            }
            else
            {

                tox = 2.0 / x;
                bkm = BesselK0(x); //Upward recurrence for allx...
                bk = BesselK1(x);
                for (k = 1; k < n; k++) // ...and here it is.
                {
                    bkp = bkm + k * tox * bk;
                    bkm = bk;
                    bk = bkp;
                }
                ans = bk;
            }
            return ans;
        }


        /// <summary>
        /// BesselI computes the modified bessel function of the first kind for order n and argument x;
        /// </summary>
        /// <param name="n">Order of Bessel Function</param>
        /// <param name="x">value of argument for which Bessel Function is to be computed</param>
        /// <returns>I_n(x)</returns>
        public static double BesselI(int n, double x)
        {
            double ACC = 40.0, BIGNO = 1.0e10, BIGNI = 1.0e-10;
            int j;
            double bi, bim, bip, tox, ans;
            if (n == 0)
            {
                ans = BesselI0(x);
            }
            else if (n == 1)
            {
                ans = BesselI1(x);
            }
            else
            {
                if (x == 0.0)
                {
                    ans = 0.0;
                }
                else //Downwards recurrence from an evenmhere com-puted. 
                {
                    tox = 2.0 / Math.Abs(x);
                    bip = ans = 0.0;
                    bi = 1.0;
                    for (j = 2 * (n + (int)Math.Sqrt(ACC * n)); j > 0; j--) // Downward recurrence from even
                    {
                        bim = bip + j * tox * bi;
                        bip = bi;
                        bi = bim;
                        if (Math.Abs(bi) > BIGNO) // Renormalize to prevent overflows.
                        {
                            ans *= BIGNI;
                            bi *= BIGNI;
                            bip *= BIGNI;
                        }
                        if (j == n) ans = bip;
                    }
                    ans *= BesselI0(x) / bi;
                    if ((x < 0.0) && (n > 0))
                    {
                        ans = -ans;
                    }
                }
            }
            return ans;
        }

        /// <summary>
        /// BesselK computes the modified bessel function of the second kind for order n and argument x;
        /// </summary>
        /// <param name="n">Order of Bessel Function</param>
        /// <param name="x">value of argument for which Bessel Function is to be computed</param>
        /// <returns>K_n(x)</returns>
        public static double BesselK(int n, double x)
        {
            int k;
            double bk, bkm, bkp, tox, ans;
            if (n == 0)
            {
                ans = BesselK0(x);
            }
            else if (n == 1)
            {
                ans = BesselK1(x);
            }
            else
            {

                tox = 2.0 / x;
                bkm = BesselK0(x); //Upward recurrence for allx...
                bk = BesselK1(x);
                for (k = 1; k < n; k++) // ...and here it is.
                {
                    bkp = bkm + k * tox * bk;
                    bkm = bk;
                    bk = bkp;
                }
                ans = bk;
            }
            return ans;
        }

        public static double BesselJ(double n, double x)
        {
            int nint = (int)n;
            double ans = 0;
            if (nint - n == 0)
            {
                ans = BesselJ(nint, x);
            }
            else
            {
                double ry;
                double rjp;
                double ryp;
                besseljy(x, n, out ans, out ry, out rjp, out ryp);
            }
            return ans;
        }

        public static double BesselY(double n, double x)
        {
            int nint = (int)n;
            double ans = 0;
            if (nint - n == 0)
            {
                ans = BesselY(nint, x);
            }
            else
            {
                double rj;
                double rjp;
                double ryp;
                besseljy(x, n, out rj, out ans, out rjp, out ryp);
            }
            return ans;
        }

        /// <summary>
        /// Returns the value Γ(x).
        /// </summary>
        /// <param name="x">value of argument for which Gamma is to be computed</param>
        /// <returns>Γ(x)</returns>
        public static double Gamma(double x)
        {
            return ((int)x - x == 0) ? Factorial((int)x - 1) : Math.Exp(LnGamma(x));
        }

        /// <summary>
        /// Returns the value Γ_p(a,x).
        /// </summary>
        /// <param name="a">value of argument for which Upper Incomplete Gamma is to be computed</param>
        /// <param name="x">start of integration</param>
        /// <returns>Γ_p(a,x)</returns>
        public static double GammaP(double a, double x)
        {
            //Returns the incomplete gamma functionP(a, x).
            if (x < 0.0 || a <= 0.0) throw new Exception("Invalid arguments in routine gammp");
            if (x < (a + 1.0))
            { //Use the series representation.
                double[] serln = gser(a, x);
                return serln[0];
            }
            else
            {
                //Use the continued fraction representation
                double[] cfln = gcf(a, x);
                return 1.0 - cfln[0]; // and take its complement.
            }
        }

        /// <summary>
        /// Returns the value Γ_q(a,x).
        /// </summary>
        /// <param name="a">value of argument for which Lower Incomplete Gamma is to be computed</param>
        /// <param name="x">end of integration</param>
        /// <returns>Γ_p(a,x)</returns>
        public static double GammaQ(double a, double x)
        {
            //Returns the incomplete gamma functionQ(a, x).
            if (x < 0.0 || a <= 0.0) throw new Exception("Invalid arguments in routine gammq");
            if (x < (a + 1.0))
            { //Use the series representation
                double[] serln = gser(a, x);
                return 1 - serln[0];
            }
            else
            { //Use the continued fraction representation.
                double[] cfln = gcf(a, x);
                return cfln[0]; // and take its complement.
            }
        }

        /// <summary>
        /// Returns the value Γ_p^-1(p, a).
        /// </summary>
        /// <param name="a">value of argument for which Upper Incomplete Gamma is to be computed</param>
        /// <param name="x">start of integration</param>
        /// <returns>Γ_p^-1(p, a)</returns>
        public static double InvGammaP(double p, double a)
        {
            //Returns x such that Γ_p(a,x) = p for an argument p between 0 and 1.
            int j;
            double x, err, t, u, pp, lna1 = 0, gln, afac = 0, a1 = a - 1;
            const double EPS = 1e-8; //Accuracy is the square of EPS.
            gln = LnGamma(a);
            if (a <= 0.0) throw new SpecialFunctionsException("a must be pos in InvGammaP");
            if (p >= 1.0) return Max(100.0, a + 100.0 * Sqrt(a));
            if (p <= 0.0) return 0.0;
            if (a > 1.0)
            {
                //Initial guess based on reference [1].
                lna1 = Ln(a1);
                afac = Exp(a1 * (lna1 - 1.0) - gln);
                pp = (p < 0.5) ? p : 1.0 - p;
                t = Sqrt(-2.0 * Ln(pp));
                x = (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t;
                if (p < 0.5) x = -x;
                x = Max(1e-3, a * Pow(1.0 - 1.0 / (9.0 * a) - x / (3.0 * Sqrt(a)), 3));
            }
            else
            {
                //Initial guess based on equations (6.2.8) and (6.2.9). 
                t = 1.0 - a * (0.253 + a * 0.12);
                if (p < t) x = Pow(p / t, 1.0 / a);
                else x = 1.0 - Ln(1.0 - (p - t) / (1.0 - t));
            }
            for (j = 0; j < 12; j++)
            {
                if (x <= 0.0) return 0.0; //x too small to compute accurately.
                err = GammaP(a, x) - p;
                if (a > 1.0) t = afac * Exp(-(x - a1) + a1 * (Ln(x) - lna1));
                else t = Exp(-x + a1 * Ln(x) - gln);
                u = err / t;
                x -= (t = u / (1.0 - 0.5 * Min(1.0, u * ((a - 1.0) / x - 1)))); //Halley’s method.
                if (x <= 0.0) x = 0.5 * (x + t);// Halve old value if x tries to go negative.
                if (Abs(t) < EPS * x) break;
            }
            return x;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="t"></param>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double HeavisideTheta(double t, double a = 0)
        {
            return t >= a ? 1 : 0;
        }

        /// <summary>
        /// Error Function
        /// </summary>
        /// <param name="x">value of argument for which the error function is to be computed</param>
        /// <returns>Erf[x]</returns>
        public static double Erf(double x)
        { //Returns the error function erf(x).
            return x < 0.0 ? -GammaP(0.5, x * x) : GammaP(0.5, x * x);
        }

        /// <summary>
        /// Error Function
        /// </summary>
        /// <param name="x">value of argument for which the error function is to be computed</param>
        /// <returns>Erf[x]</returns>
        public static double InvErf(double x)
        { //Returns a value whose error function is z.
            double a = 0.147, the_sign_of_x = Math.Sign(x), z;
            if (x == 0)
            {
                return 0;
            }
            else
            {
                var ln_1minus_x_sqrd = Math.Log(1 - x * x);
                var ln_1minusxx_by_a = ln_1minus_x_sqrd / a;
                var ln_1minusxx_by_2 = ln_1minus_x_sqrd / 2;
                var ln_etc_by2_plus2 = ln_1minusxx_by_2 + (2 / (Math.PI * a));
                var first_sqrt = Math.Sqrt((ln_etc_by2_plus2 * ln_etc_by2_plus2) - ln_1minusxx_by_a);
                var second_sqrt = Math.Sqrt(first_sqrt - ln_etc_by2_plus2);
                return second_sqrt * the_sign_of_x;
            }
        }

        public static double Laguerre(int n, double x)
        {
            double L0 = 1, L1 = -x + 1, L2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return L0;
            else if (n == 1) return L1;
            else
            {
                while (i < n)
                {
                    L2 = ((2.0 * i + 1.0 - x) * L1 - i * L0) / (i + 1);
                    L0 = L1; L1 = L2; i++;
                }
                return L2;
            }
        }

        public static double Hermite(int n, double x)
        {
            double H0 = 1, H1 = 2 * x, H2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return H0;
            else if (n == 1) return H1;
            else
            {
                while (i < n)
                {
                    H2 = 2 * x * H1 - 2 * i * H0;
                    H0 = H1; H1 = H2; i++;
                }
                return H2;
            }
        }

        public static double ChebyshevT(int n, double x)
        {
            double T0 = 1, T1 = x, T2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return T0;
            else if (n == 1) return T1;
            else
            {
                while (i < n)
                {
                    T2 = 2 * x * T1 - T0;
                    T0 = T1; T1 = T2; i++;
                }
                return T2;
            }
        }

        public static double ChebyshevU(int n, double x)
        {
            double U0 = 1, U1 = 2 * x, U2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return U0;
            else if (n == 1) return U1;
            else
            {
                while (i < n)
                {
                    U2 = 2 * x * U1 - U0;
                    U0 = U1; U1 = U2; i++;
                }
                return U2;
            }
        }

        public static double Legendre(int n, double x)
        {
            double P0 = 1, P1 = x, P2 = 0;
            int i = 1;
            if (n < 0) return -1;
            if (n == 0) return P0;
            else if (n == 1) return P1;
            else
            {
                while (i < n)
                {
                    P2 = 2 * i * P1 - P0 - (x * P1 - P0) / (i + 1);
                    P0 = P1; P1 = P2; i++;
                }
                return P2;
            }
        }

        /// <summary>
        /// Complemetary Error Function
        /// </summary>
        /// <param name="x">value of argument for which the complemetary error function is to be computed</param>
        /// <returns>Erfc[x]</returns>
        public static double Erfc(double x)
        { //Returns the complementary error function erfc(x).
            return x < 0.0 ? 1.0 + GammaP(0.5, x * x) : GammaQ(0.5, x * x);
        }

        /// <summary>
        /// Zeta Function
        /// </summary>
        /// <param name="x">value of argument for which the zeta function is to be computed</param>
        /// <returns>Z[x]</returns>
        public static double Zeta(double x)
        {
            double s = 1, sum = 1; int i = 1;
            if (x < 1)
            {
                return 1 - Zeta(2 - x);
            }
            else
            {
                while (s > 1e-20)
                {
                    i++;
                    s = Math.Pow(i, -x);
                    sum += s;
                }
                return sum;
            }
        }

        /// <summary>
        /// Beta Function
        /// </summary>
        /// <param name="a">value of the first argument of the beta function</param>
        /// <param name="b">value of the second argument of the beta function</param>
        /// <returns>B[a,b]</returns>
        public static double Beta(double a, double b)
        {
            return Gamma(a) * Gamma(b) / Gamma(a + b);
        }

        /// <summary>
        /// Incomplete Beta Function
        /// </summary>
        /// <param name="a">value of the first argument of the beta function</param>
        /// <param name="b">value of the second argument of the beta function</param>
        /// <returns>B[a, b, x]</returns>
        public static double IncBeta(double a, double b, double x)
        {
            double bt; int SWITCH = 3000;
            if (a <= 0.0 || b <= 0.0) throw new SpecialFunctionsException("Bad a or b in routine betai");
            if (x < 0.0 || x > 1.0) throw new SpecialFunctionsException("Bad x in routine betai");
            if (x == 0.0 || x == 1.0) return x;
            if (a > SWITCH && b > SWITCH) return betaiapprox(a, b, x);
            bt = Exp(LnGamma(a + b) - LnGamma(a) - LnGamma(b) + a * Ln(x) + b * Ln(1.0 - x));
            if (x < (a + 1.0) / (a + b + 2.0)) return bt * betacf(a, b, x) / a;
            else return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
        }

        /// <summary>
        /// Inverse Incomplete Beta Function
        /// </summary>
        /// <param name="a">value of the first argument of the beta function</param>
        /// <param name="b">value of the second argument of the beta function</param>
        /// <returns>B^-1[p, a, b]</returns>
        public static double InvIncBeta(double p, double a, double b)
        {
            return invbetai(p, a, b);
        }

        /// <summary>
        /// x!
        /// </summary>
        /// <param name="x">value of argument for which factorial is to be computed</param>
        /// <returns>x!</returns>
        public static double Factorial(int x)
        {
            return x < 19 ? Fvals[x] : Exp(LnGamma(x + 1));
        }

        /// <summary>
        /// Exponential Integral(It is defined as a  definite integral of the ratio between an exponential function and powers of its argument.)
        /// </summary>
        /// <param name="n">power of the exponential argument</param>
        /// <param name="x">lower limit of integration</param>
        /// <returns>ExpInt_n[x]</returns>
        public static double ExpInt(int n, double x)
        {
            int MAXIT = 100; //Maximum allowed number of iterations.
            double EULER = 0.5772156649;// Euler’s constantγ.
            double FPMIN = 1.0e-30; // Close to smallest representable floating-point number.
            double EPS = 1.0e-7; // Desired relative error, not smaller than the machine pre-cision.
            int i, ii, nm1;
            double a, b, c, d, del, fact, h, psi, ans;
            nm1 = n - 1;
            if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
                throw new Exception("bad arguments in expint");
            else
            {
                if (n == 0) ans = Math.Exp(-x) / x; //Special case.
                else
                {
                    if (x == 0.0) ans = 1.0 / nm1; // Another special case.
                    else
                    {
                        if (x > 1.0)
                        { //Lentz’s algorithm (§5.2).
                            b = x + n;
                            c = 1.0 / FPMIN;
                            d = 1.0 / b;
                            h = d;
                            for (i = 1; i <= MAXIT; i++)
                            {
                                a = -i * (nm1 + i);
                                b += 2.0;
                                d = 1.0 / (a * d + b); //Denominators cannot be zero.
                                c = b + a / c;
                                del = c * d;
                                h *= del;
                                if (Abs(del - 1.0) < EPS)
                                {
                                    ans = h * Math.Exp(-x);
                                    return ans;
                                }
                            }
                            throw new Exception("continued fraction failed in expint");
                        }
                        else
                        { //Evaluate series.
                            ans = (nm1 != 0 ? 1.0 / nm1 : -Math.Log(x) - EULER); //Setfirstterm.
                            fact = 1.0;
                            for (i = 1; i <= MAXIT; i++)
                            {
                                fact *= -x / i;
                                if (i != nm1) del = -fact / (i - nm1);
                                else
                                {
                                    psi = -EULER; //Computeψ(n).
                                    for (ii = 1; ii <= nm1; ii++) psi += 1.0 / ii;
                                    del = fact * (-Math.Log(x) + psi);
                                }
                                ans += del;
                                if (Abs(del) < Abs(ans) * EPS) return ans;
                            }
                            throw new Exception("series failed in expint");
                        }
                    }
                }
            }
            return ans;
        }

        /// <summary>
        /// Exponential Integral(It is defined as a  definite integral of the ratio between an exponential function and its argument.)
        /// </summary>
        /// <param name="x">Negative of lower limit of integration</param>
        /// <returns>EI[x] = -ExpInt[-x]</returns>
        public static double EI(double x)
        {
            int MAXIT = 100; //Maximum allowed number of iterations.
            double EULER = 0.5772156649;// Euler’s constantγ.
            double FPMIN = 1.0e-30; // Close to smallest representable floating-point number.
            double EPS = 1.0e-8; // Desired relative error, not smaller than the machine pre-cision.
            int k;
            double fact, prev, sum, term;
            if (x <= 0.0) throw new Exception("Bad argument in ei");
            if (x < FPMIN) return Math.Log(x) + EULER; //Special case: avoid failure of convergence test because of underflow. 
            if (x <= -Math.Log(EPS))
            {
                sum = 0.0; // Use power series.
                fact = 1.0;
                for (k = 1; k <= MAXIT; k++)
                {
                    fact *= x / k;
                    term = fact / k;
                    sum += term;
                    if (term < EPS * sum) break;
                }
                if (k > MAXIT) throw new Exception("Series failed in ei");
                return sum + Math.Log(x) + EULER;
            }
            else
            {// Use asymptotic series.
                sum = 0.0; //Start with second term.
                term = 1.0;
                for (k = 1; k <= MAXIT; k++)
                {
                    prev = term;
                    term *= k / x;
                    if (term < EPS) break;
                    //Since final sum is greater than one,termitself approximates the relative error.
                    if (term < prev) sum += term; //Still converging: add new term.
                    else
                    {
                        sum -= prev; //Diverging: subtract previous term and exit. 
                        break;
                    }
                }
                return Math.Exp(x) * (1.0 + sum) / x;
            }
        }

        /// <summary>
        /// Sine Integral(It is defined as a definite integral of the ratio between an sine function and its argument.)
        /// </summary>
        /// <param name="x">Lower limit of integration</param>
        /// <returns>SI[x]</returns>
        public static double SI(double x)
        {
            return cisi(x)[0];
        }

        /// <summary>
        /// Cosine Integral(It is defined as a definite integral of the ratio between an cosine function and its argument.)
        /// </summary>
        /// <param name="x">Lower limit of integration</param>
        /// <returns>CI[x]</returns>
        public static double CI(double x)
        {
            return cisi(x)[1];
        }

        /// <summary>
        /// Fresnel Integral (It is defined as a definite integral of sine of square of its argument.)
        /// </summary>
        /// <param name="x">Upper limit of the integration</param>
        /// <returns>S[x]</returns>
        public static double FresnelS(double x)
        {
            return frenel(x)[0];
        }

        /// <summary>
        /// Fresnel Integral (It is defined as a definite integral of cosine of square of its argument.)
        /// </summary>
        /// <param name="x">Upper limit of the integration</param>
        /// <returns>C[x]</returns>
        public static double FresnelC(double x)
        {
            return frenel(x)[1];
        }

        public static double EllipticF(double phi, double k)
        {
            double s = Sin(phi);
            return s * rf(SQR(Cos(phi)), (1.0 - s * k) * (1.0 + s * k), 1.0);
        }

        public static double EllipticK(double k)
        {
            double phi = Math.PI / 2;
            double s = Sin(phi);
            return s * rf(SQR(Cos(phi)), (1.0 - s * k) * (1.0 + s * k), 1.0);
        }

        public static double EllipticE(double phi, double k)
        {
            double s = Sin(phi), cc = SQR(Cos(phi)), q = (1 - s * k) * (1 + s * k);
            return s * (rf(cc, q, 1.0) - (SQR(s * k)) * rd(cc, q, 1.0) / 3.0);
        }

        public static double EllipticE(double k)
        {
            double phi = Math.PI / 2;
            double s = Sin(phi), cc = SQR(Cos(phi)), q = (1 - s * k) * (1 + s * k);
            return s * (rf(cc, q, 1.0) - (SQR(s * k)) * rd(cc, q, 1.0) / 3.0);
        }

        public static double EllipticPi(double phi, double n, double k)
        {
            double cc, enss, q, s;
            s = Sin(phi);
            enss = n * s * s;
            cc = SQR(Cos(phi));
            q = (1.0 - s * k) * (1.0 + s * k);
            return s * (rf(cc, q, 1.0) - enss * rj(cc, q, 1.0, 1.0 + enss) / 3.0);
        }

        public static double EllipticPi(double n, double k)
        {
            double cc, enss, q, s, phi = Math.PI / 2;
            s = Sin(phi);
            enss = n * s * s;
            cc = SQR(Cos(phi));
            q = (1.0 - s * k) * (1.0 + s * k);
            return s * (rf(cc, q, 1.0) - enss * rj(cc, q, 1.0, 1.0 + enss) / 3.0);
        }

        #region Private Functions

        static double rf(double x, double y, double z)
        {
            //Computes Carlson’s elliptic integral of the first kind, RF(x,y,z)
            // x, y,and z must be non-negative, and at most one can be zero.
            double ERRTOL = 0.0025, THIRD = 1.0 / 3.0, C1 = 1.0 / 24.0, C2 = 0.1, C3 = 3.0 / 44.0, C4 = 1.0 / 14.0,
                TINY = 5.0 * 2.22507e-308, BIG = 0.2 * 1.79769e+308, alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty,
                sqrtz, xt, yt, zt;
            if (Min(Min(x, y), z) < 0.0 || Min(Min(x + y, x + z), y + z) < TINY || Max(Max(x, y), z) > BIG)
                throw new SpecialFunctionsException("invalid arguments in rf");
            xt = x; yt = y; zt = z;
            do
            {
                sqrtx = Sqrt(xt);
                sqrty = Sqrt(yt);
                sqrtz = Sqrt(zt);
                alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
                xt = 0.25 * (xt + alamb);
                yt = 0.25 * (yt + alamb);
                zt = 0.25 * (zt + alamb);
                ave = THIRD * (xt + yt + zt);
                delx = (ave - xt) / ave;
                dely = (ave - yt) / ave;
                delz = (ave - zt) / ave;
            } while (Max(Max(Abs(delx), Abs(dely)), Abs(delz)) > ERRTOL);
            e2 = delx * dely - delz * delz;
            e3 = delx * dely * delz;
            return (1.0 + (C1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / Sqrt(ave);
        }

        static double rd(double x, double y, double z)
        {
            //Computes Carlson’s elliptic integral of the second kind, RD(x,y,z)
            // x and y must be non-negative, and at most one can be zero. z must be positive.
            double ERRTOL = 0.0015, C1 = 3.0 / 14.0, C2 = 1.0 / 6.0, C3 = 9.0 / 22.0, C4 = 3.0 / 26.0,
                C5 = 0.25 * C3, C6 = 1.5 * C4, TINY = 2.0 * Pow(1.79769e+308, -2.0 / 3), BIG = 0.2 * ERRTOL * Pow(2.22507e-308, -2.0 / 3),
                     alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;
            if (Min(x, y) < 0.0 || Min(x + y, z) < TINY || Max(Max(x, y), z) > BIG)
                throw new SpecialFunctionsException("invalid arguments in rd");
            xt = x; yt = y; zt = z; sum = 0.0; fac = 1.0;
            do
            {
                sqrtx = Sqrt(xt);
                sqrty = Sqrt(yt);
                sqrtz = Sqrt(zt);
                alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
                sum += fac / (sqrtz * (zt + alamb));
                fac = 0.25 * fac;
                xt = 0.25 * (xt + alamb);
                yt = 0.25 * (yt + alamb);
                zt = 0.25 * (zt + alamb);
                ave = 0.2 * (xt + yt + 3.0 * zt);
                delx = (ave - xt) / ave;
                dely = (ave - yt) / ave;
                delz = (ave - zt) / ave;
            } while (Max(Max(Abs(delx), Abs(dely)), Abs(delz)) > ERRTOL);
            ea = delx * dely;
            eb = delz * delz;
            ec = ea - eb;
            ed = ea - 6.0 * eb;
            ee = ed + ec + ec;
            return 3.0 * sum + fac * (1.0 + ed * (-C1 + C5 * ed - C6 * delz * ee)
            + delz * (C2 * ee + delz * (-C3 * ec + delz * C4 * ea))) / (ave * Math.Sqrt(ave));
        }

        static double rc(double x, double y)
        {
            //Computes Carlson’s degenerate elliptic integral, RC(x,y)
            // x must be non-negative, and y must be zero.
            // if y < 0 the cauchy principal is returned
            double ERRTOL = 0.0012, THIRD = 1.0 / 3.0, C1 = 0.3, C2 = 1.0 / 7.0, C3 = 0.375, C4 = 9.0 / 22.0,
                TINY = 5.0 * 2.22507e-308, BIG = 0.2 * 1.79769e+308,
                COMP1 = 2.236 / Sqrt(TINY), COMP2 = SQR(TINY * BIG) / 25.0,
                     alamb, ave, s, w, xt, yt;
            if (x < 0.0 || y == 0.0 || (x + Abs(y)) < TINY || (x + Abs(y)) > BIG ||
                (y < -COMP1 && x > 0.0 && x < COMP2))
                throw new SpecialFunctionsException("invalid arguments in rd");
            xt = x; yt = y;
            if (y > 0.0)
            {
                xt = x; yt = y; w = 1.0;
            }
            else
            {
                xt = x - y; yt = -y; w = Sqrt(x) / Sqrt(xt);
            }
            do
            {
                alamb = 2.0 * Sqrt(xt) * Sqrt(yt) + yt;
                xt = 0.25 * (xt + alamb);
                yt = 0.25 * (yt + alamb);
                ave = THIRD * (xt + yt + yt);
                s = (yt - ave) / ave;
            } while (Abs(s) > ERRTOL);
            return w * (1.0 + s * s * (C1 + s * (C2 + s * (C3 + s * C4)))) / Sqrt(ave);
        }

        static double rj(double x, double y, double z, double p)
        {
            //Computes Carlson’s elliptic integral of the second kind, RJ(x,y,z,p)
            // x, y,and z must be nonnegative, and at most one can be zero.p must be nonzero. 
            //Ifp<0, the Cauchy principal value is returned.
            double ERRTOL = 0.0015, C1 = 3.0 / 14.0, C2 = 1.0 / 6.0, C3 = 9.0 / 22.0, C4 = 3.0 / 26.0,
                C5 = 0.75 * C3, C6 = 1.5 * C4, C7 = 0.5 * C2, C8 = C3 + C3, TINY = Pow(5 * 2.22507e-308, 1.0 / 3), BIG = 0.3 * Pow(1.79769e+308, 1.0 / 3),
                     a = 0, alamb, alpha, ans, ave, b = 0, beta, delp, delx, dely, delz, ea, eb, ec, ed, ee, fac, pt, rcx = 0, rho, sqrtx, sqrty,
                     sqrtz, sum, tau, xt, yt, zt;
            if (Min(Min(x, y), z) < 0.0 || Min(Min(x + y, x + z), Min(y + z, Abs(p))) < TINY
                || Max(Max(x, y), Max(z, Abs(p))) > BIG)
                throw new SpecialFunctionsException("invalid arguments in rj");
            sum = 0.0; fac = 1.0;
            if (p > 0)
            {
                xt = x; yt = y; zt = z; pt = p;
            }
            else
            {
                xt = Min(Min(x, y), z);
                zt = Max(Max(x, y), z);
                yt = x + y + z - xt - zt;
                a = 1.0 / (yt - p);
                b = a * (zt - yt) * (yt - xt);
                pt = yt + b;
                rho = xt * zt / yt;
                tau = p * pt / yt;
                rcx = rc(rho, tau);
            }
            do
            {
                sqrtx = Sqrt(xt);
                sqrty = Sqrt(yt);
                sqrtz = Sqrt(zt);
                alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
                alpha = SQR(pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz);
                beta = pt * SQR(pt + alamb);
                sum += fac * rc(alpha, beta);
                fac = 0.25 * fac;
                xt = 0.25 * (xt + alamb);
                yt = 0.25 * (yt + alamb);
                zt = 0.25 * (zt + alamb);
                pt = 0.25 * (pt + alamb);
                ave = 0.2 * (xt + yt + zt + 2 * pt);
                delx = (ave - xt) / ave;
                dely = (ave - yt) / ave;
                delz = (ave - zt) / ave;
                delp = (ave - pt) / ave;
            } while (Max(Max(Abs(delx), Abs(dely)), Max(Abs(delz), Abs(delp))) > ERRTOL);
            ea = delx * (dely + delz) + dely * delz;
            eb = delx * dely * delz;
            ec = delp * delp;
            ed = ea - 3.0 * ec;
            ee = eb + 2.0 * delp * (ea - ec);
            ans = 3.0 * sum + fac * (1.0 + ed * (-C1 + C5 * ed - C6 * ee) + eb * (C7 + delp * (-C8 + delp * C4))
            + delp * ea * (C2 - delp * C3) - C2 * delp * ec) / (ave * Sqrt(ave));
            if (p <= 0.0) ans = a * (b * ans + 3.0 * (rcx - rf(xt, yt, zt)));
            return ans;
        }

        static double SQR(double x)
        {
            return x * x;
        }

        static double betacf(double a, double b, double x)
        {
            //Evaluates continued fraction for incomplete beta function by 
            //modified Lentz’s method (5.2). User should not call directly.
            int m, m2;
            double aa, c, d, del, h, qab, qam, qap, EPS = 1e-52, FPMIN = 2.22507e-308 / EPS;
            qab = a + b; //Theseq’s will be used in factors that occur in the coefficients (6.4.6). 
            qap = a + 1.0;
            qam = a - 1.0;
            c = 1.0; //First step of Lentz’s method.
            d = 1.0 - qab * x / qap;
            if (Abs(d) < FPMIN) d = FPMIN;
            d = 1.0 / d;
            h = d;
            for (m = 1; m < 10000; m++)
            {
                m2 = 2 * m;
                aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = 1.0 + aa * d; // One step (the even one) of the recur-rence. 
                if (Abs(d) < FPMIN) d = FPMIN;
                c = 1.0 + aa / c;
                if (Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
                d = 1.0 + aa * d; //Next step of the recurrence (the oddone). 
                if (Abs(d) < FPMIN) d = FPMIN;
                c = 1.0 + aa / c;
                if (Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                del = d * c;
                h *= del;
                if (Abs(del - 1.0) <= EPS) break; //Arewedone?
            }
            return h;
        }

        static double betaiapprox(double a, double b, double x)
        {
            //Evaluates Incomplete beta by quadrature 
            // User should not call directly.
            int j, m;
            double xu, t, sum, ans;
            double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
            double lnmu = Ln(mu), lnmuc = Ln(1.0 - mu);
            t = Sqrt(a * b / (SQR(a + b) * (a + b + 1.0)));
            if (x > a / (a + b))
            { //Set how far to integrate into the tail:
                if (x >= 1.0) return 1.0;
                xu = Min(1, Max(mu + 10 * t, x + 5.0 * t));
            }
            else
            {
                if (x <= 0.0) return 0.0;
                xu = Max(0, Min(mu - 10 * t, x - 5.0 * t));
            }
            sum = 0;
            double[] y, w;
            gauleg(18, out y, out w);
            for (j = 0; j < 18; j++)
            { //Gauss-Legendre.
                t = x + (xu - x) * y[j];
                sum += w[j] * Exp(a1 * (Ln(t) - lnmu) + b1 * (Ln(1 - t) - lnmuc));
            }
            ans = sum * (xu - x) * Exp(a1 * lnmu - LnGamma(a) + b1 * lnmuc - LnGamma(b) + LnGamma(a + b));
            return ans > 0 ? 1 - ans : -ans;
        }

        static double invbetai(double p, double a, double b)
        {
            //Inverse of incomplete beta function. Returns xsuch that I_x(a,b) = p 
            //for argument p between 0 and 1.
            const double EPS = 1e-8;
            double pp, t, u, err, x = 0, al, h, w, afac, a1 = a - 1, b1 = b - 1;
            int j;

            if (p <= 0) return 0;
            else if (p >= 1) return 1;
            else if (a >= 1 && b >= 1)
            {
                pp = (p < 0.5) ? p : 1 - p;
                t = Sqrt(-2 * Ln(pp));
                x = (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t;
                if (p < 0.5) x = -x;
                al = (SQR(x) - 3.0) / 6.0;
                h = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
                w = (x * Sqrt(al + h) / h) - (1.0 / (2.0 * b - 1) - 1.0 / (2.0 * a - 1.0)) * (al + 5.0 / 6.0 - 2.0 / (3.0 * h));
                x = a / (a + b * Exp(2 * w));
            }
            else
            {
                double lna = Ln(a / (a + b)), lnb = Ln(b / (a + b));
                t = Exp(a * lna) / a;
                u = Exp(b * lnb) / b;
                w = t + u;
                if (p < t / w) x = Pow(a * w * p, 1.0 / a);
                else x = 1.0 - Pow(b * w * (1.0 - p), 1.0 / b);
            }
            afac = -LnGamma(a) - LnGamma(b) + LnGamma(a + b);
            for (j = 0; j < 10; j++)
            {
                if (x == 0.0 || x == 1.0) return x; //a or b too small for accurate calcu-lation. 
                err = IncBeta(a, b, x) - p;
                t = Exp(a1 * Ln(x) + b1 * Ln(1.0 - x) + afac);
                u = err / t; //Halley:
                x -= (t = u / (1.0 - 0.5 * Min(1.0, u * (a1 / x - b1 / (1.0 - x)))));
                if (x <= 0.0) x = 0.5 * (x + t); //Bisect if xtries to go neg or>1.
                if (x >= 1.0) x = 0.5 * (x + t + 1.0);
                if (Abs(t) < EPS * x && j > 0) break;
            }
            return x;
        }

        static void gauleg(int n, out double[] x, out double[] w)
        {
            // Given the lower and upper limits of integrationx1andx2,andgivenn, 
            // this routine returns arraysx[1..n]andw[1..n]of lengthn, 
            // containing the abscissas and weights of the Gauss-Legendren-point quadrature formula.
            double EPS = 1e-10;
            int m, j, i;
            double x1 = -1, x2 = 1, z1, z, xm, xl, pp, p3, p2, p1; // High precision is a good idea for this rou-tine.
            x = new double[n]; w = new double[n];
            m = (n + 1) / 2; //The roots are symmetric in the interval, so
            //we only have to find half of them. 
            xm = 0.5 * (x1 + x2);
            xl = 0.5 * (x2 - x1);
            for (i = 0; i < m; i++)
            { //Loop over the desired roots.
                z = Math.Cos(Math.PI * (i + 1 - 0.25) / (n + 0.5));
                //Starting with the above approximation to theith root, we enter the main loop of
                //refinement by Newton’s method.
                do
                {
                    p1 = 1.0;
                    p2 = 0.0;
                    for (j = 1; j <= n; j++)
                    { //Loop up the recurrence relation to get the
                        //Legendre polynomial evaluated atz. 
                        p3 = p2;
                        p2 = p1;
                        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
                    }
                    //p1is now the desired Legendre polynomial. We next compute pp, its derivative,
                    //by a standard relation involving alsop2, the polynomial of one lower order.
                    pp = n * (z * p1 - p2) / (z * z - 1.0);
                    z1 = z;
                    z = z1 - p1 / pp; //Newton’s method.
                } while (Math.Abs(z - z1) > EPS);
                x[i] = xm - xl * z; // Scale the root to the desired interval,
                x[n - 1 - i] = xm + xl * z; // and put in its symmetric counterpart.
                w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp); // Compute the weight
                w[n - 1 - i] = w[i]; // and its symmetric counterpart.
            }
        }

        static double BesselJ0(double x)
        {
            double ax, xxx, xx, y, ans, ans1, ans2;
            ax = Abs(x);
            if (ax < 8.0)
            {
                y = x * x;
                ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                       + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
                ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                        + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
                ans = ans1 / ans2;
            }
            else
            {
                xxx = 8.0 / ax; y = xxx * xxx; xx = ax - 0.785398164;
                ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                     + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                       + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                        - y * 0.934945152e-7)));
                ans = Sqrt(0.636619772 / ax) * (Cos(xx) * ans1 - xxx * Sin(xx) * ans2);
            }
            return ans;
        }

        static Complex BesselJ0(Complex x)
        {
            Complex xxx, xx, y, ans, ans1, ans2;
            double ax = Abs(x);
            if (ax < 8.0)
            {
                y = x * x;
                ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                       + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
                ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                        + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
                ans = ans1 / ans2;
            }
            else
            {
                xxx = 8.0 / ax; y = xxx * xxx; xx = ax - 0.785398164;
                ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                     + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                       + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                        - y * 0.934945152e-7)));
                ans = Sqrt(0.636619772 / ax) * (Cos(xx) * ans1 - xxx * Sin(xx) * ans2);
            }
            return ans;
        }

        static double BesselJ1(double x)
        {
            double ax, xxx, xx, y, ans, ans1, ans2;
            ax = Abs(x);
            if (ax < 8.0)
            {
                y = x * x;
                ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
                       + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
                ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
                       + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
                ans = ans1 / ans2;
            }
            else
            {
                xxx = 8.0 / ax; y = xxx * xxx; xx = ax - 2.356194491;
                ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                      + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                ans2 = 0.04687499995 + y * (-0.2002690873e-3
                       + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
                ans = Sqrt(0.636619772 / ax) * (Cos(xx) * ans1 - xxx * Sin(xx) * ans2);
                if (x < 0)
                {
                    ans = -ans;
                }
            }
            return ans;
        }

        static Complex BesselJ1(Complex x)
        {
            Complex xxx, xx, y, ans, ans1, ans2;
            double ax = Abs(x);
            if (ax < 8.0)
            {
                y = x * x;
                ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
                       + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
                ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
                       + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
                ans = ans1 / ans2;
            }
            else
            {
                xxx = 8.0 / ax; y = xxx * xxx; xx = ax - 2.356194491;
                ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                      + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                ans2 = 0.04687499995 + y * (-0.2002690873e-3
                       + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
                ans = Sqrt(0.636619772 / ax) * (Cos(xx) * ans1 - xxx * Sin(xx) * ans2);
                if (x.Real < 0)
                {
                    ans = -ans;
                }
            }
            return ans;
        }

        static double BesselY0(double x)
        {
            dynamic xxx, xx, y, ans, ans1, ans2;
            if (x < 8.0)
            {
                y = x * x;
                ans1 = -2957821389.0 + y * (7062834065.0 + y * (-512359803.6
                       + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
                ans2 = 40076544269.0 + y * (745249964.8 + y * (7189466.438
                       + y * (47447.26470 + y * (226.1030244 + y * 1.0))));
                ans = (ans1 / ans2) + 0.636619772 * BesselJ0(x) * Ln(x);
            }
            else
            {
                xxx = 8.0 / x; y = xxx * xxx; xx = x - 0.785398164;
                ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                        + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                        + y * (-0.6911147651e-5 + y * (0.7621095161e-6 + y * (-0.934945152e-7))));
                ans = Sqrt(0.636619772 / x) * (Sin(xx) * ans1 + xxx * Cos(xx) * ans2);
            }
            return ans;
        }

        static Complex BesselY0(Complex x)
        {
            Complex xxx, xx, y, ans, ans1, ans2;
            if (Abs(x) < 8.0)
            {
                y = x * x;
                ans1 = -2957821389.0 + y * (7062834065.0 + y * (-512359803.6
                       + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
                ans2 = 40076544269.0 + y * (745249964.8 + y * (7189466.438
                       + y * (47447.26470 + y * (226.1030244 + y * 1.0))));
                ans = (ans1 / ans2) + 0.636619772 * BesselJ0(x) * Ln(x);
            }
            else
            {
                xxx = 8.0 / x; y = xxx * xxx; xx = x - 0.785398164;
                ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                        + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                        + y * (-0.6911147651e-5 + y * (0.7621095161e-6 + y * (-0.934945152e-7))));
                ans = Sqrt(0.636619772 / x) * (Sin(xx) * ans1 + xxx * Cos(xx) * ans2);
            }
            return ans;
        }

        static double BesselY1(double x)
        {
            double xxx, xx, y, ans, ans1, ans2;
            if (x < 8.0)
            {
                y = x * x;
                ans1 = x * (-0.4900604943e13 + y * (0.1275274390e13
                        + y * (-0.5153438139e11 + y * (0.7349264551e9
                        + y * (-0.4237922726e7 + y * 0.8511937935e4)))));
                ans2 = 0.2499580570e14 + y * (0.4244419664e12
                       + y * (0.3733650367e10 + y * (0.2245904002e8
                       + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));
                ans = (ans1 / ans2) + 0.636619772 * (BesselJ1(x) * Ln(x) - 1.0 / x);
            }
            else
            {
                xxx = 8.0 / x; y = xxx * xxx; xx = x - 2.356194491;
                ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                       + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                ans2 = 0.04687499995 + y * (-0.2002690873e-3
                       + y * (0.8449199096e-5 + y * (-0.88228987e-6
                       + y * 0.105787412e-6)));
                ans = Sqrt(0.636619772 / x) * (Sin(xx) * ans1 + xxx * Cos(xx) * ans2);
            }
            return ans;
        }

        static Complex BesselY1(Complex x)
        {
            Complex xxx, xx, y, ans, ans1, ans2;
            if (Abs(x) < 8.0)
            {
                y = x * x;
                ans1 = x * (-0.4900604943e13 + y * (0.1275274390e13
                        + y * (-0.5153438139e11 + y * (0.7349264551e9
                        + y * (-0.4237922726e7 + y * 0.8511937935e4)))));
                ans2 = 0.2499580570e14 + y * (0.4244419664e12
                       + y * (0.3733650367e10 + y * (0.2245904002e8
                       + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));
                ans = (ans1 / ans2) + 0.636619772 * (BesselJ1(x) * Ln(x) - 1.0 / x);
            }
            else
            {
                xxx = 8.0 / x; y = xxx * xxx; xx = x - 2.356194491;
                ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                       + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                ans2 = 0.04687499995 + y * (-0.2002690873e-3
                       + y * (0.8449199096e-5 + y * (-0.88228987e-6
                       + y * 0.105787412e-6)));
                ans = Sqrt(0.636619772 / x) * (Sin(xx) * ans1 + xxx * Cos(xx) * ans2);
            }
            return ans;
        }

        static double BesselI0(double x)
        {
            double ax, ans, y; // Accumulate polynomials in double precision.
            if ((ax = Abs(x)) < 3.75) //Polynomial fit.
            {
                y = x / 3.75;
                y *= y;
                ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                       + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
            }
            else
            {
                y = 3.75 / ax;
                ans = (Exp(ax) / Sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
                     + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
                     + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
                     + y * 0.392377e-2))))))));
            }
            return ans;
        }

        static Complex BesselI0(Complex x)
        {
            Complex y, ans;
            double ax; // Accumulate polynomials in double precision.
            if ((ax = Abs(x)) < 3.75) //Polynomial fit.
            {
                y = x / 3.75;
                y *= y;
                ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                       + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
            }
            else
            {
                y = 3.75 / ax;
                ans = (Exp(ax) / Sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
                     + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
                     + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
                     + y * 0.392377e-2))))))));
            }
            return ans;
        }

        static double BesselI1(double x)
        {
            double ax, ans, y;
            ax = Abs(x);
            if (ax < 3.75) // Polynomial fit.
            {
                y = x / 3.75;
                y *= y;
                ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                    + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
            }
            else
            {
                y = 3.75 / ax;
                ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
                ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
                ans *= (Exp(ax) / Sqrt(ax));
            }
            if (x < 0)
            {
                ans = -ans;
            }
            return ans;
        }

        static Complex BesselI1(Complex x)
        {
            Complex ans, y;
            double ax = Abs(x);
            if (ax < 3.75) // Polynomial fit.
            {
                y = x / 3.75;
                y *= y;
                ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                    + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
            }
            else
            {
                y = 3.75 / ax;
                ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
                ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
                ans *= (Exp(ax) / Sqrt(ax));
            }
            if (x.Real < 0)
            {
                ans = -ans;
            }
            return ans;
        }

        static double BesselK0(double x)
        {
            double y, ans; //Accumulate polynomials in double precision
            if (x <= 2.0) // Polynomial fit.
            {
                y = x * x / 4.0;
                ans = (-Ln(x / 2.0) * BesselI0(x)) + (-0.57721566 + y * (0.42278420
                    + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2
                    + y * (0.10750e-3 + y * 0.74e-5))))));
            }
            else
            {
                y = 2.0 / x;
                ans = (Exp(-x) / Sqrt(x)) * (1.25331414 + y * (-0.7832358e-1
                    + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
                    + y * (-0.251540e-2 + y * 0.53208e-3))))));
            }
            return ans;
        }

        static Complex BesselK0(Complex x)
        {
            Complex y, ans; //Accumulate polynomials in double precision
            if (Abs(x) <= 2.0) // Polynomial fit.
            {
                y = x * x / 4.0;
                ans = (-Ln(x / 2.0) * BesselI0(x)) + (-0.57721566 + y * (0.42278420
                    + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2
                    + y * (0.10750e-3 + y * 0.74e-5))))));
            }
            else
            {
                y = 2.0 / x;
                ans = (Exp(-x) / Sqrt(x)) * (1.25331414 + y * (-0.7832358e-1
                    + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
                    + y * (-0.251540e-2 + y * 0.53208e-3))))));
            }
            return ans;
        }

        static double BesselK1(double x)
        {
            double y, ans;
            if (x <= 2.0) // Polynomial fit.
            {
                y = x * x / 4.0;
                ans = (Ln(x / 2.0) * BesselI1(x)) + (1.0 / x) * (1.0 + y * (0.15443144
                    + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1
                    + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
            }
            else
            {
                y = 2.0 / x;
                ans = (Exp(-x) / Sqrt(x)) * (1.25331414 + y * (0.23498619
                    + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
                    + y * (0.325614e-2 + y * (-0.68245e-3)))))));
            }
            return ans;
        }

        static Complex BesselK1(Complex x)
        {
            Complex y, ans;
            if (Abs(x) <= 2.0) // Polynomial fit.
            {
                y = x * x / 4.0;
                ans = (Ln(x / 2.0) * BesselI1(x)) + (1.0 / x) * (1.0 + y * (0.15443144
                    + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1
                    + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
            }
            else
            {
                y = 2.0 / x;
                ans = (Exp(-x) / Sqrt(x)) * (1.25331414 + y * (0.23498619
                    + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
                    + y * (0.325614e-2 + y * (-0.68245e-3)))))));
            }
            return ans;
        }

        static double Chebev(double a, double b, double[] c, double x)
        {
            double d = 0.0, dd = 0.0, sv, y, y2;
            int j, m = c.Length;
            if ((x - a) * (x - b) > 0.0) throw new Exception("x not in range in routine chebev");
            y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a)); //Change of variable.
            for (j = m - 1; j >= 1; j--) // Clenshaw’s recurrence.
            {
                sv = d;
                d = y2 * d - dd + c[j];
                dd = sv;
            }
            return y * d - dd + 0.5 * c[0]; //Last step is different.
        }

        static void beschb(double x, out double gam1, out double gam2, out double gampl, out double gammi)
        {
            //Evaluates Γ1 and Γ2 by Chebyshev expansion for|x|≤1/2. Also returns 1/Γ(1 +x) and
            //1/Γ(1−x). If converting to double precision, set NUSE1=7, NUSE2=8.
            //double chebev(double a, double b, double[] c, float x);
            double xx;
            double[] c1 = {-1.142022680371168e0, 6.5165112670737e-3, 3.087090173086e-4,
                                  -3.4706269649e-6, 6.9437664e-9, 3.67795e-11, -1.356e-13};
            double[] c2 = {1.843740587300905e0,-7.68528408447867e-2,1.2719271366546e-3,
                              -4.9717367042e-6,-3.31261198e-8,2.423096e-10,-1.702e-13,-1.49e-15};
            xx = 8.0 * x * x - 1.0; // Multiplyxby 2tomake range be−1to 1,
            // and then apply transformation for eval-uating even Chebyshev series.
            gam1 = Chebev(-1.0, 1.0, c1, xx);
            gam2 = Chebev(-1.0, 1.0, c2, xx);
            gampl = gam2 - x * gam1;
            gammi = gam2 + x * (gam1);
        }

        static void besseljy(double x, double xnu, out double rj, out double ry, out double rjp, out double ryp)
        {
            double EPS = 1.0e-16;
            double FPMIN = 1.0e-30;
            double MAXIT = 10000;
            double XMIN = 2.0;
            double PI = Math.PI;
            /*
             Returns the Bessel functionsrj=Jν, ry=Yν and their derivativesrjp=Jν, ryp=Yν,for
             positivexand forxnu=ν≥0. The relative accuracy is within one or two significant digits
             ofEPS, except near a zero of one of the functions, where EPScontrols its absolute accuracy.
             FPMINis a number close to the machine’s smallest floating-point number. All internal arithmetic
             is in double precision.
             */

            int i, isign, l, nl;
            double a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2,
            fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl,
            rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
            temp, w, x2, xi, xi2, xmu, xmu2;

            nl = (x < XMIN ? (int)(xnu + 0.5) : Math.Max(0, (int)(xnu - x + 1.5)));
            xmu = xnu - nl;
            xmu2 = xmu * xmu;
            xi = 1.0 / x;
            xi2 = 2.0 * xi;
            w = xi2 / PI;    // The Wronskian.
            isign = 1;      // Evaluate CF1 by modified Lentz’s method (§5.2).
            // isignkeeps track of sign changes in the de-nominator.
            h = Math.Max(xnu * xi, FPMIN);
            b = xi2 * xnu;
            d = 0.0;
            c = h;
            for (i = 1; i <= MAXIT; i++)
            {
                b += xi2;
                d = b - d;
                if (Math.Abs(d) < FPMIN) d = FPMIN;
                c = b - 1.0 / c;
                if (Math.Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                del = c * d;
                h = del * h;
                if (d < 0.0) isign = -isign;
                if (Math.Abs(del - 1.0) < EPS) break;
            }

            if (i > MAXIT) throw new Exception("x too large in bessjy; try asymptotic expansion");
            rjl = isign * FPMIN; // Initialize Jν and Jν for downward recurrence.
            rjpl = h * rjl;
            rjl1 = rjl; //Store values for later rescaling.
            rjp1 = rjpl;
            fact = xnu * xi;
            for (l = nl; l >= 1; l--)
            {
                rjtemp = fact * rjl + rjpl;
                fact -= xi;
                rjpl = fact * rjtemp - rjl;
                rjl = rjtemp;
            }
            if (rjl == 0.0) rjl = EPS;
            f = rjpl / rjl; // Now have unnormalized Jµ and Jµ.
            if (x < XMIN) //Use series.
            {
                x2 = 0.5 * x;
                pimu = PI * xmu;
                fact = (Math.Abs(pimu) < EPS ? 1.0 : pimu / Math.Sin(pimu));
                d = -Math.Log(x2);
                e = xmu * d;
                fact2 = (Math.Abs(e) < EPS ? 1.0 : Math.Sinh(e) / e);
                beschb(xmu, out gam1, out gam2, out gampl, out gammi); //Chebyshev evaluation ofΓ1andΓ2.
                ff = 2.0 / PI * fact * (gam1 * Math.Cosh(e) + gam2 * fact2 * d); //f0.
                e = Math.Exp(e);
                p = e / (gampl * PI); //p0.
                q = 1.0 / (e * PI * gammi); //q0.
                pimu2 = 0.5 * pimu;
                fact3 = (Math.Abs(pimu2) < EPS ? 1.0 : Math.Sin(pimu2) / pimu2);
                r = PI * pimu2 * fact3 * fact3;
                c = 1.0;
                d = -x2 * x2;
                sum = ff + r * q;
                sum1 = p;
                for (i = 1; i <= MAXIT; i++)
                {
                    ff = (i * ff + p + q) / (i * i - xmu2);
                    c *= (d / i);
                    p /= (i - xmu);
                    q /= (i + xmu);
                    del = c * (ff + r * q);
                    sum += del;
                    del1 = c * p - i * del;
                    sum1 += del1;
                    if (Math.Abs(del) < (1.0 + Math.Abs(sum)) * EPS) break;
                }
                if (i > MAXIT) throw new Exception("bessy series failed to converge");
                rymu = -sum;
                ry1 = -sum1 * xi2;
                rymup = xmu * xi * rymu - ry1;
                rjmu = w / (rymup - f * rymu); //Equation (6.7.13).
            }
            else // Evaluate CF2 by modified Lentz’s method (§5.2).
            {
                a = 0.25 - xmu2;
                p = -0.5 * xi;
                q = 1.0;
                br = 2.0 * x;
                bi = 2.0;
                fact = a * xi / (p * p + q * q);
                cr = br + q * fact;
                ci = bi + p * fact;
                den = br * br + bi * bi;
                dr = br / den;
                di = -bi / den;
                dlr = cr * dr - ci * di;
                dli = cr * di + ci * dr;
                temp = p * dlr - q * dli;
                q = p * dli + q * dlr;
                p = temp;
                for (i = 2; i <= MAXIT; i++)
                {
                    a += 2 * (i - 1);
                    bi += 2.0;
                    dr = a * dr + br;
                    di = a * di + bi;
                    if (Math.Abs(dr) + Math.Abs(di) < FPMIN) dr = FPMIN;
                    fact = a / (cr * cr + ci * ci);
                    cr = br + cr * fact;
                    ci = bi - ci * fact;
                    if (Math.Abs(cr) + Math.Abs(ci) < FPMIN) cr = FPMIN;
                    den = dr * dr + di * di;
                    dr /= den;
                    di /= -den;
                    dlr = cr * dr - ci * di;
                    dli = cr * di + ci * dr;
                    temp = p * dlr - q * dli;
                    q = p * dli + q * dlr;
                    p = temp;
                    if (Math.Abs(dlr - 1.0) + Math.Abs(dli) < EPS) break;
                }
                if (i > MAXIT) throw new Exception("cf2 failed in bessjy");
                gam = (p - f) / q; // Equations (6.7.6) – (6.7.10).
                rjmu = Math.Sqrt(w / ((p - f) * gam + q));
                rjmu = Math.Sign(rjmu) * Math.Abs(rjl);
                rymu = rjmu * gam;
                rymup = rymu * (p + q / gam);
                ry1 = xmu * xi * rymu - rymup;
            }
            fact = rjmu / rjl;
            rj = rjl1 * fact; //Scale originalJν andJν
            rjp = rjp1 * fact;
            for (i = 1; i <= nl; i++) // Upward recurrence ofYν.
            {
                rytemp = (xmu + i) * xi2 * ry1 - rymu;
                rymu = ry1;
                ry1 = rytemp;
            }
            ry = rymu;
            ryp = xnu * xi * rymu - ry1;
        }

        static double[] gser(double a, double x)
        {
            const int ITMAX = 100;
            const double EPS = 3.0e-7;
            double sum, del, ap, gamser, gln;
            double[] ans = new double[2];
            int n;
            gln = LnGamma(a);
            ans[1] = gln;
            if (x <= 0.0)
            {
                if (x < 0.0) throw new Exception("x less than 0 in routine gser");
                gamser = 0.0;
                ans[0] = gamser;
                return ans;
            }
            else
            {
                ap = a;
                del = sum = 1.0 / a;
                for (n = 1; n <= ITMAX; n++)
                {
                    ++ap;
                    del *= x / ap;
                    sum += del;
                    if (Abs(del) < Abs(sum) * EPS)
                    {
                        gamser = sum * Math.Exp(-x + a * Math.Log(x) - (gln));
                        ans[0] = gamser;
                        return ans;
                    }
                }
                throw new Exception("a too large, ITMAX too small in routine gser");
            }
        }

        static double[] gcf(double a, double x)
        {
            const int ITMAX = 100;
            const double EPS = 3.0e-7, FPMIN = 1e-30;
            int i;
            double an, b, c, d, del, h, gln, gammcf;
            double[] ans = new double[2];
            gln = LnGamma(a);
            ans[1] = gln;
            b = x + 1.0 - a; //Set up for evaluating continued fraction by modified Lentz’s method (§5.2) with b0=0.
            c = 1.0 / FPMIN;
            d = 1.0 / b;
            h = d;
            for (i = 1; i <= ITMAX; i++)
            { //Iterate to convergence.
                an = -i * (i - a);
                b += 2.0;
                d = an * d + b;
                if (Abs(d) < FPMIN) d = FPMIN;
                c = b + an / c;
                if (Abs(c) < FPMIN) c = FPMIN;
                d = 1.0 / d;
                del = d * c;
                h *= del;
                if (Abs(del - 1.0) < EPS) break;
            }
            gammcf = Math.Exp(-x + a * Math.Log(x) - (gln)) * h;
            ans[0] = gammcf;
            return ans;
        }

        static double[] cisi(double x)
        {
            // Computes the cosine and sine integralsCi(x)andSi(x). Ci(0)is returned as a large negative
            // number and no error message is generated. Forx<0the routine returnsCi(−x)and you must
            // supply the −iπ yourself.
            double EPS = 6.0e-8; // Relative error, or absolute error near a zero ofCi(x).
            double EULER = 0.57721566; //Euler’s constant γ.
            int MAXIT = 100; // Maximum number of iterations allowed.
            double PIBY2 = 1.5707963; // π/2.
            double FPMIN = 1.0e-30; //Close to smallest representable floating-point number.
            double TMIN = 2.0; //Dividing line between using the series and continued fraction. 
            int i, k;
            bool odd;
            double a, err, fact, sign, sum, sumc = 0, sums, t, term, si, ci;
            Complex h, b, c, d, del; double[] ans = new double[2];
            t = Abs(x);
            if (t == 0.0)
            {// Special case.
                si = 0.0;
                ci = -1.0 / FPMIN;
                ans[0] = si; ans[1] = ci;
                return ans;
            }
            if (t > TMIN)
            { //Evaluate continued fraction by modified
                //Lentz’s method (§5.2). 
                b = new Complex(1.0, t);
                c = 1.0 / FPMIN;
                d = h = 1 / b;
                for (i = 2; i <= MAXIT; i++)
                {
                    a = -(i - 1) * (i - 1);
                    b += 2;
                    d = 1 / ((a * d) + b); //Denominators cannot be zero.
                    c = b + (a / c);
                    del = c * d;
                    h = h * del;
                    if (Abs(del.Real - 1.0) + Abs(del.Imaginary) < EPS) break;
                }
                if (i > MAXIT) throw new Exception("cf failed in cisi");
                h *= Complex.Cart(1, -t);
                ci = -(double) h.Real;
                si = PIBY2 + (double)h.Imaginary;
            }
            else
            { //Evaluate both series simultaneously.
                if (t < Math.Sqrt(FPMIN))
                { //Special case: avoid failure of convergence test because of underflow. sumc=0.0;
                    sums = t;
                }
                else
                {
                    sum = sums = sumc = 0.0;
                    sign = fact = 1.0;
                    odd = true;
                    for (k = 1; k <= MAXIT; k++)
                    {
                        fact *= t / k;
                        term = fact / k;
                        sum += sign * term;
                        err = term / Abs(sum);
                        if (odd)
                        {
                            sign = -sign;
                            sums = sum;
                            sum = sumc;
                        }
                        else
                        {
                            sumc = sum;
                            sum = sums;
                        }
                        if (err < EPS) break;
                        odd = !odd;
                    }
                    if (k > MAXIT) throw new Exception("maxits exceeded in cisi");
                }
                si = sums;
                ci = sumc + Math.Log(t) + EULER;
            }
            if (x < 0.0) si = -si;
            ans[0] = si; ans[1] = ci;
            return ans;
        }

        static double[] frenel(double x)
        {
            // Computes the Fresnel integralsS(x) andC(x) for all real x.
            double EPS = 6.0e-8; // Relative error, or absolute error near a zero ofCi(x).
            double EULER = 0.57721566; //Euler’s constant γ.
            int MAXIT = 100; // Maximum number of iterations allowed.
            double PIBY2 = 1.5707963; // π/2.
            double PI = 3.1415927; // π.
            double FPMIN = 1.0e-30; //Close to smallest representable floating-point number.
            double XMIN = 1.5; //Dividing line between using the series and continued fraction.
            int k, n;
            bool odd;
            double a, ax, fact, pix2, sign, sum, sumc, sums, term, test, s, c;
            Complex b, cc, d, h, del, cs; double[] ans = new double[2];
            ax = Abs(x);
            if (ax < Math.Sqrt(FPMIN))
            {// Special case.
                s = 0.0;
                c = ax;
                ans[0] = s; ans[1] = c;
                return ans;
            }
            else if (ax < XMIN)
            {
                sum = sums = 0.0;
                sumc = ax;
                sign = 1.0;
                fact = PIBY2 * ax * ax;
                odd = true;
                term = ax;
                n = 3;
                for (k = 1; k <= MAXIT; k++)
                {
                    term *= fact / k;
                    sum += sign * term / n;
                    test = Abs(sum) * EPS;
                    if (odd)
                    {
                        sign = -sign;
                        sums = sum;
                        sum = sumc;
                    }
                    else
                    {
                        sumc = sum;
                        sum = sums;
                    }
                    if (term < test) break;
                    odd = !odd;
                    n += 2;
                }
                if (k > MAXIT) throw new Exception("series failed in frenel");
                s = sums;
                c = sumc;
            }
            else
            { //Evaluate continued fraction by modified Lentz’s method (§5.2). 
                pix2 = PI * ax * ax;
                b = new Complex(1.0, -pix2);
                cc = new Complex(1.0 / FPMIN, 0.0);
                d = h = 1 / b;
                n = -1;
                for (k = 2; k <= MAXIT; k++)
                {
                    n += 2;
                    a = -n * (n + 1);
                    b = b + 4;
                    d = 1 / ((a * d) + b); //Denominators cannot be zero.
                    cc = b + a / cc;
                    del = cc * d;
                    h = h * del;
                    if (Abs(del.Real - 1.0) + Abs(del.Imaginary) < EPS) break;
                }
                if (k > MAXIT) throw new Exception("cf failed in frenel");
                h = new Complex(ax, -ax) * h;
                cs = new Complex(0.5, 0.5) * (1 - Complex.Cart(1, 0.5 * pix2) * h);
                c = cs.REAL();
                s = cs.IMAG();
            }
            if (x < 0.0)
            { //Use antisymmetry.
                c = -c;
                s = -s;
            }
            ans[0] = s; ans[1] = c;
            return ans;
        }

        static double Abs(double x)
        {
            return Math.Abs(x);
        }

        static double Abs(Complex x)
        {
            return x.Abs();
        }

        static double Exp(double x)
        {
            return Math.Exp(x);
        }

        static double Pow(double x, double n)
        {
            return Math.Pow(x, n);
        }

        static double Min(double x, double y)
        {
            return Math.Min(x, y);
        }

        static double Max(double x, double y)
        {
            return Math.Max(x, y);
        }

        static Complex Exp(Complex x)
        {
            return x.Exp();
        }

        static double Sqrt(double x)
        {
            return Math.Sqrt(x);
        }

        static Complex Sqrt(Complex x)
        {
            return x.Sqrt();
        }

        static double Sin(double x)
        {
            return Math.Sin(x);
        }

        static Complex Sin(Complex x)
        {
            return x.Sin();
        }

        static double Cos(double x)
        {
            return Math.Cos(x);
        }

        static Complex Cos(Complex x)
        {
            return x.Cos();
        }

        static double Sinh(double x)
        {
            return Math.Sinh(x);
        }

        static Complex Sinh(Complex x)
        {
            return x.Sinh();
        }

        static double Cosh(double x)
        {
            return Math.Cosh(x);
        }

        static Complex Cosh(Complex x)
        {
            return x.Cosh();
        }

        static double Ln(double x)
        {
            return Math.Log(x);
        }

        static Complex Ln(Complex x)
        {
            return x.Log();
        }

        #endregion

    }


    class SpecialFunctionsException : Exception
    {
        /// <summary>
        /// Exception from the Matrix Class
        /// </summary>
        /// <param name="Message">Message to be displayed</param>
        public SpecialFunctionsException(string Message)
            : base(Message)
        { }
    }
}
