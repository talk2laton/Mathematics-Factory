using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    public class Statistics
    {
        public class Normal
        {
            double mu, sig;
            public Normal(double mmu = 0, double ssig = 0)
            {
                //Constructor. Initialize withand. The default with no arguments is N.0;1/
                mu = mmu;
                sig = ssig;
                if (sig <= 0.0) throw new StatisticException("bad sig in Normaldist");
            }

            public double P(double x)
            {
                // Return probability density function.
                return (0.398942280401432678 / sig) * Exp(-0.5 * SQR((x - mu) / sig));
            }

            public double Cdf(double x)
            {
                // Return cumulative distribution function.
                return 0.5 * SpecialFunctions.Erfc(-0.707106781186547524 * (x - mu) / sig);

            }

            public double InvCdf(double P)
            {
                //Return inverse cumulative distribution function.
                if (P <= 0.0 || P >= 1.0)
                    throw new StatisticException("bad p in Normaldist");
                return -1.41421356237309505 * sig * SpecialFunctions.InvErf(1 - 2.0 * P) + mu;
            }
        }

        public class Cauchy
        {
            double mu, sig, OneOnPI = 0.318309886183791333066907387547;
            public Cauchy(double mmu = 0, double ssig = 0)
            {
                //Constructor. Initialize withand. The default with no arguments is N.0;1/
                mu = mmu;
                sig = ssig;
                if (sig < 0.0) throw new StatisticException("bad sig in Cauchydist");
            }

            public double P(double x)
            {
                // Return probability density function.
                return OneOnPI / (sig * (1.0 + SQR((x - mu) / sig)));
            }

            public double Cdf(double x)
            {
                // Return cumulative distribution function.
                return 0.5 + OneOnPI * Atan2(x - mu, sig);
            }

            public double InvCdf(double P)
            {
                //Return inverse cumulative distribution function.
                if (P <= 0.0 || P >= 1.0)
                    throw new StatisticException("bad P in Cauchydist");
                return mu + sig * Tan(Math.PI * (P - 0.5));
            }
        }

        public class StudentT
        {
            double nu, mu, sig, np, fac;
            public StudentT(double nnu, double mmu = 0, double ssig = 0)
            {
                //Constructor. Initialize withand. The default with no arguments is N.0;1/
                nu = nnu;
                mu = mmu;
                sig = ssig;

                if (sig <= 0.0 || nu <= 0.0) throw new StatisticException("bad sig, nu in Studentdist");
                np = 0.5 * (nu + 1.0);
                fac = SpecialFunctions.LnGamma(np) - SpecialFunctions.LnGamma(0.5 * nu);
            }

            public double P(double t)
            {
                // Return probability density function.
                return Exp(-np * Ln(1.0 + SQR((t - mu) / sig) / nu) + fac) / (Sqrt(Math.PI * nu) * sig);
            }

            public double Cdf(double t)
            {
                // Return cumulative distribution function.
                double p = 0.5 * SpecialFunctions.IncBeta(0.5 * nu, 0.5, nu / (nu + SQR((t - mu) / sig)));
                if (t >= mu) return 1.0 - p;
                else return p;
            }

            public double InvCdf(double P)
            {
                //Return inverse cumulative distribution function.
                if (P <= 0.0 || P >= 1.0)
                    throw new StatisticException("bad P in Studentdist");
                double x = SpecialFunctions.InvIncBeta(2.0 * Min(P, 1.0 - P), 0.5 * nu, 0.5);
                x = sig * Sqrt(nu * (1.0 - x) / x);
                return (P >= 0.5 ? mu + x : mu - x);
            }

            public double Aa(double t)
            {
                // Return the two-tailed cdf A(t/v).
                if (t < 0.0) throw new StatisticException("bad t in Studentdist");
                return 1.0 - SpecialFunctions.IncBeta(0.5 * nu, 0.5, nu / (nu + SQR(t)));
            }

            public double InvAa(double P)
            {
                //Return the inverse, namelyt such that p= A(t/v)
                if (P <= 0.0 || P >= 1.0)
                    throw new StatisticException("bad P in Studentdist");
                double x = SpecialFunctions.InvIncBeta(1.0 - P, 0.5 * nu, 0.5);
                return Sqrt(nu * (1.0 - x) / x);
            }
        }

        public class Logistic
        {
            //Logistic distribution.
            double mu, sig;
            public Logistic(double mmu = 0.0, double ssig = 1.0)
            {
                //Constructor. Initialize withand. The default with no arguments is Logistic(0,1)/
                mu = mmu;
                sig = ssig;
                if (sig <= 0.0) throw new StatisticException("bad sig in Logisticdist");
            }

            public double P(double x)
            {
                //Return probability density function.
                double e = Exp(-Abs(1.81379936423421785 * (x - mu) / sig));
                return 1.81379936423421785 * e / (sig * SQR(1.0 + e));
            }

            public double Cdf(double x)
            {
                //Return cumulative distribution function.
                double e = Exp(-Abs(1.81379936423421785 * (x - mu) / sig));
                if (x >= mu) return 1.0 / (1.0 + e); //Because we usedabsto control over-flow, we now have two cases. 
                else return e / (1.0 + e);
            }

            public double InvCdf(double p)
            {
                //Return inverse cumulative distribution function.
                if (p <= 0.0 || p >= 1.0)
                    throw new StatisticException("bad p in Logisticdist");
                return mu + 0.551328895421792049 * sig * Ln(p / (1.0 - p));
            }
        }

        public class Exponential
        {
            //Exponential distribution.
            double bet;
            public Exponential(double bbet)
            {
                //Constructor Initialize with beta
                bet = bbet;
                if (bet <= 0.0)
                    throw new StatisticException("bad bet in Expondist");
            }

            public double P(double x)
            {
                //Return probability density function.
                if (x < 0.0)
                    throw new StatisticException("bad x in Expondist");
                return bet * Exp(-bet * x);
            }

            public double Cdf(double x)
            {
                //Return cumulative distribution function.
                if (x < 0.0)
                    throw new StatisticException("bad x in Expondist");
                return 1.0 - Exp(-bet * x);
            }

            public double InvCdf(double p)
            {
                //Return inverse cumulative distribution function.
                if (p < 0.0 || p >= 1.0)
                    throw new StatisticException("bad p in Expondist");
                return -Ln(1.0 - p) / bet;
            }
        }

        public class LogNormal
        {
            //Lognormal distribution, derived from the error function Erf.
            double mu, sig;
            public LogNormal(double mmu = 0.0, double ssig = 1.0)
            {
                mu = mmu;
                sig = ssig;
                if (sig <= 0.0) throw new StatisticException("bad sig in Lognormaldist");
            }

            public double P(double x)
            {
                //Return probability density function.
                if (x < 0.0) throw new StatisticException("bad x in Lognormaldist");
                if (x == 0.0) return 0.0;
                return (0.398942280401432678 / (sig * x)) * Exp(-0.5 * SQR((Ln(x) - mu) / sig));
            }

            public double Cdf(double x)
            {
                //Return cumulative distribution function.
                if (x < 0.0) throw new StatisticException("bad x in Lognormaldist");
                if (x == 0.0) return 0.0;
                return 0.5 * SpecialFunctions.Erfc(-0.707106781186547524 * (Ln(x) - mu) / sig);
            }

            public double InvCdf(double p)
            {
                //Return inverse cumulative distribution function.
                if (p <= 0.0 || p >= 1.0) throw new StatisticException("bad p in Lognormaldist");
                return Exp(-1.41421356237309505 * sig * SpecialFunctions.InvErf(1 - 2.0 * p) + mu);
            }
        }

        public class ChiSquare
        {
            //distribution, derived from the gamma function Gamma.
            double nu, fac;
            public ChiSquare(double nnu)
            {
                //Constructor. Initialize with.
                nu = nnu;
                if (nu <= 0.0) throw new StatisticException("bad nu in Chisqdist");
                fac = 0.693147180559945309 * (0.5 * nu) + SpecialFunctions.LnGamma(0.5 * nu);
            }

            public double P(double x2)
            {
                //Return probability density function.
                if (x2 <= 0.0)
                    throw new StatisticException("bad x2 in Chisqdist");
                return Exp(-0.5 * (x2 - (nu - 2.0) * Ln(x2)) - fac);
            }

            public double Cdf(double x2)
            {
                //Return cumulative distribution function.
                if (x2 < 0.0)
                    throw new StatisticException("bad x2 in Chisqdist");
                return SpecialFunctions.GammaP(0.5 * nu, 0.5 * x2);
            }

            public double InvCdf(double P)
            {
                //Return inverse cumulative distribution function.
                if (P < 0.0 || P >= 1.0)
                    throw new StatisticException("bad p in Chisqdist");
                return 2.0 * SpecialFunctions.InvGammaP(P, 0.5 * nu);
            }
        }

        public class Gamma
        {
            //Gamma distribution, derived from the gamma function Gamma.
            double alph, bet, fac;
            public Gamma(double aalph, double bbet = 1.0)
            {
                //Constructor. Initialize with˛ and .
                if (alph <= 0.0 || bet <= 0.0) throw new StatisticException("bad alph,bet in Gammadist");
                alph = aalph;
                bet = bbet;
                fac = alph * Ln(bet) - SpecialFunctions.LnGamma(alph);
            }
            public double P(double x)
            {
                //Return probability density function.
                if (x <= 0.0) throw new StatisticException("bad x in Gammadist");
                return Exp(-bet * x + (alph - 1.0) * Ln(x) + fac);
            }

            public double Cdf(double x)
            {
                //Return cumulative distribution function.
                if (x < 0.0) throw new StatisticException("bad x in Gammadist");
                return SpecialFunctions.GammaP(alph, bet * x);
            }

            public double InvCdf(double p)
            {
                //Return inverse cumulative distribution function.
                if (p < 0.0 || p >= 1.0) throw new StatisticException("bad p in Gammadist");
                return SpecialFunctions.InvGammaP(p, alph) / bet;
            }
        }

        public class F
        {
            //Fdistribution, derived from the beta functionBeta.
            double nu1, nu2, fac;
            public F(double nnu1, double nnu2)
            {
                //Constructor. Initialize with nu1 and nu2.
                if (nu1 <= 0.0 || nu2 <= 0.0) throw new StatisticException("bad nu1,nu2 in Fdist");
                fac = 0.5 * (nu1 * Ln(nu1) + nu2 * Ln(nu2)) + SpecialFunctions.LnGamma(0.5 * (nu1 + nu2))
                    - SpecialFunctions.LnGamma(0.5 * nu1) - SpecialFunctions.LnGamma(0.5 * nu2);
            }

            public double P(double f)
            {
                //Return probability density function.
                if (f <= 0.0) throw new StatisticException("bad f in Fdist");
                return Exp((0.5 * nu1 - 1.0) * Ln(f) - 0.5 * (nu1 + nu2) * Ln(nu2 + nu1 * f) + fac);
            }

            public double Cdf(double f)
            {
                //Return cumulative distribution function.
                if (f < 0.0) throw new StatisticException("bad f in Fdist");
                return SpecialFunctions.IncBeta(0.5 * nu1, 0.5 * nu2, nu1 * f / (nu2 + nu1 * f));
            }

            public double InvCdf(double p)
            {
                //Return inverse cumulative distribution function.
                if (p <= 0.0 || p >= 1.0) throw new StatisticException("bad p in Fdist");
                double x = SpecialFunctions.InvIncBeta(p, 0.5 * nu1, 0.5 * nu2);
                return nu2 * x / (nu1 * (1.0 - x));
            }
        }

        public class Beta
        {
            //Beta distribution, derived from the beta function Beta.
            double alph, bet, fac;
            public Beta(double aalph, double bbet)
            {
                //Constructor. Initialize with alpha and beta
                alph = aalph;
                bet = bbet;
                if (alph <= 0.0 || bet <= 0.0)
                    throw new StatisticException("bad alph,bet in Betadist");
                fac = SpecialFunctions.LnGamma(alph + bet) -
                    SpecialFunctions.LnGamma(alph) - SpecialFunctions.LnGamma(bet);
            }

            public double P(double x)
            {
                //Return probability density function.
                if (x <= 0.0 || x >= 1.0)
                    throw new StatisticException("bad x in Betadist");
                return Exp((alph - 1.0) * Ln(x) + (bet - 1.0) * Ln(1.0 - x) + fac);
            }

            public double Cdf(double x)
            {
                //Return cumulative distribution function.
                if (x < 0.0 || x > 1.0) throw new StatisticException("bad x in Betadist");
                return SpecialFunctions.IncBeta(alph, bet, x);
            }

            public double InvCdf(double p)
            {
                //Return inverse cumulative distribution function.
                if (p < 0.0 || p > 1.0) throw new StatisticException("bad p in Betadist");
                return SpecialFunctions.InvIncBeta(p, alph, bet);
            }
        }

        public class Binomial
        {
            //Binomial distribution, derived from the beta functionBeta.
            int n;
            double pe, fac;
            public Binomial(int nn, double ppe)
            {
                //Constructor. Initialize withn(sample size) and p(event probability).
                n = nn;
                pe = ppe;
                if (n <= 0 || pe <= 0.0 || pe >= 1.0) throw new StatisticException("bad args in Binomialdist");
                fac = SpecialFunctions.LnGamma(n + 1.0);
            }
            public double P(int k)
            {
                //Return probability density function.
                if (k < 0) throw new StatisticException("bad k in Binomialdist");
                if (k > n) return 0.0;
                return Exp(k * Ln(pe) + (n - k) * Ln(1.0 - pe) + fac -
                    SpecialFunctions.LnGamma(k + 1.0) - SpecialFunctions.LnGamma(n - k + 1.0));
            }

            public double Cdf(int k)
            {
                //Return cumulative distribution function.
                if (k < 0) throw new StatisticException("bad k in Binomialdist");
                if (k == 0) return 0.0;
                if (k > n) return 1.0;
                return 1.0 - SpecialFunctions.IncBeta(k, n - k + 1.0, pe);
            }

            public int InvCdf(double p)
            {
                //Given argumentP, return integer nsuch thatP.< n/PP.< nC1/.
                int k, kl, ku, inc = 1;
                if (p <= 0.0 || p >= 1.0) throw new StatisticException("bad p in Binomialdist");
                k = (int)Max(0, Min(n, (n * pe))); //Starting guess near peak of density.
                if (p < Cdf(k))
                {
                    //Expand interval until we bracket.
                    do
                    {
                        k = (int)Max(k - inc, 0);
                        inc *= 2;
                    } while (p < Cdf(k));
                    kl = k; ku = k + inc / 2;
                }
                else
                {
                    do
                    {
                        k = (int)Min(k + inc, n + 1);
                    } while (p > Cdf(k));
                    ku = k; kl = k - inc / 2;
                }

                while (ku - kl > 1)
                {
                    // Now contract the interval by bisection.
                    k = (kl + ku) / 2;
                    if (p < Cdf(k)) ku = k;
                    else kl = k;
                }
                return kl;
            }
        }

        public class Poisson
        {
            //Poisson distribution, derived from the gamma function Gamma.
            double lam;
            public Poisson(double llam)
            {
                //Constructor. Initialize with.
                lam = llam;
                if (lam <= 0.0) throw new StatisticException("bad lam in Poissondist");
            }

            public double P(int n)
            {
                //Return probability density function.
                if (n < 0) throw new StatisticException("bad n in Poissondist");
                return Exp(-lam + n * Ln(lam) - SpecialFunctions.LnGamma(n + 1.0));
            }

            public double Cdf(int n)
            {
                //Return cumulative distribution function.
                if (n < 0) throw new StatisticException("bad n in Poissondist");
                if (n == 0) return 0.0;
                return SpecialFunctions.GammaQ(n, lam);
            }

            public int InvCdf(double p)
            {
                //Given argument P, return integer n such that P(< n)<=P <= P(< n + 1).
                int n, nl, nu, inc = 1;
                if (p <= 0.0 || p >= 1.0) throw new StatisticException("bad p in Poissondist");
                if (p < Exp(-lam)) return 0;
                n = (int)Max(Sqrt(lam), 5.0); //Starting guess near peak of density.
                if (p < Cdf(n))
                { //Expand interval until we bracket.
                    do
                    {
                        n = (int)Max(n - inc, 0);
                        inc *= 2;
                    } while (p < Cdf(n));
                    nl = n; nu = n + inc / 2;
                }
                else
                {
                    do
                    {
                        n += inc;
                        inc *= 2;
                    } while (p > Cdf(n));
                    nu = n; nl = n - inc / 2;
                }

                while (nu - nl > 1)
                {
                    // Now contract the interval by bisection.
                    n = (nl + nu) / 2;
                    if (p < Cdf(n)) nu = n;
                    else nl = n;
                }
                return nl;
            }
        }

        public class Triangular
        {
            double min, likely, max, h, threshold;
            public Triangular(double Min = 0, double Likely = 0.5, double Max = 1)
            {
                //Constructor. Initialize withand. The default with no arguments is N.0;1/
                min = Min;
                likely = Likely;
                max = Max;
                if (likely < min || likely > max) throw new StatisticException("bad sig in Trianglardist");
                h = 2 / (Max - Min); threshold = (Likely - Min) / (Max - Min);
            }

            public double P(double x)
            {
                // Return probability density function.
                if (x <= likely)
                {
                    return (x - min) * h / (likely - min);
                }
                else
                {
                    return (max - x) * h / (max - likely);
                }
            }

            public double Cdf(double x)
            {
                // Return cumulative distribution function.
                if (x <= likely)
                {
                    return 0.5 * Math.Pow(x - min, 2) * h / (likely - min);
                }
                else
                {
                    return 1 - 0.5 * Math.Pow(max - x, 2) * h / (max - likely);
                }

            }

            public double InvCdf(double P)
            {
                //Return inverse cumulative distribution function.
                if (P <= threshold)
                {
                    return min + Math.Sqrt(2 * P * (likely - min) / h);
                }
                else
                {
                    return max - Math.Sqrt(2 * (1 - P) * (max - likely) / h);
                }
            }
        }

        public class Uniform
        {
            double min, max, h;
            public Uniform(double Min = 0, double Max = 1)
            {
                //Constructor. Initialize withand. The default with no arguments is N.0;1/
                min = Min;
                max = Max;
                if (min > max) throw new StatisticException("bad sig in Uniformdist");
                h = 1 / (Max - Min);
            }

            public double P(double x)
            {
                // Return probability density function.
                return h;
            }

            public double Cdf(double x)
            {
                // Return cumulative distribution function.
                return (x - min) / (max - min);
            }

            public double InvCdf(double P)
            {
                //Return inverse cumulative distribution function.
                return min + P * (max - min);
            }
        }

        #region private
        static double SQR(double x)
        {
            return x * x;
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

        static double Tan(double x)
        {
            return Math.Tan(x);
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

        static double Atan2(double y, double x)
        {
            return Math.Atan2(y, x);
        }

        static Complex Ln(Complex x)
        {
            return x.Log();
        }

        static double Min(double a1, double a2)
        {
            return Math.Min(a1, a2);
        }

        static double Max(double a1, double a2)
        {
            return Math.Max(a1, a2);
        }

        #endregion

    }

    class StatisticException : Exception
    {
        /// <summary>
        /// Exception from the Matrix Class
        /// </summary>
        /// <param name="Message">Message to be displayed</param>
        public StatisticException(string Message)
            : base(Message)
        { }
    }
}
