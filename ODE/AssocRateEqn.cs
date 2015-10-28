using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AlgLib.LsFit.ODE
{
    /// <summary>
    /// An attempt to model a 1:1 Association via the Rate equation
    /// 
    /// based on example code from http://www.alglib.net/translator/man/manual.csharp.html#unit_odesolver
    /// </summary>
    class AssocRateEqn
    {

        /// <summary>
        /// returns a rate constant (RU/sec) 
        /// </summary>
        /// <param name="conc_M"></param>
        /// <param name="ka_perMSec"></param>
        /// <param name="ru"></param>
        /// <param name="rMax"></param>
        /// <returns></returns>
        static double calcAssocRate(double conc_M, double ka_perMSec, double kd_perSec, double ru, double rMax)
        {
            return ka_perMSec * conc_M * (rMax - ru) - kd_perSec * ru;
        }

        private static double _ka_perMSec = 4.19E+03; //From CAII
        private static double _kd_perSec = 8.51E-02;
        private static double _conc_M = 11.0E-06; // 11 uM
        private static double _rmax = 35.86;

        public static void ode_simpleAssoc(double[] y, double x, double[] dy, object obj)
        {
            // this callback calculates f(y[],x)=-y[0]

            System.Diagnostics.Debug.Assert(y.Length == 1 && dy.Length == 1, "y should be length 1 until we add multiple equations");

            dy[0] = calcAssocRate(_conc_M, _ka_perMSec, _kd_perSec, y[0], _rmax);
        }



        /// <summary>
        /// 
        /// A = Conc
        /// 
        /// B[0] = Rmax
        /// dB/dt = - (ka*A*B - kd*AB)
        /// 
        /// AB[0] = 0
        /// dAB/dt = (ka* A*B - kd* AB)
        /// </summary>
        /// <param name="y"></param>
        /// <param name="x"></param>
        /// <param name="dy"></param>
        /// <param name="obj"></param>
        public static void ode_LangmuirAssoc(double[] y, double x, double[] dy, object obj)
        {
            // this callback calculates f(y[],x)=-y[0]

            System.Diagnostics.Debug.Assert(y.Length == 2 && dy.Length == 2, "y should be length 2 with Langmuir eqn");

            //y[0] = B (a.k.a. rmax - R)
            //y[1] = RU

            //dy[0] = dB/dt
            //dy[1] = dRU/dt

            dy[0] = -1 * (_ka_perMSec * _conc_M * y[0] - _kd_perSec * y[1]);
            dy[1] = _ka_perMSec * _conc_M * y[0] - _kd_perSec * y[1];
        }

        public static void DoSomething()
        {
            double[] y = new double[] { _rmax, 0 }; //start with _rMax sites remaining, and 0 RU already bound
            double[] x = Enumerable.Range(1, 120).Select(i => (double)i / 2).ToArray();
            double eps = 0.00001;
            double h = 0;
            alglib.odesolverstate s;
            int m;
            double[] xtbl;
            double[,] ytbl;
            alglib.odesolverreport rep;
            alglib.odesolverrkck(y, x, eps, h, out s);
            alglib.odesolversolve(s, ode_LangmuirAssoc, null);
            alglib.odesolverresults(s, out m, out xtbl, out ytbl, out rep);


            System.Diagnostics.Debug.Assert(Math.Abs(10.8377 - ytbl[30, 1]) <= 1e-5, "At 15.0 sec, we should have 10.83 RU of bound CAII");

            System.Console.WriteLine("{0}", m); // EXPECTED: 4
            System.Console.WriteLine("{0}", alglib.ap.format(xtbl, 2)); // EXPECTED: [0, 1, 2, 3]
            System.Console.WriteLine("{0}", alglib.ap.format(ytbl, 2)); // EXPECTED: [[1], [0.367], [0.135], [0.050]]
            //System.Console.ReadLine();
        }
    }
}
