using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AlgLib.LsFit.ODE
{
    /// <summary>
    /// An attempt to model a Mass Transport Model with the 3 part rate equations
    ///     NOTE: The modeling of the [A] (analyte concentration in solution, at slide surface) is very fast and seems to require a stiff ode solver
    ///           while the other 2 parts of non-stiff. Kevin tried using ode45 and ode15 and others.
    ///           We're interested in such a 3 part equation struggles to predict mass transport data when coded as 3-part ODE
    /// 
    /// 
    /// REF: https://bioptix.sharepoint.com/A4softwareupdates/A4SoftwareDevelopment/_layouts/OneNote.aspx?id=%2FA4softwareupdates%2FA4SoftwareDevelopment%2FDocuments%2FA4%20Data%20Artifacts%2C%20Test%20Plans%2C%20and%20Resulting%20Software%20Changes&wd=target%28BiOptix%20Analysis%20Tool.one%7C1339185A-20B7-4CD1-BB19-F1EAF84D1C1F%2F1%3A1%20binding%20with%20mass%20transfer%7CF7E6CF73-C7F1-4123-B401-FD57BF2EBD08%2F%29
    /// A(solution) = Conc
    /// 
    /// A[0] = 0
    /// dA/dt = kt*(Conc-A) - (ka* A* B - kd* AB)
    /// 
    /// B[0] = Rmax
    /// dB/dt = - (ka* A* B - kd* AB)
    /// 
    /// AB[0] = 0
    /// dAB/dt = (ka* A* B - kd* AB)
    /// 
    /// Total response:
    /// AB + RI
    /// 
    /// 
    /// based on example code from http://www.alglib.net/translator/man/manual.csharp.html#unit_odesolver
    /// </summary>
    class MassTransportThreeRateEqn
    {
        /// <summary>
        /// enum to help book keeping in y[] and dy[] indices
        /// </summary>
        internal enum Eqn : int { A=0, B=1, AB=2 }

        /// <summary>
        /// 
        /// A[0] = 0
        /// dA/dt = kt*(Conc-A) - (ka* A* B - kd* AB)
        /// 
        /// B[0] = Rmax
        /// dB/dt = - (ka* A* B - kd* AB)
        /// 
        /// AB[0] = 0
        /// dAB/dt = (ka* A* B - kd* AB)
        /// </summary>
        public static void ode_massTransportAssoc(double[] y, double x, double[] dy, object obj)
        {
            System.Diagnostics.Debug.Assert(y.Length == 3 && dy.Length == 3, "y should be length 3 with mass transport 3-part eqn");
            int A = (int)Eqn.A;
            int B = (int)Eqn.B;
            int AB = (int)Eqn.AB;

            dy[A] = _kt * (_conc_M - y[A]) - (_ka_perMSec * y[A] * y[B] - _kd_perSec * y[AB]);

            dy[B] = -1 * (_ka_perMSec * y[A] * y[B] - _kd_perSec * y[AB]);

            dy[AB] = -1 * dy[B];
        }

        private static double _ka_perMSec = 1e3;
        private static double _kd_perSec = 1e-3;
        private static double _rmax = 100;
        private static double _conc_M = 1e-4;

        private static double _kt = 1e8;

        public static void DoSomething_TakeForever()
        {
            double[] y = new double[3];
            y[(int)Eqn.A] = 0;
            y[(int)Eqn.B] = _rmax;
            y[(int)Eqn.AB] = 0;
            double[] x = Enumerable.Range(1, 120).Select(i => (double)i / 2).ToArray();
            double eps = 0.00001;
            double h = 0;
            alglib.odesolverstate s;
            int m;
            double[] xtbl;
            double[,] ytbl;
            alglib.odesolverreport rep;
            alglib.odesolverrkck(y, x, eps, h, out s);
            alglib.odesolversolve(s, ode_massTransportAssoc, null);
            alglib.odesolverresults(s, out m, out xtbl, out ytbl, out rep);


            //System.Diagnostics.Debug.Assert(Math.Abs(10.8377 - ytbl[30, 1]) <= 1e-5, "At 15.0 sec, we should have 10.83 RU of bound CAII");

            System.Console.WriteLine("{0}", m); // EXPECTED: 4
            System.Console.WriteLine("{0}", alglib.ap.format(xtbl, 2)); // EXPECTED: [0, 1, 2, 3]
            System.Console.WriteLine("{0}", alglib.ap.format(ytbl, 2)); // EXPECTED: [[1], [0.367], [0.135], [0.050]]
            //System.Console.ReadLine();
        }
    }
}
