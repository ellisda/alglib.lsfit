using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AlgLib.LsFit.ODE
{
    /// <summary>
    /// This is example ODE Solver code from http://www.alglib.net/translator/man/manual.csharp.html#unit_odesolver
    /// </summary>
    class Example_ODE
    {
        public static void ode_function_1_diff(double[] y, double x, double[] dy, object obj)
        {
            // this callback calculates f(y[],x)=-y[0]
            //dy[0] = -y[0];

            // returning -y[0] is same as computing f`(x) = -1 * f(x)  -- solution to 1st deriv of f(x) = exp(-x) -> f`(x) = -1 * exp(-x)

            //other: if f(x) = 2x + 12, then f'(x) = 2
            //dy[0] = 2;


            //other2: if f(x) = 5x^2, f'(x) = 10x
            dy[0] = 10 * x;
        }

        public static void DoSomething()
        {
            double[] y = new double[] { 1 };
            double[] x = new double[] { 0, 1, 2, 3 };
            double eps = 0.00001;
            double h = 0;
            alglib.odesolverstate s;
            int m;
            double[] xtbl;
            double[,] ytbl;
            alglib.odesolverreport rep;
            alglib.odesolverrkck(y, x, eps, h, out s);
            alglib.odesolversolve(s, ode_function_1_diff, null);
            alglib.odesolverresults(s, out m, out xtbl, out ytbl, out rep);
            System.Console.WriteLine("{0}", m); // EXPECTED: 4
            System.Console.WriteLine("{0}", alglib.ap.format(xtbl, 2)); // EXPECTED: [0, 1, 2, 3]
            System.Console.WriteLine("{0}", alglib.ap.format(ytbl, 2)); // EXPECTED: [[1], [0.367], [0.135], [0.050]]
            //System.Console.ReadLine();
        }
    }
}
