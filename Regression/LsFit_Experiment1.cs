using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BiOptix.AlgLib
{
    /// <summary>
    /// Exapmle code to excercise the Nonliear Least Squares solver in AlgLib's LsFt subpackage, using the "F Only" method (no analytic gradient or hessian)
    /// 
    /// ref: http://www.alglib.net/translator/man/manual.csharp.html#example_lsfit_d_nlf
    /// </summary>
    class LsFit_Experiment1
    {
        /// <summary>
        /// called for each point in each iteration of the levenberg, calculates model points based on each input point x and the parameters c
        /// </summary>
        /// <param name="c">parameter vector</param>
        /// <param name="x">a single data point from the x data matrix, where n points exist in m dimensions, here x is a m-length row vector (note: when dealing with time series data, m is frequently 1)</param>
        /// <param name="func"></param>
        /// <param name="obj"></param>
        public static void function_cx_1_func(double[] c, double[] x, ref double func, object obj)
        {
            // this callback calculates f(c,x)=exp(-c0*sqr(x0))
            // where x is a position on X-axis and c is adjustable parameter
            func = System.Math.Exp(-c[0] * x[0] * x[0]);
        }

        /// <summary>
        /// Run the example
        /// </summary>
        public static void Run_lsfit_d_nlf()
        {
            //
            // In this example we demonstrate exponential fitting
            // by f(x) = exp(-c*x^2)
            // using function value only.
            //
            // Gradient is estimated using combination of numerical differences
            // and secant updates. diffstep variable stores differentiation step 
            // (we have to tell algorithm what step to use).
            //
            double[,] x = new double[,] { { -1 }, { -0.8 }, { -0.6 }, { -0.4 }, { -0.2 }, { 0 }, { 0.2 }, { 0.4 }, { 0.6 }, { 0.8 }, { 1.0 } };
            double[] y = new double[] { 0.223130, 0.382893, 0.582748, 0.786628, 0.941765, 1.000000, 0.941765, 0.786628, 0.582748, 0.382893, 0.223130 };
            double[] c = new double[] { 0.3 };
            double epsf = 0;
            double epsx = 0.000001;
            int maxits = 0;
            int info;
            alglib.lsfitstate state;
            alglib.lsfitreport rep;
            double diffstep = 0.0001;

            //
            // Fitting without weights
            //
            alglib.lsfitcreatef(x, y, c, diffstep, out state);
            alglib.lsfitsetcond(state, epsf, epsx, maxits);
            alglib.lsfitfit(state, function_cx_1_func, null, null);
            alglib.lsfitresults(state, out info, out c, out rep);
            System.Console.WriteLine("{0}", info); // EXPECTED: 2
            System.Console.WriteLine("{0}", alglib.ap.format(c, 1)); // EXPECTED: [1.5]

            //
            // Fitting with weights
            // (you can change weights and see how it changes result)
            //
            double[] w = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            alglib.lsfitcreatewf(x, y, w, c, diffstep, out state);
            alglib.lsfitsetcond(state, epsf, epsx, maxits);
            alglib.lsfitfit(state, function_cx_1_func, null, null);
            alglib.lsfitresults(state, out info, out c, out rep);
            System.Console.WriteLine("{0}", info); // EXPECTED: 2
            System.Console.WriteLine("{0}", alglib.ap.format(c, 1)); // EXPECTED: [1.5]
            System.Console.ReadLine();
            //return 0;
        }
    }
    

}
