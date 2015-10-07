using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BiOptix.AlgLib
{
    /// <summary>
    /// An attempt to use LSFitFit to solve for 5PL parameters with Concentration Analysis data and return paremeter error for purpose of confidence intervals
    ///
    /// base on LsFit_Experiment1 from ref: http://www.alglib.net/translator/man/manual.csharp.html#example_lsfit_d_nlf
    /// </summary>
    class LsFit_Experiment2
    {

        /// <summary>
        /// calculate 5PL at value x with parameters c
        /// </summary>
        /// <param name="c"></param>
        /// <param name="x">vector whose first and only value is the x value that we use to calculate y, f(x, a,b,c,d,g) = y </param>
        /// <param name="func"></param>
        /// <param name="obj"></param>
        public static void f_5PL(double[] c, double[] x, ref double func, object obj)
        {
            System.Diagnostics.Debug.Assert(c.Length == 5, "5PL needs 5 parameters");
            func = alglib.logisticcalc5(x[0], c[0], c[1], c[2], c[3], c[4]);
        }

        /// <summary>
        /// Customimzed example of lsfit_d_nlf, with some concentration analysis data and the 5PL function, boundary constrained to ensure 5PL is defined
        /// </summary>
        public static void Run_lsfit_5PL()
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
            //double[,] x = new double[,] { { -1 }, { -0.8 }, { -0.6 }, { -0.4 }, { -0.2 }, { 0 }, { 0.2 }, { 0.4 }, { 0.6 }, { 0.8 }, { 1.0 } };
            //double[] y = new double[] { 0.223130, 0.382893, 0.582748, 0.786628, 0.941765, 1.000000, 0.941765, 0.786628, 0.582748, 0.382893, 0.223130 };

            //Data from ABS (data that DoubleReferenceGUI starts with, ch3, DMSO-disabled, standard blanking, Bound Avg [10,30]
            double[,] x_conc_M = new double[,]{
            {1.23E-06},
            {3.70E-06},
            {1.10E-05},
            {1.10E-05},
            {3.30E-05},
            {3.30E-05},
            {1.00E-04}};

            double[] y_RU = new double[]{
            0.32110643,
            2.335124169,
            5.0964102,
            5.361168514,
            8.562463415,
            8.5834102,
            11.32491796
            };

            double[] c_5PL = new double[] { -10, 0.5, 0.003, 20, 1.0};

            //double[] c_5PL = new double[] /*VALUES FROM Existin Algorithm*/ { -2.22332885169764, 0.629493994723955, 9.07171520030884E-05, 13.4940430048216, 2.73435657052643 };

            //Apparrently the 5PL calc function requires c to be positive
            double[] c_5PL_low = new double[]{double.NegativeInfinity, double.NegativeInfinity, 1e-14, double.NegativeInfinity, double.NegativeInfinity};
            double[] c_5PL_High = new double[]{double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity};

            //double epsf = 0;
            //double epsx = 0.000001;
            int maxits = 0;
            int info;
            alglib.lsfitstate state;
            alglib.lsfitreport rep;
            //double diffstep = 0.0001;

            double epsf = 0;
            double epsx = 0;//0.000000001;
            double diffstep = 0.0001;

            //
            // WITH A BAD INITIAL GUESS and WITHOUT the Restarts done in LogisticFit5, we get high error on all parameters
            //
            alglib.lsfitcreatef(x_conc_M, y_RU, c_5PL, diffstep, out state);
            alglib.lsfitsetcond(state, epsf, epsx, maxits);

            alglib.lsfitsetbc(state, c_5PL_low, c_5PL_High);

            alglib.lsfitfit(state, f_5PL, null, null);
            alglib.lsfitresults(state, out info, out c_5PL, out rep);

            //Make a Couple of predictions with the BAD PARAMETERS
            double fTestBad = 0; //exp: ~5.2

            f_5PL(c_5PL, new double[] { 1.10E-05 }, ref fTestBad, null);


            //Now Get the LogisticFit Sol'n and check it's error 
            double[] c_best = Fit_5PL(x_conc_M, y_RU);
            double[] c_bestConfirmed = Fit_5PL(x_conc_M, y_RU);


            alglib.lsfitcreatef(x_conc_M, y_RU, c_best, diffstep, out state);
            alglib.lsfitsetcond(state, epsf, epsx, maxits);

            alglib.lsfitsetbc(state, c_5PL_low, c_5PL_High);

            alglib.lsfitfit(state, f_5PL, null, null);
            alglib.lsfitresults(state, out info, out c_bestConfirmed, out rep);


            //Get the Parameter Error values
            double[] c_error = rep.errpar;

            //5PL Solution from DoubleReferenceGUI (a,b,c,d,g,r2) = -2.22332885169764	0.629493994723955	9.07171520030884E-05	13.4940430048216	2.73435657052643	0.999601385180452
            double[] c_FromOtherProg = new double[] { -2.22332885169764, 0.629493994723955, 9.07171520030884E-05, 13.4940430048216, 2.73435657052643 };
            double r2_FromOtherProgr = 0.999601385180452;

            for (int i = 0; i < c_best.Length; i++)
            {
                System.Diagnostics.Debug.Assert(Math.Abs(c_best[i] - c_bestConfirmed[i]) < 1e-7, "logisticfit sol'n should be unchanged after extra lsfit");
                System.Diagnostics.Debug.Assert(Math.Abs(c_best[i] - c_FromOtherProg[i]) < 1e-7, "logisticfit sol'n should match other prog");
            }
            System.Diagnostics.Debug.Assert(Math.Abs(rep.r2 - r2_FromOtherProgr) < 1e-7, "logisticfit sol'n should match other prog");


            System.Console.WriteLine("{0}", info); // EXPECTED: 2
            System.Console.WriteLine("{0}", alglib.ap.format(c_5PL, 1)); // EXPECTED: [1.5]
        }


        /// <summary>
        /// Use the LogisticFit helper to handle restarts and find optimum solution, that we then hand to LSFitFit to get error estimates
        /// </summary>
        /// <param name="xMatrix"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        static double[] Fit_5PL(double[,] xMatrix, double[] y)
        {

            System.Diagnostics.Debug.Assert(xMatrix.GetLowerBound(0) == 0, "we're going to create vector from all items in matrix, assuming it is a 'n x 1' matrix");
            double[] x = (from double xv in xMatrix
                          select xv).ToArray();

            alglib.lsfitreport fit;
            double a, b, c, d, g;
            alglib.logisticfit5(x, y, x.Length, out a, out b, out c, out d, out g, out fit);

            return new double[] { a, b, c, d, g };
        }

    }
    

}
