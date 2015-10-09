using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BiOptix.AlgLib
{
    /// <summary>
    /// Meant to replicate Kevin's work in matlab with lsqnonlin and nlparci
    /// REF: https://bioptix.sharepoint.com/A4softwareupdates/A4SoftwareDevelopment/_layouts/OneNote.aspx?id=%2fA4softwareupdates%2fA4SoftwareDevelopment%2fDocuments%2fA4%20Data%20Artifacts%2c%20Test%20Plans%2c%20and%20Resulting%20Software%20Changes&wd=target%28BiOptix%20Analysis%20Tool.one%7c1339185A-20B7-4CD1-BB19-F1EAF84D1C1F%2fMatlab%20-%20full%20routine%20-%20Analytical%20and%20ODE%7cAD99BFB8-1E6F-4E0C-8677-DE820A083430%2f%29
    /// </summary>
    class LsFit_TimeSeriesKinetics
    {

        /// <summary>
        /// I'm tempted to try and add other pieces like t0 and r0 into the x vector, since
        ///   it says "Callbacks for  parameterized  functions,  i.e.  for  functions  which depend on two vectors: P and Q.  Gradient  and Hessian are calculated with respect to P only."
        /// REF: definiaiton of `ndimensional_pfunc` REF: C:\BiOptix\svn2git\test\arapahoe-old\branch-Release_A1\WpfDataBrowser\WpfDataBrowser\AlgLib\ap.cs
        /// </summary>
        /// <param name="c"></param>
        /// <param name="x"></param>
        /// <param name="func"></param>
        /// <param name="obj"></param>
        public static void func(double[] c, double[] x, ref double func, object obj)
        {
            double ka = c[0];
            double kd = c[1];
            double rmax = c[2];

            //iPoint = (int)x[1];
            //iCurve = iPoint / 241;

            //int iCurveX = (int)x[1] / 241;
            //iCurve = iCurveX;
            //System.Diagnostics.Debug.Assert(iCurveX == iCurve, "new old should be same");
            //double conc_m = curveConc_M[iCurve];
            //System.Diagnostics.Debug.WriteLine("i={0}, curve_i={1}, curve_m={2:E3}", iPoint, iCurve, conc_m);


            double conc_m = curveConc_M[iCurve];
            //double conc_m = x[2];

            double t = x[0];
            getConcAndTimeFromSequenceNumber((int)x[0], ref t, ref conc_m);

            const double tR0 = 60; //time at which point we will predict R0 value for beginning of dissociation
            const double tDissocStart = 60; //time at which point we say the dissociation has been happening for 0 sec
            double tDissocDuration = t - tDissocStart;

            //Console.WriteLine("Evaluating iPoint={0:0000} iCurve={1} X_sec={2:000.000}", iPoint, iCurve, t);
            //Console.WriteLine("Evaluating t={3:000.000}, conc_M={4}", iPoint, iCurve, x[0], t, conc_m);

            System.Diagnostics.Debug.Assert(t == x[1], "time are same");
            System.Diagnostics.Debug.Assert(conc_m == x[2], "conc are same");

            if (t <= tR0)
            {
                func = KineticModel_AssocOneToOne(t, conc_m * 1e9, ka, kd, rmax);
                //Store this Model value because it becomes a constant in the dissociation region where we decay from a given y value

                R0_Model = func;
            }
            else if (t >= tDissocStart)
            {
                func = KineticModel_DissocOneToOne(tDissocDuration, R0_Model, kd);
            }

            //DESIGN REVIEW - Keep state to know when we switch to next concentration
            iPoint++;
            if(iPoint == totalPoints)
            {
                iPoint = 0;
                iCurve = 0;
            }

            else if (iPoint > 0 && (iPoint % 241) == 0)
                iCurve++;
        }

        private static int iCurve = 0;
        private static int iPoint = 0;
        private static double[] curveConc_M = null;
        private static int totalPoints;


        private static double R0_Model = 0;


        /// <summary>
        /// 
        /// </summary>
        /// <param name="c_params"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public double computeModel(double[] c_params, double[] x)
        {
            double ka = c_params[0];
            double kd = c_params[1];
            double rmax = c_params[2];

            //DESIGN NOTE: Since our transitions are not instantaneous, we may need to predict the R0 value before valve switch and begin counting dissociation duration AFTER valve switch
            
            //DESIGN REVIEW - For now, we'll set the two times below to be equal (implies instantaneous valve switch)
            const double tR0 = 60; //time at which point we will predict R0 value for beginning of dissociation
            const double tInjStart = 60; //time at which point we say the dissociation has been happening for 0 sec

            double r0 = 0;// Design review - how to hold on to this value except for stateful parameter? (static feels like global)

            return double.NaN;
        }


        /// <summary>
        /// Repeating Kevin's Matlab experiments with AlgLib minlm and lsFitFit to get standard error for fitted parameters
        /// </summary>
        public static void Run_lsfit_timeSeries()
        {
            //double[,] ff = new double[,] { _conc1RU, _conc2RU, _conc3RU, _conc4RU };


            
            ////Data from ABS (data that DoubleReferenceGUI starts with, ch3, DMSO-disabled, standard blanking, Bound Avg [10,30]
            //double[,] x_conc_M = new double[,] { { 1.23E-06 }, { 3.70E-06 }, { 1.10E-05 }, { 1.10E-05 }, { 3.30E-05 }, { 3.30E-05 }, { 1.00E-04 } };
            //double[] y_RU = new double[] { 0.32110643, 2.335124169, 5.0964102, 5.361168514, 8.562463415, 8.5834102, 11.32491796 };

            //double[] c_5PL = new double[] { -10, 0.5, 0.003, 20, 1.0 };

            ////double[] c_5PL = new double[] /*VALUES FROM Existin Algorithm*/ { -2.22332885169764, 0.629493994723955, 9.07171520030884E-05, 13.4940430048216, 2.73435657052643 };

            ////Apparrently the 5PL calc function requires c to be positive
            //double[] c_5PL_low = new double[] { double.NegativeInfinity, double.NegativeInfinity, 1e-14, double.NegativeInfinity, double.NegativeInfinity };
            //double[] c_5PL_High = new double[] { double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity };

            double[] c_OneToOne = new double[] { 0.000005, 0.000005, 100 };

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


            var x_Time = getXMatrix();
            var y_RU = getYVector();
            curveConc_M = _concentrations_M;
            totalPoints = y_RU.Length;


            var timer = System.Diagnostics.Stopwatch.StartNew();

            //
            // WITH A BAD INITIAL GUESS and WITHOUT the Restarts done in LogisticFit5, we get high error on all parameters
            //
            alglib.lsfitcreatef(x_Time, y_RU, c_OneToOne, diffstep, out state);
            alglib.lsfitsetcond(state, epsf, epsx, maxits);

            //alglib.lsfitsetbc(state, c_5PL_low, c_5PL_High);

            alglib.lsfitfit(state, func, null, null);
            alglib.lsfitresults(state, out info, out c_OneToOne, out rep);

            string c_string = string.Format("{0:0.####E+0}, {1:0.####E+0}, {2:0.####E+0}", c_OneToOne[0], c_OneToOne[1], c_OneToOne[2]);
            string cErr_string = string.Format("{0:0.####E+0}, {1:0.####E+0}, {2:0.####E+0}", rep.errpar[0], rep.errpar[1], rep.errpar[2]);
            Console.WriteLine("{0:00} iterations, {1:ss\\:fffffff} elapsed.", rep.iterationscount, timer.Elapsed);
            Console.WriteLine("\tC\t={0}", c_string);
            Console.WriteLine("\tErr\t={0}", cErr_string);

            foreach (var cpi in Enumerable.Range(0, rep.covpar.GetUpperBound(0) + 1))
            {
                foreach (var cpj in Enumerable.Range(0, rep.covpar.GetUpperBound(1) + 1))
                    Console.Write("({1},{2}) {0:0.####E+0}", rep.covpar[cpi, cpj], cpi, cpj);
                Console.WriteLine();
            }
            
            //Get Two-Tail T Value (StudentDistribution) for 95% Confidence Interval
            //DESGIN NOTE: It seems that we're getting back one tail from the function below, so for two tail we have to ask for 97.5% condifdence interval (so we can use it as both plus and minus)
            double p = 0.975;
            int df = y_RU.Length - c_OneToOne.Length; //number of points minus number of parameters (ref: http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?reg_how_standard_errors_are_comput.htm)
            var t = alglib.invstudenttdistribution(df, p);

            //Calculate Confidence Interval
            Console.WriteLine("\t TInverse({0}, {1:0.000}) = {2:0.####}", df, p, t);
            var confIntervals = from i in Enumerable.Range(0, c_OneToOne.Length)
                                select new double[]{c_OneToOne[i], c_OneToOne[i] - t*rep.errpar[i], c_OneToOne[i] + t*rep.errpar[i]};

            //
            foreach (var conf in confIntervals)
                Console.WriteLine("\t {0:0.####E+0} - 95% Conf Interval [{1:0.####E+0}, {2:0.####E+0}]", conf[0], conf[1], conf[2]);
        }


        /// <summary>
        /// set the conc and time fields based on the sequence number of the data point
        /// </summary>
        /// <param name="iPoint"></param>
        /// <param name="t"></param>
        /// <param name="conc_M"></param>
        private static void getConcAndTimeFromSequenceNumber(int iPoint, ref double t, ref double conc_M)
        {
            iPoint /= 3;
            int iCurve = iPoint / 241;
            conc_M = _concentrations_M[iCurve];
            t = (iPoint % 241) / 2.0;
        }

        private static double[,] getXMatrix()
        {
            var xVec = getXVector();
            double[,] x = new double[xVec.Length, 3];
            for(int i = 0; i < xVec.Length; i++)
            {
                //x[i, 0] = xVec[i];
                x[i, 0] = 3 * i; //Store Sequence Number
                x[i, 1] = xVec[i];
                x[i, 2] = _concentrations_M[i / 241];
            }

            //DESIGN REVIEW - Maybe I can add extra timing and concentration metadata as other x-dimensions (would be cleaner)
            return x;
        }

        /// <summary>
        /// helper to get the time as a Nx1 matrix of time data (where x is our single dimen
        /// </summary>
        /// <returns></returns>
        private static double[] getXVector()
        {
            //t = [0:0.5:120]';               %time series
            //tstep=60;                       %at time 60, the concentration goes to zero

            double[] time = Enumerable.Range(0, 241).Select(v => (double)v / 2).ToArray();

            double[] allTime = time.Concat(time).Concat(time).Concat(time).Concat(time).ToArray();

            return allTime;
        }

        /// <summary>
        /// helper to get the RU data in a single vector
        /// </summary>
        /// <returns></returns>
        private static double[] getYVector()
        {
            var allRU = _conc1_RU.Concat(_conc2_RU).Concat(_conc3_RU).Concat(_conc4_RU).Concat(_conc5_RU).ToArray();

            //double[,] y_RU = new double[1, allRU.Length];
            //for (int i = 0; i < allRU.Length; i++)
            //    y_RU[0, i] = allRU[i];
            //
            //return y_RU;

            return allRU;
        }

        private static double[] _concentrations_M = new double[] { 1.00E-04, 1.10E-05, 3.70E-06, 3.30E-05, 1.23E-06};
        private static double[] _conc1_RU = new double[] { -0.0136, 6.3584, 16.2994, 21.8584, 25.4714, 27.5324, 28.7144, 29.8434, 30.3324, 30.5034, 30.9154, 31.4124, 31.5294, 31.7144, 31.9854, 31.8514, 31.7584, 31.8774, 32.0564, 32.2614, 32.2014, 32.0034, 32.0224, 32.2054, 32.2354, 31.9644, 31.8704, 31.8934, 32.1294, 32.3034, 32.0514, 32.1634, 32.1834, 32.1934, 32.0114, 31.6184, 31.9274, 32.1864, 32.0474, 31.9494, 32.0134, 31.9724, 31.9734, 31.9004, 32.0534, 32.1154, 31.9314, 31.8954, 31.9624, 32.1244, 32.2184, 31.9914, 31.7904, 32.0384, 31.9044, 31.7244, 31.9584, 31.8104, 32.0354, 31.9574, 31.7934, 32.0714, 32.1424, 32.0064, 31.8634, 31.9084, 31.8424, 31.7004, 31.4694, 31.5634, 31.5514, 31.6684, 31.7574, 31.4724, 31.5224, 31.4874, 31.6244, 31.6124, 31.5704, 31.7324, 31.6654, 31.7464, 31.8034, 31.7014, 31.6604, 31.6714, 31.7024, 31.6744, 31.7094, 31.4964, 31.4834, 31.8204, 32.0334, 31.8964, 31.7894, 31.8254, 31.8494, 31.7414, 31.7174, 31.5904, 31.6484, 31.8264, 31.5544, 31.5704, 31.6874, 31.5964, 31.5224, 31.5274, 31.6484, 31.6424, 31.5094, 31.5404, 31.6504, 31.8614, 31.7164, 31.6834, 31.6124, 31.4454, 31.0434, 30.8004, 29.6254, 28.3014, 27.2174, 25.9134, 24.1684, 22.9674, 21.8414, 20.6004, 19.4864, 18.7624, 18.0604, 16.9014, 15.8254, 14.9804, 14.1554, 13.7344, 13.1514, 12.1064, 11.4044, 11.1044, 10.6504, 10.0294, 9.5604, 9.2294, 8.8724, 8.5134, 8.1054, 7.4824, 7.1244, 7.1554, 6.7144, 6.3814, 6.0734, 5.7694, 5.5634, 5.3174, 5.0074, 4.7434, 4.8594, 4.8634, 4.5184, 4.1304, 3.9244, 3.7284, 3.4754, 3.3964, 3.1904, 2.9144, 2.8274, 2.8804, 2.8334, 2.5264, 2.2774, 2.1444, 2.2794, 2.3094, 2.2084, 1.9924, 1.8424, 1.8854, 1.7404, 1.5594, 1.5244, 1.4294, 1.2784, 1.2754, 1.3404, 1.2074, 1.2194, 1.0614, 0.9604, 1.0434, 0.8704, 0.9734, 0.9064, 0.8684, 0.9804, 1.0754, 0.8654, 0.5624, 0.6504, 0.6654, 0.7784, 0.9064, 0.7564, 0.5664, 0.5584, 0.6144, 0.6184, 0.6174, 0.5804, 0.4424, 0.4244, 0.3054, 0.1364, 0.2144, 0.2264, 0.3254, 0.4764, 0.5254, 0.3954, 0.4354, 0.7294, 0.5914, 0.4174, 0.4844, 0.6834, 0.8144, 0.6654, 0.5934, 0.6014, 0.5744, 0.5094, 0.7164, 0.7244, 0.5394, 0.4744, 0.5844, 0.6084, 0.5274, 0.5184 };
        private static double[] _conc2_RU = new double[] { -0.1871, 0.9619, 3.2249, 4.0889, 5.0879, 6.0689, 7.0659, 8.3759, 9.1129, 10.1029, 10.9059, 11.6169, 12.3609, 12.8839, 13.6609, 14.2579, 14.6419, 14.9669, 15.1829, 15.7379, 16.1099, 16.1549, 16.4019, 16.6609, 16.8769, 16.8289, 16.8469, 16.9359, 17.2319, 17.5479, 17.4269, 17.4429, 17.5539, 17.8379, 17.6479, 17.3319, 17.7109, 18.0649, 17.9799, 17.8059, 17.7669, 17.8369, 17.8839, 18.0099, 18.2079, 18.1479, 17.9889, 18.0019, 18.2839, 18.4079, 18.0219, 17.9829, 18.0839, 18.3059, 18.4079, 18.1439, 18.1109, 18.1079, 18.3139, 18.0919, 18.0169, 18.3319, 18.6729, 18.5349, 18.2679, 18.2929, 18.3719, 18.3769, 18.2709, 18.2919, 18.4649, 18.5199, 18.2119, 18.2299, 18.4669, 18.1919, 18.1429, 18.0809, 18.1219, 18.2559, 18.2549, 18.4589, 18.5249, 18.5589, 18.5189, 18.5649, 18.6079, 18.4209, 18.4159, 18.4749, 18.3959, 18.3629, 18.4329, 18.3209, 18.3979, 18.4809, 18.3549, 18.3749, 18.3929, 18.4199, 18.4939, 18.1979, 18.1779, 18.3719, 18.5089, 18.4909, 18.3569, 18.2489, 18.3619, 18.4409, 18.2379, 18.3029, 18.4459, 18.5269, 18.3709, 18.3119, 18.4549, 18.3809, 18.2089, 18.0919, 16.9989, 16.4789, 15.7349, 15.0839, 14.1029, 13.4819, 13.0219, 12.3019, 11.7389, 11.3239, 10.8449, 10.3229, 9.6409, 8.9169, 8.5059, 8.1359, 7.6159, 7.1589, 6.8669, 6.7519, 6.6769, 6.4149, 6.0809, 5.9649, 5.7869, 5.4279, 5.2289, 5.0169, 4.6139, 4.4539, 4.4529, 4.4799, 4.0879, 3.7179, 3.7389, 3.6649, 3.3939, 3.2439, 3.2419, 3.1449, 2.9969, 2.8649, 2.7559, 2.5049, 2.4039, 2.4749, 2.2339, 2.0759, 2.0979, 2.2829, 2.2559, 1.8439, 1.7819, 1.7899, 1.7439, 1.6559, 1.6519, 1.6829, 1.5479, 1.5389, 1.5149, 1.3569, 1.2149, 1.2359, 1.2119, 1.0619, 1.1089, 1.1429, 1.1159, 1.1919, 1.1519, 0.9439, 0.8879, 1.1479, 1.0339, 0.7439, 0.6899, 0.8909, 0.9389, 0.6579, 0.4789, 0.6299, 0.8889, 0.8319, 0.7039, 0.8179, 0.9339, 1.0269, 0.9949, 0.8229, 0.7149, 0.7779, 0.6539, 0.5979, 0.8129, 0.6469, 0.5779, 0.5299, 0.4399, 0.5309, 0.6339, 0.6689, 0.7379, 0.8029, 0.7769, 0.6579, 0.7479, 0.8709, 0.9209, 0.8619, 0.7199, 0.6779, 0.8309, 1.0549, 0.6509, 0.5689, 0.7779, 0.8529, 0.8019, 0.7969, 0.6889 };
        private static double[] _conc3_RU = new double[] { 0.946, 1.216, 1.736, 2.153, 3.083, 3.899, 4.329, 4.737, 5.075, 5.731, 5.828, 5.997, 6.441, 6.674, 7.07, 7.559, 7.755, 7.742, 7.941, 8.42, 8.596, 8.608, 8.892, 9.004, 8.969, 8.893, 9.123, 9.342, 9.185, 9.378, 9.696, 9.904, 9.929, 9.974, 9.782, 9.718, 10.064, 10.28, 10.141, 10.016, 10.002, 10.015, 10.114, 10.174, 10.166, 10.339, 10.655, 10.541, 10.549, 10.817, 10.578, 10.629, 10.703, 10.861, 10.974, 10.74, 10.863, 10.746, 10.817, 10.643, 10.392, 10.521, 10.905, 10.997, 10.809, 10.683, 10.635, 10.585, 10.472, 10.52, 10.539, 10.7, 10.669, 10.625, 10.648, 10.754, 10.845, 10.609, 10.692, 10.787, 10.649, 10.653, 10.789, 10.923, 10.915, 11.087, 11.202, 10.905, 10.755, 10.796, 10.941, 10.881, 10.845, 10.949, 10.938, 11.011, 11.109, 10.907, 10.772, 10.802, 11.098, 11.078, 10.865, 10.977, 11.046, 11.302, 11.14, 11.003, 11.026, 11.047, 11.119, 10.95, 11.114, 11.436, 11.185, 11.044, 11.143, 11.297, 10.952, 10.254, 9.577, 8.5182, 8.6692, 8.5682, 7.6502, 7.2432, 7.0252, 6.5732, 6.2602, 6.2422, 6.1102, 5.7282, 5.1662, 4.7242, 4.7022, 4.7722, 4.5352, 4.0272, 3.8702, 3.9302, 3.8592, 3.5322, 3.3212, 3.3522, 3.2022, 3.0802, 3.1582, 3.0572, 2.6272, 2.5362, 2.4852, 2.4492, 2.4412, 2.3472, 2.1592, 2.1302, 1.9602, 1.8792, 2.2642, 2.2742, 1.9472, 1.8802, 1.7372, 1.5032, 1.5602, 1.4322, 1.3202, 1.2862, 1.3182, 1.5412, 1.5352, 1.1472, 1.0882, 0.9102, 0.8742, 1.0632, 1.1542, 1.1242, 1.0472, 1.2102, 1.2062, 0.9482, 1.0472, 1.1922, 0.8822, 0.7752, 1.1432, 1.2732, 1.0902, 1.1022, 1.0992, 1.1012, 0.9922, 0.8762, 0.9252, 0.9392, 0.9012, 0.9302, 0.7032, 0.4202, 0.5242, 0.7722, 0.9512, 0.8522, 0.7952, 0.7572, 0.8082, 0.6912, 0.5162, 0.6802, 0.6262, 0.3482, 0.4452, 0.5972, 0.6422, 0.5202, 0.3672, 0.4612, 0.5692, 0.6452, 0.6212, 0.6942, 0.8372, 0.7692, 0.8462, 0.7462, 0.7882, 0.8832, 0.8272, 0.7462, 0.8062, 0.8622, 0.6292, 0.8142, 0.7082, 0.4142, 0.6792, 0.9802, 0.8922, 0.7002, 0.6662 };
        private static double[] _conc4_RU = new double[] { 0.0831, 1.7871, 4.5041, 6.8541, 10.4541, 13.1391, 14.9221, 17.1481, 18.8941, 20.1871, 21.1611, 21.8461, 22.7881, 23.5231, 24.0921, 24.4521, 24.7561, 24.9151, 25.1241, 25.7121, 25.8381, 25.7561, 25.8821, 26.1491, 26.2811, 26.0201, 26.0781, 26.1461, 26.3621, 26.6521, 26.5501, 26.5431, 26.6051, 26.5841, 26.3131, 26.1041, 26.4221, 26.6311, 26.6121, 26.6251, 26.4231, 26.2511, 26.3261, 26.4361, 26.4471, 26.5411, 26.6661, 26.3791, 26.4641, 26.8741, 26.9051, 26.6481, 26.3261, 26.5481, 26.6191, 26.4661, 26.5071, 26.4161, 26.6301, 26.4931, 26.2451, 26.5221, 26.7221, 26.6691, 26.5321, 26.5161, 26.4121, 26.5321, 26.5221, 26.4411, 26.4361, 26.6481, 26.6881, 26.7251, 26.8761, 26.4491, 26.4401, 26.4061, 26.4711, 26.7271, 26.4591, 26.4271, 26.5081, 26.5391, 26.4931, 26.6081, 26.6561, 26.5821, 26.6951, 26.6971, 26.6271, 26.7631, 26.8271, 26.4241, 26.4701, 26.8341, 26.6281, 26.6071, 26.7141, 26.5531, 26.7541, 26.7551, 26.5841, 26.8801, 26.9301, 26.9811, 26.9031, 26.5611, 26.7921, 26.7781, 26.4111, 26.4741, 26.6851, 27.0601, 26.9511, 26.7331, 26.7401, 26.7391, 26.3981, 26.1191, 25.2061, 24.1096, 22.9856, 22.0316, 20.6946, 19.7866, 19.1096, 18.1766, 17.3556, 16.5086, 15.6026, 14.7876, 13.9066, 13.1756, 12.5576, 12.0586, 11.5926, 10.8706, 10.2986, 10.0646, 9.8396, 9.2216, 8.6756, 8.4006, 7.9856, 7.5016, 7.2966, 7.1236, 6.7366, 6.5056, 6.3366, 6.1836, 5.7406, 5.5466, 5.4076, 5.1056, 4.9526, 4.7596, 4.7156, 4.5476, 4.1466, 3.9956, 3.8046, 3.4606, 3.2586, 3.2456, 3.1716, 2.8096, 2.7696, 2.9636, 3.0186, 2.9246, 2.7076, 2.4136, 2.3136, 2.2376, 2.3216, 2.2636, 1.9656, 2.0166, 2.1596, 2.0426, 1.9856, 1.8286, 1.6096, 1.7626, 1.8986, 1.7926, 1.7416, 1.7536, 1.6216, 1.3396, 1.2776, 1.3876, 1.1916, 1.1566, 1.1126, 1.2596, 1.4186, 1.1656, 0.9446, 1.1916, 1.4746, 1.3376, 1.2086, 1.2886, 1.3906, 1.2216, 1.0626, 1.2456, 1.2226, 1.1866, 1.0636, 0.9106, 0.9156, 0.5036, 0.4656, 0.6096, 0.5616, 0.8346, 0.8766, 0.7646, 0.9646, 0.9426, 0.7106, 0.7106, 0.8726, 0.9836, 1.0116, 0.8396, 0.8016, 0.7326, 0.7646, 1.0836, 0.7926, 0.7826, 0.9046, 1.0176, 1.1456, 1.1116, 0.8706 };
        private static double[] _conc5_RU = new double[] { 2.373, 2.583, 2.207, 0.939, 0.722, 0.894, 1.091, 1.238, 1.202, 1.52, 1.698, 1.778, 2.178, 2.523, 2.909, 3.09, 3.033, 3.105, 3.256, 3.633, 3.808, 3.596, 3.526, 3.896, 4.191, 3.95, 3.831, 3.808, 3.842, 4.047, 4.087, 4.136, 4.247, 4.399, 4.201, 4.041, 4.34, 4.536, 4.423, 4.296, 4.209, 4.315, 4.372, 4.602, 4.773, 4.633, 4.656, 4.411, 4.449, 4.607, 4.619, 4.652, 4.499, 4.623, 4.799, 4.883, 4.984, 4.795, 4.875, 4.611, 4.476, 4.642, 4.866, 5.091, 5.017, 5.092, 4.912, 4.905, 4.96, 4.762, 4.661, 4.797, 4.958, 4.953, 4.82, 4.624, 4.705, 4.51, 4.603, 4.986, 4.928, 5.099, 5.192, 5.086, 5.064, 5.125, 5.157, 4.935, 4.91, 4.873, 4.768, 4.902, 5.161, 5.177, 5.01, 4.919, 4.935, 4.941, 4.91, 4.949, 5.134, 5.132, 5.042, 5.121, 5.198, 5.263, 5.244, 5.165, 5.281, 5.233, 5.067, 5.048, 5.102, 5.254, 4.81, 4.637, 4.938, 5.144, 4.542, 3.65, 4.139, 2.9454, 3.7424, 3.9584, 3.5494, 3.3644, 3.0804, 2.9254, 2.9284, 2.9854, 2.9654, 2.8254, 2.5784, 2.1914, 2.0604, 2.2074, 2.2114, 1.9974, 1.9594, 2.0854, 1.9704, 1.6744, 1.5984, 1.7224, 1.8624, 1.9114, 1.8744, 1.6764, 1.5784, 1.7564, 1.6614, 1.3774, 1.4784, 1.5174, 1.4164, 1.4544, 1.4694, 1.5444, 1.7024, 1.8084, 1.5404, 1.3094, 1.3564, 1.3374, 1.1794, 1.1584, 1.1244, 1.1194, 1.2104, 1.1554, 1.0594, 1.1154, 1.0224, 0.9584, 1.1874, 1.0044, 1.1064, 1.0864, 0.8604, 0.9294, 0.9384, 0.9204, 0.9784, 0.8224, 0.6104, 0.5804, 0.7584, 0.6934, 0.7114, 0.8674, 0.9434, 0.9414, 0.6884, 0.8364, 0.8064, 0.8184, 0.7854, 0.7504, 0.8524, 0.6944, 0.7894, 0.8904, 0.7254, 0.6544, 0.7414, 0.8804, 1.0614, 0.9034, 0.7334, 0.8214, 0.9394, 0.9204, 0.7824, 0.7574, 0.7344, 0.5274, 0.4324, 0.7374, 0.8544, 0.6704, 0.5694, 0.7204, 1.0444, 1.0324, 0.8194, 0.7854, 1.0804, 1.2214, 1.1094, 1.0204, 1.0934, 1.0794, 1.2184, 1.4424, 0.9564, 0.8054, 1.0274, 1.1074, 1.0484, 1.0104, 0.9624 };

        //private double tInv(int k, double p)
        //{

        //    //http://forum.alglib.net/viewtopic.php?f=2&t=149
        //    return alglib.invstudenttdistribution(k, p);

        //    return -1* (alglib.invstudenttdistribution()
        //}


        #region helper methods from Old Data Viewer
        private static double KineticModel_AssocOneToOne(double t, double sampleConcentrationNM, double kOn, double kOff, double S)
        {
            //double RMax = fc.Kon * (fcDetails.SampleConcentrationNM * 1000000000) * .......
            double Rc = kOn * sampleConcentrationNM * Math.Pow(10, -9) * S
                / (kOn * sampleConcentrationNM * Math.Pow(10, -9) + kOff);
            double Kob = kOn * sampleConcentrationNM * Math.Pow(10, -9) + kOff;
            //Func<double> Rc = () => fc.; //(Func<double>)delegate(sender){return 5.0;};

            //for (double i = 0.5; i < fcDetails.AssociationEnd - fcDetails.InjectionStart; i += 0.5)
            return Rc * (1 - Math.Exp(-1 * Kob * t)); //DE[4-21-11] Pg 200 equation:  Y = (Bmax * X) / (X + Kd)
            //kPlot.PlotXYAppend(i, (1 - Math.Exp(-1 * i / 100))); //DE[4-21-11] Pg 200 equation:  Y = (Bmax * X) / (X + Kd)
        }

        private static double KineticModel_DissocOneToOne(double t, double Ri, double kOff)
        {
            return Ri * Math.Exp(-1 * kOff * t);
        }
        #endregion

    }


}
