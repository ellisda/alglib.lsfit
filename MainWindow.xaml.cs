using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace AlgLib.LsFit
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            //BiOptix.AlgLib.LsFit_Experiment1.Run_lsfit_d_nlf();

            //BiOptix.AlgLib.LsFit_Experiment2.Run_lsfit_5PL();


            //AlgLib.LsFit.ODE.Example_ODE.DoSomething();

            //AlgLib.LsFit.ODE.AssocRateEqn.DoSomething();

            AlgLib.LsFit.ODE.MassTransportThreeRateEqn.DoSomething_TakeForever();

            BiOptix.AlgLib.LsFit_TimeSeriesKinetics.Run_lsfit_timeSeries();

            InitializeComponent();
        }
    }
}
