using System;
using System.Runtime.InteropServices;
using FFTWSharp;

namespace FFTWSharp_test
{
    public class FFTWtest
    {
       
        ~FFTWtest()
        {
        }

        static void Main(string[] args)
        {
            const int sampleSize = 1024;

            double[] din = new double[sampleSize * 2];
            double[] dout = new double[sampleSize * 2];

            din[0] = 1;

            fftw_complexarray mdin = new fftw_complexarray(din);
            fftw_complexarray mdout = new fftw_complexarray(dout);

            fftw_plan plan = fftw_plan.dft_1d(sampleSize, mdin, mdout, fftw_direction.Forward, fftw_flags.Estimate);
            plan.Execute();

            plan = fftw_plan.dft_1d(sampleSize, mdout, mdin, fftw_direction.Forward, fftw_flags.Estimate);
            plan.Execute();

            System.Numerics.Complex[] o = mdin.GetData_Complex();

        }
    }
}