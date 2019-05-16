using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using FFTWSharp;

namespace PhaseVocoder
{
    /* Short-Time Fourier Transform class 
     * Author: Geir K.Nilsen, geir.kjetil.nilsen@gmail.com, 2014-2017
     */

    public class STFT
    {
        public int hopSize;
        public int winLength;
        private Complex[] win;

        public STFT(int winLength)
        {
            setWinLength(winLength);
            this.hopSize = (int)((double)winLength / 4.0);
        }

        public void setWinLength(int winLength)
        {
            this.winLength = winLength;

            /* Periodic Hanning window */

            win = new Complex[winLength];

            for(int sample = 0; sample < winLength; sample++)
            {
                win[sample] = new Complex(0.5 * (1.0 - Math.Cos(2.0*Math.PI*(double)sample/(double)winLength)), 0.0);
            }                      
        }

        public TFR forward(double[] input_signal)
        {
            int ns = input_signal.Length;
            int nw = winLength;

            fftw_complexarray mdin = new fftw_complexarray(winLength); 
            fftw_complexarray mdout = new fftw_complexarray(winLength);
            fftw_plan plan;
            
            int nf = (int)Math.Floor((double) (ns - nw) / (double)hopSize);

            Complex[][] fframe = new Complex[nf + 1][];
            for (int sample = 0; sample < fframe.Length; sample++)
                fframe[sample] = new Complex[nw];

            Complex[] buf = new Complex[nw];
            Complex factor = new Complex(2.0/3.0, 0.0);

            int hopInd = 0;
            for(int hop = 0; hop <= nf*hopSize; hop+=hopSize)
            {
                for (int sample = 0; sample < nw; sample++)
                {
                    buf[sample] = new Complex(input_signal[hop + sample], 0.0) * factor * win[sample];
                }

                mdout.SetData(fframe[hopInd]);
                mdin.SetData(buf);
                plan = fftw_plan.dft_1d(winLength, mdin, mdout, fftw_direction.Forward, fftw_flags.Estimate);
                plan.Execute();

                fframe[hopInd++] = mdout.GetData_Complex();
            }

            return new TFR(fframe);
        }

        public double[] inverse(TFR tfr, ref double[] output_signal_prev, bool last)
        {
            int nf = tfr.getTFRAsDoubleArray().Length - 1;
            int nw = tfr.getTFRAsDoubleArray()[0].Length;
            int ns = nf * hopSize + nw;
            double[] output_signal_tmp = new double[ns];
            double[] output_signal = new double[output_signal_tmp.Length - 3 * hopSize];

            fftw_complexarray mdin = new fftw_complexarray(nw);
            fftw_complexarray mdout = new fftw_complexarray(nw);
            fftw_plan plan;

            if (output_signal_prev == null)
                output_signal_prev = new double[3 * hopSize];
            for (int sample = 0; sample < hopSize * 3; sample++)
                output_signal_tmp[sample] = output_signal_prev[sample];

            int hopInd = 0;
            Complex[] buf = new Complex[nw];

            try
            {
                for (int hop = 0; hop <= nf * hopSize; hop += hopSize)
                {
                    mdout.SetData(buf);
                    mdin.SetData(tfr.getTFRAsDoubleArray()[hopInd++]);

                    plan = fftw_plan.dft_1d(winLength, mdin, mdout, fftw_direction.Backward, fftw_flags.Estimate);
                    plan.Execute();

                    buf = mdout.GetData_Complex();

                    for (int sample = 0; sample < nw; sample++)
                        output_signal_tmp[sample + hop] = output_signal_tmp[sample + hop] + (buf[sample].Real / nw) * win[sample].Real;
                    
                }

                for (int sample = 0; sample < hopSize * 3; sample++)
                    output_signal_prev[sample] = output_signal_tmp[output_signal_tmp.Length - hopSize * 3 + sample];

                for (int sample = 0; sample < output_signal.Length; sample++)
                    output_signal[sample] = output_signal_tmp[sample];

                if (last)
                    return output_signal_tmp;
                else
                    return output_signal;
            }
            catch (Exception)
            {
                if (last)
                    return output_signal_tmp;
                else
                    return output_signal;
            }
        }

        public double[] inverse(TFR tfr)
        {
            int nf = tfr.getTFRAsDoubleArray().Length - 1;
            int nw = tfr.getTFRAsDoubleArray()[0].Length;
            int ns = nf * hopSize + nw;
            double[] output = new double[ns];

            fftw_complexarray mdin = new fftw_complexarray(nw); 
            fftw_complexarray mdout = new fftw_complexarray(nw);
            fftw_plan plan;

            int hopInd = 0;
            Complex[] buf = new Complex[nw];

            try
            {
                for (int hop = 0; hop <= nf * hopSize; hop += hopSize)
                {
                    mdout.SetData(buf);
                    mdin.SetData(tfr.getTFRAsDoubleArray()[hopInd++]);

                    plan = fftw_plan.dft_1d(winLength, mdin, mdout, fftw_direction.Backward, fftw_flags.Estimate);
                    plan.Execute();

                    buf = mdout.GetData_Complex();

                    for (int sample = 0; sample < nw; sample++)
                        output[sample + hop] = output[sample + hop] + (buf[sample].Real / nw) * win[sample].Real;
                }

                return output;
            }
            catch (Exception)
            {
                return output;
            }
        }
    }
}
