using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace Resample
{
    /// <summary>
    /// Class for uniformly resampling signals by a rational ratio p/q.
    /// Supports real time resampling.
    /// </summary>
    public class Resampler
    {
        private uint q, p;
        private double[,] filter;
        //Filter length after upsampling (actual filter length divided by p)
        private uint flen;
        private uint nstages;
        private const double defaultAlpha = 5.0;
        private const uint defaultN = 10;

        private double[] state;
        private uint l_c = 0, k_c = 0;

        /// <summary>
        /// Creates a new Resampler object which resamples the signal by a factor p/q.
        /// Default parameters are used for the lowpass filter.
        /// The default filter is a Kaiser-windowed sinc, where the alpha parameter of the 
        /// Kaiser window is 5.0. THe cutoff frequency of the filter (which is applied to the signal
        /// after upsampling by p) is 1/max(p, q) of Nyquist. The filter length is 20*cutoff+1.
        /// </summary>
        /// <param name="p">Numerator of the resampling factor</param>
        /// <param name="q">Denominator of the resampling factor</param>
        public Resampler(uint p, uint q)
        {
            InitKaiser(p, q, defaultN, defaultAlpha, 1.0);
        }

        /// <summary>
        /// Creates a new Resampler object with a custom filter. The cutoff of 
        /// the filter should be at about 1/max(p, q) of the Nyquist frequency.
        /// See also <see cref="Resampler(uint, uint)"/>
        /// </summary>
        /// <param name="p">The numerator of the resampling factor</param>
        /// <param name="q">The denominator of the resampling factor</param>
        /// <param name="filter">FIR filter inpulse response</param>
        public Resampler(uint p, uint q, double[] filter)
        {
            Init(p, q, filter);
        }


        /// <summary>
        /// Creates a Resampler object with default settings except for n, which specifies the length
        /// of the anti-aliasing filter. The actual length of the filter will be 2*n*max(p, q) + 1.
        /// The filter length provides a compromise between lower computational resources and phase delay
        /// (delay only relevant for ResampleContinuous) on one hand versus  better alias suppression and 
        /// less attenuation of high frequencies on the other. The default value of n is 10. 
        /// See also <see cref="Resampler(uint, uint)"/>
        /// </summary>
        /// <param name="p">Numerator of the resampling ratio</param>
        /// <param name="q">Denominator of the resampling ratio</param>
        /// <param name="n">See constructor description</param>
        public Resampler(uint p, uint q, uint n)
        {
            InitKaiser(p, q, n, defaultAlpha, 1.0);
        }

        /// <summary>
        /// Creates a Resampler object with a Kaiser-windowed sinc filter based on the given parameters.
        /// See also <see cref="Resampler(uint, uint)"/>
        /// </summary>
        /// <param name="p">Numerator of the resampling ratio</param>
        /// <param name="q">Denominator of the resampling ratio</param>
        /// <param name="cutoff"> Cutoff frequency of the filter, relative to the Nyquist frequency. Set below 1 for 
        /// better alias suppression and above 1 for less high frequency attenuation.
        /// The default cutoff, used by the other constructors, is 1.0.
        /// </param>
        /// <param name="n">See <see cref="Resampler(uint, uint, uint)"/>. The default is 10.</param>
        /// <param name="alpha">Alpha parameter for the Kaiser window. THe default is 5.0.</param>
        public Resampler(uint p, uint q, double cutoff, uint n, double alpha)
        {
            InitKaiser(p, q, n, alpha, cutoff);
        }

        // Magic happens here
        private void Init(uint p, uint q, double[] f)
        {
            this.q = q;
            this.p = p;
            nstages = 1;
            uint i = q % p;

            //Number of stages is p/GCD(p, q). This is a quick and easy way to find it
            while (i != 0)
            {
                nstages++;
                i = (i + q) % p;
            }
            //Ceiling of f.Length/p
            flen = DivideCeiling((uint)f.Length, p);
            filter = new double[nstages, flen];

            i = 0;
            uint l = 0;
            //Loop of magic
            for (l = 0; l < nstages; l++)
            {
                uint k = 0;
                for (uint j = i; j < f.Length; j += p)
                {
                    filter[l, k] = f[j];
                    k++;
                }
                i = (i + q) % p;
            }
            //Flip it because it simplifies convolution
            flipud(filter);

            //For real-time resampling, we need to store the current state
            state = new double[flen - 1];
        }

       
        /// <summary>
        ///  Initialize using a Kaiser-windowed sinc as the filter
        /// </summary>
        private void InitKaiser(uint p, uint q, uint complexity, double alpha, double relative_cutoff)
        {
            uint maxpq = Math.Max(p, q);
            double standard_cutoff = 1.0 / maxpq;
            double[] f = FilterBuilder.KaiserFilter(relative_cutoff*standard_cutoff, (int)(2 * maxpq * complexity + 1), alpha);
            for (int i = 0; i < f.Length; i++)
                f[i] *= p;
            Init(p, q, f);
        }


        /// <summary>
        /// Corresponds to Matlab's flipud function
        /// </summary>
        private void flipud(double[,] x)
        {
            int m = x.GetLength(0);
            int n = x.GetLength(1); 
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n/2; j++)
                {
                    double temp = x[i, j];
                    x[i, j] = x[i, n - j - 1];
                    x[i, n - j - 1] = temp;
                }
            }
        }


        /// <summary>
        /// Add zeros before and after a signal
        /// </summary>
        private double[] ZeroPad(double[] x, uint before, uint after)
        {
            double[] y = new double[x.Length + before + after];
            for (int i = 0; i < x.Length; i++)
            {
                y[i + before] = x[i];
            }
            return y;
        }


        /// <summary>
        /// This is the main loop. This is where all the magic and all the work happens.
        /// The number of multiply-adds is approximately 2*n*outputLength if p > q and 2*n*x.Length otherwise. The "n" here
        /// is the filter length parameter, 10 for default settings.
        /// If you want to optimize this class, implementing this function in native, non-managed code is the way to go.
        /// </summary>
        /// <param name="x">Input signal</param>
        /// <param name="outputLength">Number of samples of the output signal to be produced</param>
        /// <param name="k">Index of where to start in the upsampled signal on input.
        /// Upon return, k is the index into the upsampled signal where we should
        /// start for the next signal block. </param>
        /// <param name="l">Index of the filter stage where we should start. Upon exit, index
        /// of the filter stage where we should start for the next block of input samples.</param>
        /// <returns>The resampled signal</returns>
        private double[] ResampleLoop(double[] x, uint outputLength, ref uint k, ref uint l)
        {
            double[] y = new double[outputLength];
            for (int i = 0; i < outputLength; i++)
            {
                uint m = k / p;
                for (int j = 0; j < flen; j++)
                {
                    y[i] += x[m + j] * filter[l, j];
                }
                l++;
                l %= nstages;
                k = k + q;
            }
            return y;
        }

        /// <summary>
        /// Resamples the signal x, compensating for delay caused by the lowpass filter.
        /// Do not use this function for real-time processing, use only if the whole signal is 
        /// known at call-time.
        /// </summary>
        /// <param name="x">The input signal</param>
        /// <returns>The resampled output signal</returns>
        public double[] Resample(double[] x)
        {
            uint l = 0;
            uint k = 0;
            uint n = (uint)Math.Floor((double)x.Length * p / q);
            x = ZeroPad(x, flen / 2, flen / 2);
            double[] y = ResampleLoop(x, n, ref k, ref l);
            return y;
        }



        /// <summary>
        /// Puts the samples from the previous state first in the signal,
        /// and stores the last samples from the signal as the next state.
        /// </summary>
        private double[] AppendAndUpdateState(double[] x)
        {
            double[] y = new double[x.Length + flen - 1];
            for (int i = 0; i < flen - 1; i++)
            {
                y[i] = state[i];
            }
            for (int i = 0; i < x.Length; i++)
            {
                y[i + flen - 1] = x[i];
            }
            for (int i = 0; i < flen - 1; i++)
            {
                state[i] = y[x.Length + i];
            }
            return y;
        }


        /// <summary>
        /// Resets the continuous reasmpling filter state.
        /// </summary>
        public void ResetFilter()
        {
            state = new double[flen - 1];
            l_c = 0;
            k_c = 0;
        }
        

        /// <summary>
        /// Returns the ceiling of num/denom
        /// </summary>
        private uint DivideCeiling(uint num, uint denom)
        {
            return (num + denom - 1) / denom;
        }

        /// <summary>
        /// Resampling function which can be called over and over with new contiguous blocks of input signal.
        /// This can be used for real time processing, or for very long signals.
        /// Unlike Resample, this function does not compensate for delay caused by the filtering.
        /// The delay caused by the filter will be n output samples if q > p, 
        /// and n*p/q output samples if p > q. See <see cref="Resampler(uint, uint, uint)"/> for a description of n.
        /// </summary>
        /// <param name="x">Input signal</param>
        /// <returns>Output signal</returns>
        public double[] ResampleContinuous(double[] x)
        {
            uint l = l_c; //l is the stage number in the filter
            uint k = k_c; //k is the start index into the upsampled signal
            uint inLength = (uint)x.Length;
            uint outLength = DivideCeiling(inLength * p - k, q); //Number of output samples depends on previous runs through k
            x = AppendAndUpdateState(x); //Add previous state to the beginning of the signal, and store current state
            double[] y = ResampleLoop(x, outLength, ref k, ref l);
            l_c = l; //Store for next call
            k_c = k - (inLength * p); //Index next time is the current index minus number of upsampled samples in this call
            return y;
        }
    }

    /// <summary>
    /// Static functions for creating linear phase FIR-filters. Doesn't support much at the moment.
    /// </summary>
    public static class FilterBuilder
    {

        /// <summary>
        /// sin(x)/x if x != 0, 1 otherwwise
        /// </summary>
        public static double Sinc(double x)
        {
            if (x == 0.0)
                return 1.0;
            return Math.Sin(x) / x;
        }

        /// <summary>
        /// Creates a rectangular-windowed sinc lowpass filter.
        /// </summary>
        /// <param name="cutoff">Cutoff frequency wrt. the Nyquist frequency</param>
        /// <param name="n">Number of taps</param>
        /// <returns>Filter coefficients</returns>
        public static double[] LowpassSinc(double cutoff, int n)
        {
            double[] ret = new double[n];
            double delay = (double)(n - 1) / 2.0;
            for (int i = 0; i < n; i++)
            {
                double x = ((double)i - delay) * Math.PI * cutoff;
                ret[i] = Sinc(x);
            }
            return ret;
        }

        /// <summary>
        /// Hamming window. What do you want me to say?
        /// </summary>
        /// <param name="n">Number of points</param>
        /// <returns>Window coefficients</returns>
        public static double[] HammingWindow(int n)
        {
            double[] win = new double[n];
            for (int i = 0; i < n; i++)
            {
                win[i] = 0.54 - 0.46 * Math.Cos(2 * Math.PI * i / (n - 1));
            }
            return win;
        }

        /// <summary>
        /// Multiplies the given window by a sinc lowpass filter
        /// to produce a windowed sinc filter.
        /// </summary>
        /// <param name="cutoff">Cutoff frequency wrt. the Nyquist frequency</param>
        /// <param name="window">Coefficients of the window</param>
        /// <returns>The windowed sinc filter coefficients</returns>
        public static double[] WindowedSinc(double cutoff, double[] window)
        {
            int length = window.Length;
            double[] f = new double[length];
            double[] h = LowpassSinc(cutoff, length);
            double sum = 0.0;
            for (int i = 0; i < length; i++)
            {
                f[i] = h[i] * window[i];
                sum += f[i];
            }
            sum = 1.0 / sum;
            for (int i = 0; i < length; i++)
            {
                f[i] *= sum;
            }
            return f;
        }

        /// <summary>
        /// Creates a windowed sinc lowpass filter using a Hamming window
        /// </summary>
        /// <param name="cutoff">Cutoff frequency wrt. the Nyquist frequency</param>
        /// <param name="length">Filter length - number of taps</param>
        /// <returns>The filter coefficients</returns>
        public static double[] HammingFilter(double cutoff, int length)
        {
            return WindowedSinc(cutoff, HammingWindow(length));
        }
        /// <summary>
        /// Creates a windowed sinc lowpass filter using a Kaiser window
        /// </summary>
        /// <param name="cutoff">Cutoff frequency wrt. the Nyquist frequency</param>
        /// <param name="length">Filter length - number of taps</param>
        /// <param name="alpha">Alpha parameter for the Kaiser window</param>
        /// <returns>The filter coefficients</returns>
        public static double[] KaiserFilter(double cutoff, int length, double alpha)
        {
            return WindowedSinc(cutoff, KaiserWindow(length, alpha));
        }

        /// <summary>
        /// Approximation to the modified bessel function of the first kind with alpha = 0
        /// </summary>
        public static double MBessel0(double x)
        {
            x = x / 2.0;
            double fact = 1.0;
            double ret = 0;
            double pow = 1.0;
            for (uint i = 0; i < 20; i++)
            {
                ret += pow / (fact * fact);
                pow = pow * x * x;
                fact = fact * (i + 1);
            }
            return ret;
        }

        /// <summary>
        /// Kaiser window. OK?
        /// </summary>
        /// <param name="n">Number of output coefficents</param>
        /// <param name="alpha">Alpha parameter. Look it up.</param>
        /// <returns>The coefficents of the window funtion</returns>
        public static double[] KaiserWindow(int n, double alpha)
        {
            double beta = alpha * Math.PI;
            double[] win = new double[n];
            double denom = 1.0 / MBessel0(beta);
            double delay = (n - 1.0) / 2;
            for (int i = 0; i < n; i++)
            {
                double x = i - delay;
                double arg = 1.0 - Math.Pow((2.0 * x) / (n - 1.0), 2.0);
                win[i] = MBessel0(beta * Math.Sqrt(arg));
                win[i] *= denom;
            }
            return win;
        }
    }
}
