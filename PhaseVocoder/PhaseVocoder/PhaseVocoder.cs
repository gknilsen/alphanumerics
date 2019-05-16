using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Numerics;

namespace PhaseVocoder
{
    /* Phase Vocoder class capable of pitch shifting and tempo changing of signals 
     *
     * Author: Geir K.Nilsen, geir.kjetil.nilsen@gmail.com, 2014-2017
     */
    public class PhaseVocoder
    {
        private int Qt, Pt, Qp, Pp;
        private bool robotoFX = false;
        private double[] prev_input_signal = null;            
        private double[] prev_output_signal = null;
        private double[] fi = null;
        private double[] dfi = null;
        private STFT stft = null;
        private TFR tfr = null;
        Resample.Resampler resampler;
        private enum FRAMETYPE {START, MID, END};


        /*
         * Pt/Qt = tempo change factor. A tempo change factor < 1 will lead to slower tempo. A tempo change factor > 1 will lead to faster tempo
         * Pp/Qp = pitch shift factor. A pitch shift factor < 1 will lead to reduced/lower pitch. A pitch shift factor > 1 will lead to increased/raised pitch.
         * 
         */

        public PhaseVocoder(int Pt, int Qt, int Pp, int Qp)
        {
            this.Qt = Qt;
            this.Pt = Pt;
            this.Qp = Qp;
            this.Pp = Pp;
            this.robotoFX = false;
            tfr = new TFR();
            stft = new STFT(128);
            resampler = new Resample.Resampler((uint)Qp, (uint)Pp);
        }

        public PhaseVocoder(int Pt, int Qt, int Pp, int Qp, int winLength)
        {
            this.Qt = Qt;
            this.Pt = Pt;
            this.Qp = Qp;
            this.Pp = Pp;
            this.robotoFX = false;
            tfr = new TFR();
            stft = new STFT(winLength);
            resampler = new Resample.Resampler((uint)Qp, (uint)Pp);
        }

        public PhaseVocoder(int Pt, int Qt, int Pp, int Qp, bool robotoFX, int winLength)
        {
            this.Qt = Qt;
            this.Pt = Pt;
            this.Qp = Qp;
            this.Pp = Pp;
            this.robotoFX = robotoFX;
            tfr = new TFR();
            stft = new STFT(winLength);
            resampler = new Resample.Resampler((uint)Qp, (uint)Pp);
        }

        public void reset()
        {
        }

        public double[] apply(double[] input_signal)
        {
            tfr = stft.forward(input_signal);
            tfr.applyTempoChange((double)Pt/(double)Qt * (double)Qp/(double)Pp, false, this.robotoFX, ref fi, ref dfi); /* Tempo change step */
            double[] output_signal = stft.inverse(tfr);

            /* Pitch shift step. Apply function resample(Qp, Pp) on output_signal to change its sample rate by a factor Qp/Pp (multiplied by) */

            output_signal = resampler.Resample(output_signal);

            return output_signal;
        }
        
        /* For real-time usage, ensure that the output blocks are continuous */
        public double[] applyContinuous(double[] input_signal)
        {
            return applyContinuous(input_signal, this.Pt, this.Qt, this.Pp, this.Qp, this.robotoFX);
        }

        /* For real-time usage, ensure that the output blocks are continuous */
        public double[] applyContinuous(double[] input_signal, int Pt, int Qt, int Pp, int Qp)
        { 
            return applyContinuous(input_signal, Pt, Qt, Pp, Qp, this.robotoFX);
        }

        /* For real-time usage, ensure that the output blocks are continuous */
        public double[] applyContinuous(double[] input_signal, int Pt, int Qt, int Pp, int Qp, bool robotoFX)
        {            
            List<double[]> frame = prepFrame(input_signal, prev_input_signal);

            for (int frameNum = 0; frameNum < frame.Count; frameNum++)
            {
                tfr = stft.forward(frame[frameNum]);
                tfr.applyTempoChange((double)Pt / (double)Qt * (double)Qp / (double)Pp, true, robotoFX, ref fi, ref dfi); /* Tempo change step */
                if(input_signal == null)
                    frame[frameNum] = stft.inverse(tfr, ref prev_output_signal, true);
                else
                    frame[frameNum] = stft.inverse(tfr, ref prev_output_signal, false);
            }

            /* Combine blocks to produce single output block */
            double[] output_signal;
            if (frame.Count == 1)
            {
                if (input_signal == null)
                    output_signal = bondFrames(frame, FRAMETYPE.START);
                else
                    output_signal = bondFrames(frame, FRAMETYPE.END);
            }
            else
                output_signal = bondFrames(frame, FRAMETYPE.MID);

            if (input_signal != null)
            {
                prev_input_signal = new double[input_signal.Length];
                for (int sample = 0; sample < input_signal.Length; sample++)
                    prev_input_signal[sample] = input_signal[sample];
            }

            /* Pitch shift step. Apply OUTPUT-BLOCK-CONTINUOUS function resampleContinuous(Qp, Pp) on output_signal to change its sample rate by a factor Qp/Pp (multiplied by) */
            output_signal = resampler.ResampleContinuous(output_signal);

            return output_signal;
        }

        public double[] getSpectrum()
        {
            Complex[] tmp = tfr.getTFRAsDoubleArray()[2];
            double[] spectrum = new double[tmp.Length];
            for (int s = 0; s < tmp.Length; s++)
                spectrum[s] = tmp[s].Magnitude;

            return spectrum;
        }

        public double[] getPhaseTimeSlice(int freq)
        {
            Complex[][] tmp = tfr.getTFRAsDoubleArray();

            double[] slice = new double[tmp.Length];
            for (int x = 0; x < tmp.Length; x++)
                for (int y = 0; y < tmp[x].Length; y++)
                    slice[x] = tmp[x][freq].Phase;

            return slice;
        }

        private double[] bondFrames(List<double[]> block, FRAMETYPE frameType)
        {
            List<double> output;
            switch(frameType)
            {
                case FRAMETYPE.START:
                    output = block[0].ToList();                  
                    return output.ToArray();
                    break;

                case FRAMETYPE.END:
                    output = block[0].ToList();
                    return output.ToArray();
                    break;

                case FRAMETYPE.MID:
                    List<double> output1 = block[0].ToList();
                    List<double> output2 = block[1].ToList();
                    output1.AddRange(output2);
                    return output1.ToArray();
                    break;
            }

            return null;
        }

        private List<double[]> prepFrame(double[] signal, double[] prevSignal)
        {
            /* no previous input signal */
            if (prevSignal == null)
            {
                int mid = signal.Length / 2;
                int ol = stft.hopSize * 3;
                int start = 0;
                int stop = mid + ol;

                double[] block = new double[stop - start];
                for (int sample = start; sample < stop; sample++)
                {
                    block[sample - start] = signal[sample];
                }

                List<double[]> list = new List<double[]>();
                list.Add(block);
                return list;

            }
            else
            {
                /* Previous block, but the current frame is null => last call */

                if (signal == null)
                {
                    int mid = prevSignal.Length / 2;
                    int ol = stft.hopSize * 3;
                    int start = mid;
                    int stop = prevSignal.Length;

                    double[] block = new double[stop - start];
                    for (int sample = start; sample < stop; sample++)
                    {
                        block[sample - start] = prevSignal[sample];
                    }

                    List<double[]> list = new List<double[]>();
                    list.Add(block);
                    return list;
                }
                else
                {
                    /* Previous block and current signal => middle of successive calls */
                    int mid = signal.Length / 2;
                    int ol = stft.hopSize * 3;
                    int start = mid;
                    int stop = signal.Length;

                    double[] subblock1 = new double[stop - start + ol];

                    for (int sample = start; sample < stop; sample++)
                    {
                        subblock1[sample - start] = prevSignal[sample];
                    }

                    for (int sample = 0; sample < ol; sample++)
                    {
                        subblock1[subblock1.Length - ol + sample] = signal[sample];
                    }

                    List<double[]> list = new List<double[]>();
                    list.Add(subblock1);

                    start = 0;
                    stop = mid + ol;

                    double[] subblock2 = new double[stop - start];

                    for (int sample = start; sample < stop; sample++)
                    {
                        subblock2[sample - start] = signal[sample];
                    }
                    list.Add(subblock2);

                    return list;
                }
            }
        }

    }
}

