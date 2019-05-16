using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace PhaseVocoder
{
    /* Time-frequency representation class 
     *
     * Author: Geir K.Nilsen, geir.kjetil.nilsen@gmail.com, 2014-2017
     */

    public class TFR
    {
        private Complex[][] tfr;

        public TFR()
        {
            this.tfr = null;
        }

        public TFR(Complex[][] tfr)
        {
            this.tfr = tfr;
        }

        public void applyTempoChange(double tempoChangeFactor, bool continuous, bool robotoFX, ref double[] fi, ref double[] dfi)
        {
            // Set initial phase to the phase of 1) the first frame plus 2) the phase of the last frame of the previous TFR.
            // 1) Phase of first frame: starting point for phase cumulation (standard phase vocoder)
            // 2) Phase of the last from of the previous TFR: This trick ensures that we will have correct phase transitions when running in continuous mode.

            double[] initPhase = this.getPhase()[0];
            double[] phase = new double[initPhase.Length];

            double[] p = this.getPhase()[0];

            if (continuous)
            {
                if (fi == null)
                {
                    for (int sample = 0; sample < initPhase.Length; sample++)
                        phase[sample] = initPhase[sample];
                }
                else
                {
                    for (int sample = 0; sample < initPhase.Length; sample++)
                        phase[sample] = fi[sample] + (p[sample] - dfi[sample]);
                }
            }
            else
            {
                for (int sample = 0; sample < initPhase.Length; sample++)
                    phase[sample] = initPhase[sample];
            }

            int nf = this.getTFRAsDoubleArray().Length - 1;
                        
            Complex[][] pframe = new Complex[(int) (((double)nf / tempoChangeFactor) + 1)][];

            /* Avoid boundary problem when interpolating by appending extra zero frame */

            Complex[][] buf = new Complex[tfr.Length + 1][];
            for (int frame = 0; frame < tfr.Length; frame++)
                buf[frame] = tfr[frame];
            buf[tfr.Length] = new Complex[tfr[0].Length];
            for (int sample = 0; sample < tfr[0].Length; sample++)
                buf[tfr.Length][sample] = new Complex();

            Complex[] frame1 = null;
            Complex[] frame2 = null;

            double[] mag = new double[buf[0].Length];
            double[] phase_adv = new double[buf[0].Length];

            double frameInterPolant = 0.0;
            for(int frameInd = 0; frameInd < pframe.Length; frameInd++)
            {              
                frame1 = buf[(int) (Math.Floor(frameInterPolant))];
                frame2 = buf[(int) (int)Math.Floor(frameInterPolant) + 1];

                double scale = frameInterPolant - Math.Floor(frameInterPolant);
                
                pframe[frameInd] = new Complex[frame1.Length];

                for(int sample = 0; sample < frame1.Length; sample++)
                {
                    /* Magnitude interpolation */
                    mag[sample] = (1.0 - scale) * frame1[sample].Magnitude + scale * frame2[sample].Magnitude;

                    pframe[frameInd][sample] = mag[sample] * Complex.Exp(new Complex(0.0, phase[sample]));

                    phase_adv[sample] = frame2[sample].Phase - frame1[sample].Phase;

                    if(!robotoFX)
                        phase[sample] = phase[sample] + phase_adv[sample];
                }

                frameInterPolant += tempoChangeFactor;
            }

            // Save last phase frame for next call.
            fi = new double[phase.Length];
            dfi = new double[phase.Length];
            for (int sample = 0; sample < pframe[pframe.Length - 1].Length; sample++)
            {
                fi[sample] = pframe[pframe.Length - 1][sample].Phase;
                dfi[sample] = tfr[tfr.Length - 1][sample].Phase;
            }                     

            tfr = pframe;
        }

        public Complex[][] getTFRAsDoubleArray()
        {
            return tfr;
        }

        
        public double[][] getMagnitude()
        {
            double[][] mag = new double[tfr.Length][];
            for (int x = 0; x < tfr.Length; x++)
            {
                mag[x] = new double[tfr[x].Length];

                for (int y = 0; y < tfr[x].Length; y++)
                    mag[x][y] = tfr[x][y].Magnitude;
            }
            return mag;
        }
                
        public double[][] getPhase()
        {
            double[][] phase = new double[tfr.Length][];
            for (int x = 0; x < tfr.Length; x++)
            {
                phase[x] = new double[tfr[x].Length];

                for (int y = 0; y < tfr[x].Length; y++)
                    phase[x][y] = tfr[x][y].Phase;
            }
            return phase;
        }
    }

}
