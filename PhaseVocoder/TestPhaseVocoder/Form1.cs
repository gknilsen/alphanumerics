using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using PhaseVocoder;
 

namespace TestPhaseVocoder
{
    public partial class Form1 : Form
    {
        /* Author: Geir K.Nilsen, geir.kjetil.nilsen@gmail.com, 2014-2017 */

        public Form1()
        {
            InitializeComponent();

            /* Create a test signal: sinusoidal, 8192 samples long, frequency 440 Hz, sampling rate 8192, = duration 1 s */
            double[] x = new double[8192];
            for (int sample = 0; sample < 8192; sample++)
            {
                x[sample] = Math.Sin(2.0 * Math.PI * 440.0 * (double)sample / 8192.0);
            }

            PhaseVocoder.PhaseVocoder pvoc = null;

            /* TEST1. Test 2x tempo change on a 1s 440 Hz sine wave = 0.5s 440 Hz sine wave */
            pvoc = new PhaseVocoder.PhaseVocoder(2, 1, 1, 1, 128);
            double[] y = pvoc.apply(x);

            /* TEST2. Test 0.5x tempo change on a 1s 440 Hz sine wave = 2s 440 Hz sine wave */
            pvoc = new PhaseVocoder.PhaseVocoder(1, 2, 1, 1, 128);
            y = pvoc.apply(x);

            /* TEST3. Test 0.5x pitch shift on a 1s 440 Hz sine wave = 1s 220 Hz sine wave */
            pvoc = new PhaseVocoder.PhaseVocoder(1, 1, 1, 2, 128);
            y = pvoc.apply(x);

            /* TEST4. Test 2x pitch shift on a 1s 440 Hz sine wave = 1s 880 Hz sine wave */
            pvoc = new PhaseVocoder.PhaseVocoder(1, 1, 2, 1, 128);
            y = pvoc.apply(x);

            /* TEST5. Test 2x tempo change on a 10s 440Hz sine wave (block by block) = 5s 440 Hz sine wave in continuous blocks */
            pvoc = new PhaseVocoder.PhaseVocoder(2, 1, 1, 1, 128);

            List<double> inputList = new List<double>();
            List<double> outputList = new List<double>();

            double[] output;
            for (int block = 0; block < 10; block++)
            {
                double[] input = new double[1024];
                for (int sample = 0; sample < 1024; sample++)
                {
                    if (block * 1024 + sample >= x.Length)
                        break;

                    input[sample] = x[block * 1024 + sample];
                }

                inputList.AddRange(input.ToList());

                output = pvoc.applyContinuous(input);

                outputList.AddRange(output.ToList());
            }

            output = pvoc.applyContinuous(null);
            outputList.AddRange(output.ToList());

            /* Now compare inputList and outputList */

            /* TEST6. Test 0.5x tempo change on a 10s 440Hz sine wave (block by block) = 20s 440 Hz sine wave in continuous blocks */
            pvoc = new PhaseVocoder.PhaseVocoder(1, 2, 1, 1, 128);

            for (int block = 0; block < 10; block++)
            {
                double[] input = new double[1024];
                for (int sample = 0; sample < 1024; sample++)
                {
                    if (block * 1024 + sample >= x.Length)
                        break;

                    input[sample] = x[block * 1024 + sample];
                }

                inputList.AddRange(input.ToList());

                output = pvoc.applyContinuous(input);

                outputList.AddRange(output.ToList());
            }

            output = pvoc.applyContinuous(null);
            outputList.AddRange(output.ToList());

            /* Now compare inputList and outputList */

            /* TEST7. Test 2x pitch shift on a 1s 440 Hz sine wave = 1s 880 Hz sine wave in continuous blocks */
            pvoc = new PhaseVocoder.PhaseVocoder(1, 1, 2, 1, 128);

            for (int block = 0; block < 10; block++)
            {
                double[] input = new double[1024];
                for (int sample = 0; sample < 1024; sample++)
                {
                    if (block * 1024 + sample >= x.Length)
                        break;

                    input[sample] = x[block * 1024 + sample];
                }

                inputList.AddRange(input.ToList());

                output = pvoc.applyContinuous(input);

                outputList.AddRange(output.ToList());
            }

            output = pvoc.applyContinuous(null);
            outputList.AddRange(output.ToList());

            /* Now compare inputList and outputList */


            /* TEST8. Test 0.5x pitch shift on a 1s 440 Hz sine wave = 1s 220 Hz sine wave */
            pvoc = new PhaseVocoder.PhaseVocoder(1, 1, 1, 2, 128);

            for (int block = 0; block < 10; block++)
            {
                double[] input = new double[1024];
                for (int sample = 0; sample < 1024; sample++)
                {
                    if (block * 1024 + sample >= x.Length)
                        break;

                    input[sample] = x[block * 1024 + sample];
                }

                inputList.AddRange(input.ToList());

                output = pvoc.applyContinuous(input);

                outputList.AddRange(output.ToList());
            }

            output = pvoc.applyContinuous(null);
            outputList.AddRange(output.ToList());

            /* Now compare inputList and outputList */


        }
    }
}
