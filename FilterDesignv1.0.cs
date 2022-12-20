using System;
using System.Numerics;
using System.Diagnostics;

// filter the data
static double[] Filter(double filterM, double[] vFreqs, double[] vMags, out double[] h)
{
    double gridNum;
    double[] x;
    double[] w;
    double[] window;

    // Work with filter length instead of filter order
    filterM = filterM + 1;

    if (filterM < 1024)
        gridNum = 512;
    else
        gridNum = Math.Pow(2, Math.Ceiling(Math.Log(filterM) / Math.Log(2)));

    double filterM_round = Math.Round(filterM);

    double filterM_half = (filterM_round + 1) / 2.0;

    x = new Double[(int)(filterM_half)];
    w = new Double[(int)(filterM_half)];
    window = new Double[(int)(filterM_half * 2 - 1)];

    for (int i = 0; i < filterM_half; i++)
    {
        x[i] = i / (filterM_round - 1);
        w[i] = 0.54 - 0.46 * Math.Cos(2 * Math.PI * x[i]);
        window[i] = w[i];
        window[window.Length - i - 1] = w[i];
    }

    if (filterM != window.Length)
    {
        Console.WriteLine("nn and wind length are not the same");
    }

    int vFreqsLength = vFreqs.Length;
    int vMagsLength = vMags.Length;

    if (vFreqsLength != vMagsLength)
    {
        Console.WriteLine("vFreqs and vMags must have the same length");
    }

    double eps = 2.2204e-16;

    if (Math.Abs(vFreqs[0]) > eps || Math.Abs(vFreqs[vFreqsLength - 1] - 1) > eps)
    {
        Console.WriteLine("Invalid Range");
    }

    vFreqs[0] = 0;
    vFreqs[vFreqsLength - 1] = 1;

    if (filterM > 2.0 * gridNum)
    {
        Console.WriteLine("Invalid Npt");
    }

    gridNum = gridNum + 1;
    double[] H = new Double[(int)gridNum];

    int nb = 1;

    H[0] = vMags[0];
    int temp = 0;

    for (int i = 0; i < (vFreqsLength - 1); i++)
    {
        int ne = (int)Math.Floor(vFreqs[i + 1] * gridNum);
        if (nb < 0 || ne > gridNum)
        {
            Console.WriteLine("Signal Error");
        }

        double[] j = new Double[(int)ne];

        for (int m = temp; m < ne; m++)
        {
            j[m] = m + 1;
        }

        double[] inc = new Double[ne];

        if (nb == ne)
        {
            for (int k = temp; k < ne; k++)
            {
                inc[k] = 0;
            }
        }

        else
        {
            for (int k = temp; k < ne; k++)
            {
                inc[k] = (j[k] - nb) / (ne - nb);
            }
        }

        for (int l = temp; l < ne; l++)
        {
            H[l] = inc[l] * vMags[i + 1] + (1 - inc[l]) * vMags[i];
        }
        nb = ne + 1;
        temp = ne;
    }

    double dt = 0.5 * (filterM - 1);

    Complex[] Ynew1 = new Complex[(int)gridNum];
    for (int i = 0; i < gridNum; i++)
    {
        Ynew1[i] = new Complex(H[i], 0);
    }

    Complex[] Ynew = new Complex[(int)gridNum];
    Complex[] rad = new Complex[(int)gridNum];
    for (int i = 0; i < gridNum; i++)
    {
        rad[i] = new Complex(0, -dt * Math.PI * i / (gridNum - 1));
        Ynew[i] = Ynew1[i] * Complex.Exp(rad[i]);
    }


    Complex[] HnewF = new Complex[2 * (int)gridNum - 2];
    for (int i = 0; i < gridNum; i++)
    {
        HnewF[i] = Complex.Conjugate(Ynew[i]);
        if (i <= gridNum - 2)
            HnewF[2 * (int)gridNum - 3 - i] = Ynew[i + 1];
    }

    double[] Xr = new double[HnewF.Length];
    double[] Xi = new double[HnewF.Length];
    double[] Yr = new double[HnewF.Length];

    for (int i = 0; i < HnewF.Length; i++)
    {
        Xr[i] = HnewF[i].Real;
        Xi[i] = Complex.Conjugate(HnewF[i]).Imaginary;
    }
    JimNumerics.JimFFT.Inverse(Xr, Xi, out Yr);

    h = new double[(int)filterM];
    for (int i = 0; i < filterM; i++)
    {
        h[i] = Yr[i] * window[i];
    }
    return h;
}

static double[] OneDimensionFilter(double[] h, double[] signal, out double[] y)
{
    y = new double[signal.Length];
    for (int i = 0; i < signal.Length; i++)
    {
        for (int j = 0; j < h.Length; j++)
        {
            if (i - j < 0) continue;
            y[i] += h[j] * signal[i - j];
        }
    }
    return y;
}

static double[] testFilter(double[] b, out double[] result)
{
    double[] signal = { 250, 500, 1000, 2000, 3000, 4000, 6000, 8000 };
    int fs = 22050;
    int lengtht = fs * 2 + 1;
    double[] t = new double[lengtht];
    for (int i = 0; i < lengtht; i++)
    {
        if (i == 0)
            t[i] = 0;
        else
            t[i] = t[i - 1] + (1.0/fs);
    }
    double[] data250 = new double[lengtht];
    double[] data500 = new double[lengtht];
    double[] data1k = new double[lengtht];
    double[] data2k = new double[lengtht];
    double[] data3k = new double[lengtht];
    double[] data4k = new double[lengtht];
    double[] data6k = new double[lengtht];
    double[] data8k = new double[lengtht];
    result = new double[8];

    for (int j = 0; j < lengtht; j++)
    {
        data250[j] = Math.Cos(2 * Math.PI * signal[0] * t[j]);
        data500[j] = Math.Cos(2 * Math.PI * signal[1] * t[j]);
        data1k[j] = Math.Cos(2 * Math.PI * signal[2] * t[j]);
        data2k[j] = Math.Cos(2 * Math.PI * signal[3] * t[j]);
        data3k[j] = Math.Cos(2 * Math.PI * signal[4] * t[j]);
        data4k[j] = Math.Cos(2 * Math.PI * signal[5] * t[j]);
        data6k[j] = Math.Cos(2 * Math.PI * signal[6] * t[j]);
        data8k[j] = Math.Cos(2 * Math.PI * signal[7] * t[j]);
    }

    double[] y250 = new double[lengtht];
    double[] y500 = new double[lengtht];
    double[] y1k = new double[lengtht];
    double[] y2k = new double[lengtht];
    double[] y3k = new double[lengtht];
    double[] y4k = new double[lengtht];
    double[] y6k = new double[lengtht];
    double[] y8k = new double[lengtht];

    OneDimensionFilter(b, data250, out y250);
    OneDimensionFilter(b, data500, out y500);
    OneDimensionFilter(b, data1k, out y1k);
    OneDimensionFilter(b, data2k, out y2k);
    OneDimensionFilter(b, data3k, out y3k);
    OneDimensionFilter(b, data4k, out y4k);
    OneDimensionFilter(b, data6k, out y6k);
    OneDimensionFilter(b, data8k, out y8k);

    double[] y250N = new double[fs];
    double[] y500N = new double[fs];
    double[] y1kN = new double[fs];
    double[] y2kN = new double[fs];
    double[] y3kN = new double[fs];
    double[] y4kN = new double[fs];
    double[] y6kN = new double[fs];
    double[] y8kN = new double[fs];

    for (int i = 0; i < fs; i++)
    {
        y250N[i] = y250[fs + i];
        y500N[i] = y500[fs + i];
        y1kN[i] = y1k[fs + i];
        y2kN[i] = y2k[fs + i];
        y3kN[i] = y3k[fs + i];
        y4kN[i] = y4k[fs + i];
        y6kN[i] = y6k[fs + i];
        y8kN[i] = y8k[fs + i];
    }

    double dB250 = 0;
    double dB500 = 0;
    double dB1k = 0;
    double dB2k = 0;
    double dB3k = 0;
    double dB4k = 0;
    double dB6k = 0;
    double dB8k = 0;

    double dB250N = 0;
    double dB500N = 0;
    double dB1kN = 0;
    double dB2kN = 0;
    double dB3kN = 0;
    double dB4kN = 0;
    double dB6kN = 0;
    double dB8kN = 0;

    for (int i = 0; i < fs; i++)
    {
       dB250 = Math.Max(Math.Abs(y250N[i]), dB250);
       dB500 = Math.Max(Math.Abs(y500N[i]), dB500);
       dB1k = Math.Max(Math.Abs(y1kN[i]), dB1k);
       dB2k = Math.Max(Math.Abs(y2kN[i]), dB2k);
       dB3k = Math.Max(Math.Abs(y3kN[i]), dB3k);
       dB4k = Math.Max(Math.Abs(y4kN[i]), dB4k);
       dB6k = Math.Max(Math.Abs(y6kN[i]), dB6k);
       dB6k = Math.Max(Math.Abs(y8kN[i]), dB8k);
    }

    dB250N = 20 * Math.Log10(dB250);
    dB500N = 20 * Math.Log10(dB500);
    dB1kN = 20 * Math.Log10(dB1k);
    dB2kN = 20 * Math.Log10(dB2k);
    dB3kN = 20 * Math.Log10(dB3k);
    dB4kN = 20 * Math.Log10(dB4k);
    dB6kN = 20 * Math.Log10(dB6k);
    dB6kN = 20 * Math.Log10(dB8k);
   

    result[0] = dB250N;
    result[1] = dB500N;
    result[2] = dB1kN;
    result[3] = dB2kN;
    result[4] = dB3kN;
    result[5] = dB4kN;
    result[6] = dB6kN;
    result[7] = dB8kN;

    return result;
}


// Create input signal with 1s at each frequency
double fs = 44100;
double t = 0;
double[] y1 = new double[(int)fs + 1];
double[] y2 = new double[(int)fs + 1];
double[] y3 = new double[(int)fs + 1];
double[] y4 = new double[(int)fs + 1];
double[] y5 = new double[(int)fs + 1];
double[] y6 = new double[(int)fs + 1];
double[] y7 = new double[(int)fs + 1]; 
double[] y8 = new double[(int)fs + 1];

for (int i = 0; i <= (int)fs; i++)
{
    y1[i] = Math.Cos(2 * Math.PI * 250 * t);
    y2[i] = Math.Cos(2 * Math.PI * 500 * t);
    y3[i] = Math.Cos(2 * Math.PI * 1000 * t);
    y4[i] = Math.Cos(2 * Math.PI * 2000 * t);
    y5[i] = Math.Cos(2 * Math.PI * 3000 * t);
    y6[i] = Math.Cos(2 * Math.PI * 4000 * t);
    y7[i] = Math.Cos(2 * Math.PI * 6000 * t);
    y8[i] = Math.Cos(2 * Math.PI * 8000 * t);
    t = t + (1.0 / fs);
}

double[] inputData = new double[(int)fs * 8 + 8];
y1.CopyTo(inputData, 0);
y2.CopyTo(inputData, (int)fs * 1 + 1);
y3.CopyTo(inputData, (int)fs * 2 + 2);
y4.CopyTo(inputData, (int)fs * 3 + 3);
y5.CopyTo(inputData, (int)fs * 4 + 4);
y6.CopyTo(inputData, (int)fs * 5 + 5);
y7.CopyTo(inputData, (int)fs * 6 + 6);
y8.CopyTo(inputData, (int)fs * 7 + 7);

// Print the input data 
Console.WriteLine("Input Data: ");
Console.WriteLine(inputData);

// Design filter
double[] fHz = { 0, 250, 500, 1000, 2000, 3000,
    4000, 6000, 8000, fs/2.0 };
double[] vaudiogram = { 6, 0, 0, -20, -30, -40, -50, -60 };

// Filter Order M
double filterM = 1000;

// Calculate the Vector of N frequencies
double[] vFreqs = new double[fHz.Length];
for (int i = 0; i < fHz.Length; i++)
{
    vFreqs[i] = fHz[i] / (fs / 2.0);
}

// Calculate the Vector of Magnitudes
double[] vA = new double[vaudiogram.Length];
for (int i = 0; i < vaudiogram.Length; i++)
{
    vA[i] = Math.Pow(10, (double)vaudiogram[i] / 20.0);
}

double[] vMags = new double[vA.Length + 2];
vMags[0] = vA[0];
vMags[vMags.Length - 1] = vA[vA.Length - 1];

for (int i = 1; i < vMags.Length - 1; i++)
{
    vMags[i] = vA[i - 1];
}

// Calculate the Filter Coefficients
double[] h = new double[(int)filterM + 1];
Filter(filterM, vFreqs, vMags, out h);

double[] y = new double[h.Length];
OneDimensionFilter(h,inputData,out y);

Console.WriteLine("A:");
