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
