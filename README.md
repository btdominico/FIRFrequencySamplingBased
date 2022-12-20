# FIR Frequency Sampling-Based 

## Example on how to design a filter 
```c#
double fs = 44100;
double[] fHz = {0,250,500,1e3,2e3, 3e3, 4e3, 6e3, 8e3, fs/2};
double[] vaudiogram = {6, 0, 0, -20, -30, -40, -50, -60};

// Filter Order M
double filterM = 1000;

// Calculate the Vector of N frequencies
double[] vFreqs = new double[fHz.Length];
for (int i = 0; i < fHz.Length; i++)
    vFreqs[i] = fHz[i] / (fs / 2.0);

// Calculate the Vector of Magnitudes
double[] vA = new double[vaudiogram.Length];
for (int i = 0; i < vaudiogram.Length; i++)
    vA[i] = Math.Pow(10, (double)vaudiogram[i] / 20.0);

double[] vMags = new double[vA.Length + 2];
vMags[0] = vA[0];
vMags[vMags.Length - 1] = vA[vA.Length - 1];

for (int i = 1; i < vMags.Length - 1; i++)
    vMags[i] = vA[i - 1];

// Calculate the Filter Coefficients
double[] h = new double[(int)filterM + 1];
Filter(filterM, vFreqs, vMags, out h);

// Print the Filter Coefficients
Console.WriteLine("Filter Coefficients: ");
for (int i = 0; i <= filterM; i++)
    Console.WriteLine(h[i]);
```

## Example on how to create a test signal composed of eight one-second signals from 250Hz to 8kHz to validate the accuracy of the filter designed 
```c#
// Define variables
int fs = 44100;
double[] x = new double[fs*8]; // input vector
double[] fHz = {250, 500, 1e3, 2e3, 3e3, 4e3, 6e3, 8e3};

// Create the input signal x of 1s of each frequency in fHz
double t = 0;   // current time 
for (int iFreq = 0; iFreq < 8; iFreq++, t += 1.0/fs)
    for (int i = 0; i < fs; i++)
        x[i+iFreq*fs] = Math.Cos(2*Math.PI*fHz[iFreq]*t);

// Filter x using the filter h created in Fig. 5
double[] y = new double[x.Length]; // output vector
for (int i = 0; i < x.Length; i++)
{
    for (int j = 0; j < h.Length; j++)
    {
        if (i - j < 0) continue;
        y[i] += h[j] * x[i-j];
    }
}
Console.WriteLine("Filtered output Data: ");
for (int i = 0; i <= filterM; i++)
    Console.WriteLine(y[i]);

```
