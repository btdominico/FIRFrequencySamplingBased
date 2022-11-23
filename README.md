# FIR Frequency Sampling-Based 

## Example on how to use the code
```c#
double fs = 44100;
double[] fHz = { 0, 250, 500, 1000, 2000, 3000, 
    4000, 6000, 8000, fs/2.0 };
double[] vaudiogram = { 10, 0, 0, -20, -30, -40, -50, -60 };

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
double[] fir = new double[(int)filterM + 1];
Filter(filterM, vFreqs, vMags, out fir);

// Print out the Filter Coefficients
Console.WriteLine("Filter Coefficients: ");
for (int i = 0; i <= filterM; i++)
{
    Console.WriteLine(fir[i]);
}
```

## Example on how to test the code
```c#
double[] resultF = new double[8]; 
testFilter(fir, out resultF);
```
