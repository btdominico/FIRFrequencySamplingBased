# FIR Frequency Sampling-Based 

## Example on how to use the code
```c#
double fs = 44100;
double[] fHz = { 0, 250, 500, 1000, 2000, 3000, 4000, 6000, 8000, fs/2.0 };
double[] vaudiogram = { 10, 0, 0, -20, -30, -40, -50, -60 };
double filterN = 1000;

double[] vA = new double[vaudiogram.Length];
for (int i = 0; i < vaudiogram.Length; ++i)
{
    vA[i] = Math.Pow(10, (double)vaudiogram[i] / 20.0);
}

double[] vANew = new double[vA.Length + 2];
vANew[0] = vA[0];
vANew[vANew.Length - 1] = vA[vA.Length - 1];

for (int i = 1; i < vANew.Length - 1; i++)
{
    vANew[i] = vA[i - 1];
}

double[] fnormal = new double[fHz.Length];
for (int j = 0; j < fHz.Length; ++j)
{
    fnormal[j] = fHz[j] / (fs / 2.0);
}

double[] fir = new double[(int)filterN + 1];
Filter(filterN, fnormal, vANew, out fir);
```

## Example on how to build the diagram using the code
```c#
double[] resultF = new double[8]; 
testFilter(fir, out resultF);
```
