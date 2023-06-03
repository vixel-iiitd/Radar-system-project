# Radar-system-project
❖ Objective
Range-Doppler ambiguity graphs of a target using Golay transmitted signal
- Perform matched filtering for range estimation in the time domain and frequency
domain
- Perform 2D range-Doppler processing for a moving target
- Repeat exercise for multiple targets (three)
- Perform cell averaged CFAR-based detection of target parameters (amplitude,
range, and Doppler)
- Repeat the above exercises for Doppler resilient Golay sequences

❖ Methodology

➢ Transmitted and Received Signals
For the Golay transmitted signal, 512 samples of Golay sequence in 1 PRI have
been taken. The signal is then zero-padded to have a duty cycle of 50%, resulting
in 1024 samples in 1 PRI. For the received signal, the transmitted signal is shifted
according to the number of samples corresponding to the time delay.

➢ Range Estimation
We know that, td (time delay) = 2rtgt / c (where rtgt = range of target)
Given, dt = 1/BW
Sample delay, nd = round(td/dt)
- Range Estimation in time domain
Time domain Matched filtering was performed using the xcorr function.
- Range Estimation in frequency domain
We know that convolution in the time domain is equivalent to multiplication in the
frequency domain. Firstly, the received and transmitted signals are converted from
time to frequency domain, using the Fast Fourier Transform fft function. Then the
received signal is multiplied to the conjugate of the transmitted signal. Finally,
Inverse Fast Fourier Transform is performed using the ifft function.

➢ 2D Range Doppler Processing
For 2D Range Doppler Processing, the considered signal is a 2D array having M
rows and N columns, in our case M = 1024 and N = 2048. The number of columns
N represents the number of PRI in one CPI, also known as the slow time samples.
The number of rows M represents the fast time samples in 1 PRI.
Radar System Project 1
Range delay and doppler delay are computed to construct the received signal.
Range delay is incorporated by shifting the signal according to the number of
samples corresponding to the time delay. Doppler delay is incorporated by
multiplying each column with e
j2π.fd.xTpri where fd is the doppler frequency
corresponding to the target's velocity, and x representing the column number
starting from 0 to N-1.

➢ 2D Range Doppler For single target
2D Range Doppler processing is performed by doing Matched filtering across
rows in a column and doing Doppler Processing across columns in a row. Doppler
Processing was performed using Fast Fourier Transform(FFT). We know that for
Matched filtering, if the input arguments consist of n samples, we get 2n-1 samples
in the output. We also know that for FFT, if the input argument consists of n
samples, we get n samples in the output. Thus, the dimensions after 2D Range
Doppler processing were 2047 x 2048.

➢ 2D Range Doppler For Multiple targets
For 2D Range Doppler for multiple targets, we first fixed the target ranges, and
velocities, like for 3 targets we made two array of target ranges and velocities and
then we looped till 1 to 3, and in each iteration we are finding the 2d range
doppler map for each target, we are superimposing(adding) the each 2d range
doppler map to a matrix initialized to the same dimension and set to zero initially.
Then Finally after the loop we are using imagesc to plot the superimposed 2d
range doppler map.

➢ CFAR
For CFAR, we first formed the CFARDetector2D object using a phased toolbox
command, with parameters such as Training Band Size = [5,4], Guard Band Size =
[4,4], false alarm rate probability as 1e-5. Now to find the object spots in the range
doppler image, we have to make cut windows using the bands we have defined.
For forming windows, we have looped over the start and ending of rows and
columns to append all the possible windows into 2d array cutidx. Now we use the
object we had created to find the detections map using the cuts we have created.
After that, we plot the detections map using helperDetectionMaps by MATLAB to
plot the detections.
