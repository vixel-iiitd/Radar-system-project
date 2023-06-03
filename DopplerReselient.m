


%Instructions:
%Golay_sequence.mat should be in same directory
clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global constant
c = 3e8;

%Radar parameters
% Duration of pulse in one tpri (seconds)
pulse_width = 2e-6; 
%Number of packets in one pulse (no unit)
subpulse_n = 512; 
%Duration of one sub pulse (seconds)
subpulse_width = pulse_width/subpulse_n; 
%Duty cycle, T_pri/pulse_width (no unit)
duty_cycle = 0.5; 
%Carrier frequency (Hz)
fc = 60e9; 
% Bandwidth (Hz)
BW = 1.76e9;
% Pulse repetition interval (in seconds)
tpri = pulse_width/duty_cycle; 
%Sampling the signal with sampling time as dt (in seconds)
dt = 1/BW;
% Max unambiguous range (m)
r_max = c*tpri/2;

% Target parameters
% Distance of the target from radar (m)
target_range = 60; 
% Time delay of received signal for the target (seconds)
td = (2*target_range)/c; 

%Number of samples per PRI (no unit)
nsamples = round(tpri/dt); 
%Number of samples in one sub pulse (no unit)
nsamples_perpulse = round(pulse_width/dt); 
% No. of samples corresponding to time delay for the target(no unit)
nd = round(td/dt); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Perform matched filtering for range estimation in the %%%%%%%%%%%%
%%%%%%%%%%%%%%%% time domain and frequency domain %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("Doppler_resilent_Golay_sequence.mat");
subpulse_code = xmat(:,1);
%Plotting the pseudo random number
figure; %Adding separate figure for plot
plot(subpulse_code); %Plotting the pseudo random number
xlabel('Sample Number'); %Assigning X axis of plot
ylabel('Amplitude'); %Assigning Y axis of plot
title('Given Golay Sequence'); %Assigning title to the plot

%Forming Transmitted signal (volts)
t = linspace(0, pulse_width, nsamples); %Sampling time duration from 0 to tpri, sampling at nsamples


tx_signal = zeros(1,nsamples); %Initializing transmitted signal with size of nsamples
for i = 1:subpulse_n %interating on each sub pulse 1 to 20
   for j = 1 : nsamples_perpulse %iterating on pulse to check if lies with the range of current pulse
       if t(j) >= (i-1)*subpulse_width && t(j)<i*subpulse_width %Checking if current sub pulse should be in that range of transmitted signal
           tx_signal(j) = subpulse_code(i); %If True the assigning the value of  sub pulse code at that position
       end
   end
end

%Forming Received Signal (volts)
rx_signal = zeros(1,nsamples); % Initializing received signal for near target


%rx_signal(nd:nd+samples_perpulse-1) = tx_signal(1:samples_perpulse);

for i = nd:nd+nsamples_perpulse-1 %Iterating on position from nd to nd nsamples, to assigning the shifting transmitted signal which will be out received signal
   rx_signal(i) = tx_signal(i-nd+1); %assigning the value of same as transmitted signal but from sample delay
end

%Plotting Transmitted signal 
figure; 
plot(t, tx_signal); 
xlabel('Time (s)'); 
ylabel('Amplitude (Volts)'); 
title('Transmitted Signal'); 
% Plotting Received Signal
figure;
plot(t,rx_signal,'k'); 
xlabel('Time (s)');
ylabel('Amplitude (Volts)');
title('Received Signal for Target at distance 60m'); 

% Time domain matched filtering
mf_t = xcorr(rx_signal,tx_signal); 

% Plotting
figure; %Making new figure to make new plot
plot(linspace(-r_max,r_max,length(mf_t)),abs(mf_t),'k'); %Plotting the output of matched filter
xlabel('Range (m)'); %Assigning X axis of plot
ylabel('Signal after matched filtering (volts)'); %Assigning Y axis of plot
title('Matched filter output(Time domain) wrt Range'); %Title of the plot


%%Frequency Domain Matched filter

tx_f = fft(tx_signal);
rx_f = fft(rx_signal);
sout_f = rx_f.*(conj(tx_f));
sout = ifft(sout_f);
figure; %Making new figure to make new plot
plot(linspace(0,r_max,length(sout)),abs(sout),'k'); %Plotting the output of matched filter
xlabel('Range (m)'); %Assigning X axis of plot
ylabel('Signal after matched filtering (volts)'); %Assigning Y axis of plot
title('Frequency domain Matched Filter output'); %Title of the plot



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Perform 2D range-Doppler processing for a moving target %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dt = 1/BW; %Sampling time which 1/BW, BW is given to us
n_pulse = 2048; %Slow time samples, number of PRI in 1 CPI
m_pulse = 512; %Fast time, number of pulses in 1 PRI

duty_cycle = 0.5; %Duty cycle
new_m = round(m_pulse + (m_pulse * (1-duty_cycle))/duty_cycle); %Finding number samples in 1PRI
tpri = new_m*dt; %Finding PRI, from sampling time

velocity = 50; %Velocity of the target with respect to the radar
target_range = 15; %Initial position of the target with respect to the radar
td = (2*target_range)/c;  %Time delay of received signal
nd = round(td/dt);%Converting time delay to sample delay

wavelength = c/fc;%calculating wavelength of the signal transmitted
fd = (2*velocity)/wavelength; %Doppler frequency, 


tx_mat = xmat(1:m_pulse,1:n_pulse); %Taking desired number of fast time and slow time samples
tx = padding(tx_mat,duty_cycle); %Padding transmitted signal according to the duty cycle.
rx_shift = rangeDelay(tx,nd,m_pulse); %Recieved signal after range delay from the radar,
rx = dopplerDelay(rx_shift,fd,tpri); %from doppler frequency finding doppler delay for each pri
%After creating transmitted and recieved signal, we are finally ready to do
%doppler processing, for that we have made a function rangeDoppler, where
%we are giving transmitted and recieved signal and it will give us 2d range
%doppler map for desired target.
range_doppler = rangeDoppler(tx,rx);



%Now after that we are creating xaxis(velocity axis) and y-axis(range axis)
fdmax = 1/(2*tpri);%finding maximum doppler frequency
dfd = 1/(n_pulse*tpri); %finding step size of doppler axis
%creating doppler axis, from -fdmax/2 to fdmax/2, then we are converting it
%to velocity by using wavelength.
vrange = linspace(-fdmax/2,fdmax/2,size(range_doppler,1)).*(wavelength);
%Now we are finding the yaxis which is range axis
%For that we are finding maximum Range,
Rmax = (c * tpri)/2; 
dr = c/(2*BW);%Range Step size
rrange = (0:1:1024).*(dr); %range axis, 


%Plotting the range_doppler we have calculated above for the target
figure;
imagesc(vrange, rrange,(abs(range_doppler))); %ploting
xlabel('Speed (in m/s)'); %labeling x axis
ylabel('Range (in m)'); %labeling y axis
title('Range-Doppler Map For Single Target'); %giving title
colorbar; %adding colorbar to the plot



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Perform 2D range-Doppler processing for a multiple moving target %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


v = [100,200,150]; %Defining velocities for 3 targets
r = [10,20,30]; %Defining Ranges for 3 targets
tx_mat = xmat(1:m_pulse,1:n_pulse); %Taking required signal
tx = padding(tx_mat,duty_cycle); %Padding transmitted signal according to duty cycle
[m,n] = size(tx); %getting final size after padding
range_doppler3 = zeros(m,n); %Initialzing a empty matrix, for adding range doppler for multiple targets

%Now we are looping over all the targets for finding 2d range dopppler map
for i = 1:3
   velocity = v(i); %getting velocity
   target_range = r(i); %getting range
   td = (2*target_range)/c; %calculating time delay
   nd = round(td/dt); %calculating sample delay

   rx_shift = rangeDelay(tx, nd,m_pulse); %range delay shift
   wavelength = c/fc; %calculating wavelength
  
   fd = (2*velocity)/wavelength; %finding the doppler frequency for the given target
   rx = dopplerDelay(rx_shift,fd,tpri); %calculating doppler delay
   range_doppler = rangeDoppler(tx,rx); %finally calculating rangeDoppler
   
   %Now we are adding the range_doppler for current target to the global
   %variable we have creating which is range_doppler3
   range_doppler3 = range_doppler3+ abs(range_doppler);
end


%Plotting the range_doppler3 we have calculated above for three targets
figure;
imagesc(vrange, rrange,(abs(range_doppler3))); %ploting
xlabel('Speed (in m/s)'); %labeling x axis
ylabel('Range (in m)'); %labeling y axis
title('Range-Doppler Map For Multiple Target'); %giving title
colorbar; %adding colorbar to the plot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Perform cell averaged CFAR-based detection of target parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CFAR

%False alarm probability
p = 5e-4;
%creating object phased
detector = phased.CFARDetector2D('TrainingBandSize',[4,3], ...
    'ThresholdFactor','Auto','GuardBandSize',[2,3], ...
    'ProbabilityFalseAlarm',p,'Method','SOCA','ThresholdOutputPort',true);

%making windows
Ngc = detector.GuardBandSize(2);
Ngr = detector.GuardBandSize(1);
Ntc = detector.TrainingBandSize(2);
Ntr = detector.TrainingBandSize(1);
cutidx = [];
colstart = Ntc + Ngc + 1;
colend = new_m - ( Ntc + Ngc);
rowstart = Ntr + Ngr + 1;
rowend = n_pulse - ( Ntr + Ngr);
for m = colstart:colend
    for n = rowstart:rowend
        cutidx = [cutidx,[n;m]];
    end
end

[dets,th] = detector(abs(range_doppler3),cutidx);

rng_grid = rrange;
dop_grid = vrange;
dopplerIndx = [colstart, colend];
rowIndx = [rowstart, rowend];


helperDetectionsMap(abs(range_doppler3), rng_grid,dop_grid,rangeIndx, dopplerIndx, dets);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function helperDetectionsMap(resp,rng_grid,dop_grid,rangeIndx,dopplerIndx,detections)
   figure
   detectionMap = zeros(size(resp));
   detectionMap(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) = ...
     reshape(double(detections),rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);
   h = imagesc(dop_grid,rng_grid,detectionMap);
   xlabel('Doppler (Hz)'); ylabel('Range (m)'); title('Range Doppler CFAR Detections');
   h.Parent.YDir = 'normal';
end


function y = rangeDoppler(tx,rx)
   
   %This function helps to find the 2d range doppler the target
   %using transmitted and recieved signal

   %Finding the shape of the signal
   [m,n] = size(tx);

   %Initiliazing output matrix
   y = zeros(m,n);
    
   %Now firstly we are doing fft across the rows, i.e, slow time
   for i= 1:m
    y(i,:) = fftshift(fft(rx(i,:)));
   end
    
   %Now we are doing matched filter across the columns, i.e, fast time
   for i = 1:n
       txx = tx(:,i); %transmitted column
       rxx = y(:,i); %recieved column
       y(:,i) = ifft(fft(rxx).*conj(fft(txx)))./1024; % doing matched filter in frequency domain
   end
end


function y = dopplerDelay(x, fd,tpri)

   %In this function we are adding doppler delay to the signal 
   % to find the received signal
   
   %firstly initializing the dimensions 
   [m,n] = size(x);
   %Initializing zeros matrix to store the doppler delay output
   y = zeros(m,n);
    
   %Loop over all the columns to find the doppler delay as per the equation
   for i = 1:n
       y(:,i) = x(:,i).*(exp(1j*2*pi*fd*(i-1)*tpri));
   end
end


function y = rangeDelay(x,nd,m_pulse)
   %In this function we are adding range delay to the signal 
   % to find the received signal
   
   %firstly initializing the dimensions 
   [m,n] = size(x);
   %Initializing zeros matrix to store the range  delay output
   y = zeros(m,n);


   for i = 1:n

       %target delay for the first column
       ndx = nd;

       for j = 1:m_pulse
           y(ndx,i) = x(j,i); %assining value from according to the delayy
           ndx=ndx+1; %increaing delay pointer
       end
   end
end

function y = rangeDelay_velocity(x,target_range,tpri,dt,velocity,c,m_pulse)
   
   %This is the function we have not used but, we have created to
   %as per the case we have given that if we have a moving target
   %the target delay for each PRI would be different, i.e, the target
   %range from the radar would be different which will be equal to 
   % r = ro + nvtpri, where ro is the initial target position and
   %v is the velocity of the target and n is the nth pri, starting from
   % 0 to N-1, according for each pri the time delay and sample delay is
   % calculated
   
   
   %firstly initializing the dimensions 
   [m,n] = size(x);
   %Initializing zeros matrix to store the range  delay output
   y = zeros(m,n);
    

   %Looping over all columns and finding the range delay for each column

   for i = 1:n
       %Calculating the new target position from the radar
       target_rng = target_range + velocity*tpri*(i-1);
       %finding the time delay
       td = (2*target_rng)/c;
       %finding the sample delay
       nd = round(td/dt);
       ndx = nd;
       for j = 1:m_pulse
           y(ndx,i) = x(j,i); %assigning values according to the delay
           ndx=ndx+1; %increaing the pointer
       end
   end
end

function y = padding(x, duty_cycle)
   
   %getting the dimension of the signal
   [m,n] = size(x);
   %calculating new dimension after apply duty cycle padding
   new_m = round(m + (m * (1-duty_cycle))/duty_cycle);
   %initializing matrix with new dimension
   y = zeros(new_m,n);
   %adding the signal to newly created matrix, so that it looks it is
   %padded according to the duty cycle
   y(1:m,:) = x;
end
