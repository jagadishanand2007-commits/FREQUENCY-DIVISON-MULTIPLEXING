# FREQUENCY-DIVISON-MULTIPLEXING
### Aim:
Study of TDM pulse amplitude modulation/ demodulation with transmitter block (clock) and channel identification
information linked directly to the receivers.

### Apparatus Required:
1) Experimental kit DCL-02
2) Connecting chord
3) Power supply
4) 20 MHz dual trace oscilloscope.
   
### Theory:
In PAM, PPM the pulse is present for a short duration and for most of the time between the two pulses no signal is
present. This free space between the pulses can be occupied by pulses from other channels. This is known as Time
Division Multiplexing. Thus, time division multiplexing makes maximum utilization of the transmission channel.
Each channel to be transmitted is passed through the low pass filter. The outputs of the low pass filters are connected
to the rotating sampling switch (or) commutator.
It takes the sample from each channel per revolution and rotates at the rate of f s. Thus the sampling frequency
becomes fs the single signal composed due to multiplexing of input channels. These channels signals are then passed
through low pass reconstruction filters. If the highest signal frequency present in all the channels is fm, then by
sampling theorem, the sampling frequency fs must be such that fs≥2fm. Therefore, the time space between successive
samples from any one input will be Ts=1/fs, and Ts 1/2fm.

### Procedure:
1) Refer to the block diagram and carry out the following connections and switch setting.
2) Connect power supply in proper polarity DCL-02 and s it switch it on.
3) Connect 250hz, 500hz, 1khz and 2khz sine wave signals from the function generator to the multiplexer input
channels CH0, CH1,CH3 by means of connecting chords.
4) Connect the multiplexer output txd of the transmitted section to the de multiplexer input rxd of the receiver
section.
5) Connect the output of the receiver section CH0, CH1, CH2,CH3 to the IN0, IN1, IN2. IN3 of the filter section.
6) Connect the sampling clock TX CLK channel identification clock TXSYNC of the transmitter section to the
corresponding RXCLK, RXSYNC of the receiver section respectively.
7) Set the amplitude of the input sine wave desired.
8) Take the observations as mentioned below.
   
### Kit Diagram:
<img width="1245" height="700" alt="image" src="https://github.com/user-attachments/assets/6487a1dc-1f93-4ffc-befc-60cf907c7192" />

### Model Graph:
<img width="615" height="808" alt="image" src="https://github.com/user-attachments/assets/1d1f8d48-f08b-40fa-900c-0cb962a7e316" />

### CODE
```
clc; clear; close;

fs = 100000;
T = 0.02;             
t = 0:1/fs:T-1/fs;   
N = length(t);          

m1 = sin(2*%pi*100*t);
m2 = sin(2*%pi*200*t);
m3 = sin(2*%pi*300*t);
m4 = sin(2*%pi*400*t);
m5 = sin(2*%pi*500*t);
m6 = sin(2*%pi*600*t);

fc = [6000 12000 18000 24000 30000 36000]; 

s1 = m1 .* cos(2*%pi*fc(1)*t);
s2 = m2 .* cos(2*%pi*fc(2)*t);
s3 = m3 .* cos(2*%pi*fc(3)*t);
s4 = m4 .* cos(2*%pi*fc(4)*t);
s5 = m5 .* cos(2*%pi*fc(5)*t);
s6 = m6 .* cos(2*%pi*fc(6)*t);

fdm = s1 + s2 + s3 + s4 + s5 + s6;

F = fft(fdm);
freq = (0:N-1)*(fs/N);

fshift = freq;
idx = find(freq > fs/2);
if ~isempty(idx) then
    fshift(idx) = freq(idx) - fs;
end

msg_bw = 1500;
halfband = msg_bw + 300; 

recovered = zeros(6, N);

for k = 1:6
    BPmask = zeros(1, N);
    BPmask = (abs(fshift - fc(k)) <= halfband) | (abs(fshift + fc(k)) <= halfband);
    Fbp = F .* BPmask;
    xbp = fft(Fbp, -1) / N;
    down = 2 .* xbp .* cos(2*%pi*fc(k)*t);
    D = fft(down);
    LPmask = (abs(fshift) <= msg_bw);
    Dlp = D .* LPmask;
    base = fft(Dlp, -1) / N;
    recovered(k, :) = real(base);
end

figure(1);
subplot(3,2,1); plot(t, m1); title("Original m1 (100 Hz)");
subplot(3,2,2); plot(t, m2); title("Original m2 (200 Hz)");
subplot(3,2,3); plot(t, m3); title("Original m3 (300 Hz)");
subplot(3,2,4); plot(t, m4); title("Original m4 (400 Hz)");
subplot(3,2,5); plot(t, m5); title("Original m5 (500 Hz)");
subplot(3,2,6); plot(t, m6); title("Original m6 (600 Hz)");
xlabel("Time (s)");

figure(2);
plot(t, fdm);
title("FDM Composite Signal (time domain)");
xlabel("Time (s)");

figure(3);
subplot(3,2,1); plot(t, recovered(1,:)); title("Recovered m1");
subplot(3,2,2); plot(t, recovered(2,:)); title("Recovered m2");
subplot(3,2,3); plot(t, recovered(3,:)); title("Recovered m3");
subplot(3,2,4); plot(t, recovered(4,:)); title("Recovered m4");
subplot(3,2,5); plot(t, recovered(5,:)); title("Recovered m5");
subplot(3,2,6); plot(t, recovered(6,:)); title("Recovered m6");
xlabel("Time (s)");

disp("✅ FFT-based demux complete. Compare originals vs recovered visually or compute MSE:");
for k=1:6
    mse = sum((recovered(k,:) - eval("m"+string(k))).^2)/N;
    disp("Channel "+string(k)+" MSE = "+string(mse));
end

```

### OUTPUT WAVEFORM:

<img width="748" height="604" alt="image" src="https://github.com/user-attachments/assets/9ef25734-95bd-4b82-a6ae-c01ee5993fe9" />

<img width="696" height="566" alt="image" src="https://github.com/user-attachments/assets/fcfeb166-a469-486b-8d65-79fd91b2ecf6" />

<img width="739" height="572" alt="image" src="https://github.com/user-attachments/assets/69de2460-9d98-4f91-b2f1-400da0ced430" />


### TABULATION:
<img width="720" height="1280" alt="image" src="https://github.com/user-attachments/assets/02256583-ea72-4c37-a11f-5b361e417cfe" />

<img width="1280" height="720" alt="image" src="https://github.com/user-attachments/assets/e4abed9a-2ad2-45f1-9595-7849e1f708e5" />



### Result:
Thus the time division multiplexing is done experimentally and output is verified.
