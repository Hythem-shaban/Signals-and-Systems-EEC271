% Transmitter
clear;
clc;
% Import Audio File
[x, fs]=audioread('Imagine_Dragons_-_Thunder.wav');

% Period
time_initial  = input('Enter Starting time in audio file : ');   
time_final = input('Enter Ending time in audio file : ');         

% If statement used to check if input time_final less than or equal
% time_final and if so the audio file will be played from zero to the end

if( time_initial >= time_final)
    time_initial = 0;
    time_final = length(x)/fs;
end

% Using only time of input Audio file
x=x(1+(time_initial)*fs:(time_final)*fs,1).';

% making time initial to zero so that delta can be applied
        time_final = time_final - time_initial;
        time_initial = 0;
    
% Playing time of Audio file
  sound(x,fs);

% Time domain Representation
t=linspace(time_initial,time_final,fs*(time_final-time_initial));

figure
plot(t,x);
xlabel('Time in sec');
ylabel('Input Signal in time domain x')
title('Input Signal Representation in time domain');

% Frequency domain Representaion
X= fftshift(fft(x));
Fvec = linspace(-fs/2,fs/2,length(X));

% Magnitude Specturm Ploting
figure
plot(Fvec,abs(X));
xlabel('Frequency in Hz');
ylabel('Magnitude of input signal abs(X)')
title('Magnitude Spectrum');

% Phase Specturm Ploting
figure
plot(Fvec,angle(X));
xlabel('Frequency in Hz');
ylabel('phase of input signal angle(X)')
title('Phase Spectrum');
%%

% Channel

% option is variable used to store the choice of user from menu function
% which return value corresponding to the position of option
% if Option 1 was chosen , option = 1, if option 2 was chosen, option = 2 and etc
option = menu('Choose Channel Impusle resonose','Option 1','Option 2','Option 3','Option 4');

switch (option)
    
    case 1 
      % Impulse response of the channel (delta function)  
      % dirac function return signal dc signal from 0 to inf
      % then index is used to check where is the inf in h
      % then this index in h is then assigned to be equal 1
      % then h become the delta function
      
     h = delta(time_initial,time_final,0,t);
      
      % Impulse Response of channel Plotting 
        figure;
        stem(t, h);
        xlabel('Time in sec');
        ylabel('Impulse Response h');
        title(['Impulse Response of channel ' num2str(option)]);
        
    case 2    
      % Impulse response of the channel (exponential function)
      h = exp(-2*pi*5000*t);
      
      % Impulse Response of channel Plotting 
        figure;
        stem(t, h);
        xlabel('Time in sec');
        ylabel('Impulse Response h');
        title(['Impulse Response of channel ' num2str(option)]);
      
    case 3
      % Impulse response of the channel (exponential function)
      h = exp(-2*pi*1000*t);
      
      % Impulse Response of channel Plotting 
        figure;
        stem(t, h);
        xlabel('Time in sec');
        ylabel('Impulse Response h');
        title(['Impulse Response of channel ' num2str(option)]);

    case 4
      % Impulse response of the channel (two delta functions)   
      h = 2*delta(time_initial,time_final,0,t) + 0.5*delta(time_initial,time_final,1,t) ;
    
      % Impulse Response of channel Plotting 
        figure;
        stem(t, h);
        xlabel('Time in sec');
        ylabel('Impulse Response h');
        title(['Impulse Response of channel ' num2str(option)]);    
end

% Performing convolution between input signal and impulse respose to get
% y signal by using conv function
y = conv(x,h);

% t_new is calculated as number of samples of y will be equal to that of x
% plus that of h plus one so samples(t_new) = 2*samples(t)+1 
% but the rest of y after time_final equal to zero so we shorten y to
% length of x so that signal transimtted have same duration as orignal
% input
if(option == 4)
    y = y(1:length(x)+fs);

else
     y = y(1:length(x));   
end
t_new = linspace(0, length(y)/fs, length(y));

% Ploting Signal after passing through channel
figure
plot(t_new, y);
xlabel('Time in sec');
ylabel('Output Signal of the Channel y');
title(['Output Signal of the Channel ',num2str(option)]);

% Frequency domain Representaion
Y= fftshift(fft(y));
Fvec = linspace(-fs/2,fs/2,length(Y));

% Magnitude Specturm Ploting
figure
plot(Fvec,abs(Y));
xlabel('Frequency in Hz');
ylabel('Magnitude of input signal abs(X)')
title(['Magnitude Spectrum Output Signal of the Channel ',num2str(option)]);

% Phase Specturm Ploting
figure
plot(Fvec,angle(Y));
xlabel('Frequency in Hz');
ylabel('phase of input signal angle(X)')
title(['Phase Spectrum Output Signal of the Channel ',num2str(option)]);
%%

% Noise

% Adding Noise to output signal of the channel
Noise = menu('Add Noise','Yes','No');

switch (Noise)

    case 1
       % Sigma Input Section
       sigma = input('Enter value of sigma = ');

       % Noise Decleration
       Z = sigma*randn(1,length(y));

       y_N = y + Z;
        
    case 2
        
        y_N = y;
        
end

% Play Sound
sound(y_N,fs);

% Time domain Representation
figure
plot(t_new,y_N);
xlabel('Time in sec');
ylabel('Output Signal after adding Noise y_N');
title('Output Signal Representaion in Time domain');

% Frequency domain Representaion
Y_N=fftshift(fft(y_N));
Fvec_new = linspace(-fs/2,fs/2,length(Y_N));

% Magnitude spectrum ploting
figure
plot(Fvec_new,abs(Y_N));
xlabel('Frequency in Hz');
ylabel('Magnitude of output Signal after adding Noise Y_N');
title('Magnitude Spectrum after adding noise');

% Phase spectrum ploting
figure
plot(Fvec_new,angle(Y_N));
xlabel('Frequency in Hz');
ylabel('Phase of output Signal after adding Noise Y_N');
title('Phase Spectrum after adding noise');

%%

% Receiver

% y_N is signal received after noise is added
% y_R is received signal at Receiver in time domain
y_R = y_N ;

% y_R is converted to frequency domain to apply filtering
Y_R = fftshift(fft(y_R));

% Fvec_new represents all frequencies in signal Y_R
% Note: it is Fvec_new as there is Fvec previously used in reperesenting
% input signal in frequency domain at the transmitter 
Fvec_new = linspace(-fs/2,fs/2,length(Y_R));

% sample_per_hertz is calculated which will be used in filtering the signal
sample_per_hertz = (length(y_N)/fs);

% the filter will allow signal of frequency ranging from
% -3.4khz to 3.4khz while prevent the other frequencies from passing so
% samples from start to -3.4 khz and samples from 3.4 khz to end become zeros
% round is used to eliminate the probability of existing fraction index as
% that will result in error.
Y_R (1:round(sample_per_hertz*(fs/2 - 3.4*10^3))) = 0;
Y_R (round(sample_per_hertz*(fs/2 + 3.4*10^3))+ 1:end)=0;

% Representing Received Signal in frequency domain after filtering

% Magnitude Spectrum after filtering
figure
plot(Fvec_new,abs(Y_R));
title('Magnitude Spectrum after filtering')
xlabel('Frequency in (HZ)');
ylabel('Received signal Magnitude  abs(Y_R)');

% Phase Specturm after filtering
figure
plot(Fvec_new,angle(Y_R));
title('phase Spectrum after filtering')
xlabel('Frequency in (HZ)');
ylabel('Received signal Phase angle(Y_R)');

% Returning signal back to time domain
y_R = real(ifft(ifftshift(Y_R)));

% Playing Received signal after filtering
% as it is noticed only the voice is heard as this filter pass frequency of
% voice which is from -3.4khz to 3.4khz and reject the other frequecnies
sound(y_R,fs);

% Representing Received Signal in time domain after filtering
figure
plot(t_new,y_R);
title('Received Signal in time domain after filtering')
xlabel('time in (sec)');
ylabel('Received signal y_R');
