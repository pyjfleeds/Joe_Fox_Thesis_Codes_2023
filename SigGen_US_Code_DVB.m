clear all
close all
clc
%Instrument Connection

% Find a serial port object.
obj1 = instrfind('Type', 'serial', 'Port', 'COM3', 'Tag', '');

% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = serial('COM3');
else
    fclose(obj1);
    obj1 = obj1(1);
end
fopen(obj1);
fprintf(obj1, 'OUTPUT OFF'); % Turn off Agilent
%% Set frequency, amplitude in mV, Duty Cycle, PRF and Duration
Frequency = 2.25e6 ; % Frequency in Hz
MI = 0.6;
PNP = 1e3*(MI * sqrt(Frequency*1e-6)); %PNP in kPa
Amplitude = 1e-3*((PNP/2.43) - 21.68); %  Sig Gen Amplitude in V calibrate for V323-SM
DutyCycle = 1 ; % Duty Cycle in Percent
PRF = 1e3; % Pulse Repetition Frequency in Hz
Duration = 5; % Total Duration in Seconds

AmplifiedAmp = (Amplitude * 446.6835921509635) / 1000 ; % + 53 dB Amplified Voltage in V
BurstPeriod = 1 / PRF  ;
BurstCount = ((DutyCycle)/100) * Frequency / PRF ;
FrequencyString = sprintf( 'FREQ %d', Frequency) ; 
AmplitudeString = sprintf('AMPL %d', Amplitude);
CountString = sprintf('BSTCOUNT %d', BurstCount);
PeriodString = sprintf('BSTPER %d', BurstPeriod);
fprintf(obj1, 'BST NCYC');
fprintf(obj1, FrequencyString); % set frequency
fprintf(obj1, 'WAVE SINE'); % set wave type
fprintf(obj1 , 'AMPUNIT VPP'); %set amplitude units
fprintf(obj1 ,  AmplitudeString) % set amplitude using units above
fprintf(obj1, CountString) %set burst count
fprintf(obj1, PeriodString); % Set Burst Period

%% Turn on, pause, turn off
if Amplitude > 0.4 || DutyCycle > 1
    display('Input Voltage > 400 mV! or Duty Cycle > 1 %')
    fprintf(obj1, 'OUTPUT OFF');
else
    input('Ready?');
    for j=1:Duration-1
    display(num2str(j));
    pause(1);
    end
    display('ON!');
    fprintf(obj1, 'OUTPUT ON'); %turn on,
    
    pause(Duration);
   
    fprintf(obj1, 'OUTPUT OFF');
end

 display('Done!');
 