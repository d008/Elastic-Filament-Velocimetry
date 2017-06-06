%% HEADER
clear all
clc
N = 240;    %  Number of points/seconds of data collection
density = 1; %g/mL of fluids

data = zeros(N,3);
formatString = ' Point  : %d\t\n Time(s): %0.1f\t\n Mass(g): %0.3f\t\n R(Ohms): %0.3f\n';
%formatString = '\t%d \t|\t %0.1f s\t|\t %0.3f g\t|\t %0.3f Ohms \n';

fprintf('\nConnecting to Instruments...\n')

% Find a VISA-USB object.
LCR = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0957::0x0909::MY54202264::0::INSTR', 'Tag', '');
% Find a serial port object.
SCALE = instrfind('Type', 'serial', 'Port', 'COM9', 'Tag', '');

% Create the VISA-USB object if it does not exist
% otherwise use the object that was found.
if isempty(LCR)
    LCR = visa('AGILENT', 'USB0::0x0957::0x0909::MY54202264::0::INSTR');
else
    fclose(LCR);
    LCR = LCR(1);
end
% Connect to instrument object, LCR.
fopen(LCR);

% Create the serial port object if it does not exist
% otherwise use the object that was found.
%USB PORT ON - SBI , AUTO-PRINT w/o stability
if isempty(SCALE)
    SCALE = serial('COM9');
else
    fclose(SCALE);
    SCALE = SCALE(1);
end
% Connect to instrument object, LCR.
fopen(SCALE);

fprintf(LCR, ':FETCh:SMONitor:IAC? ');
IAC = str2double(fscanf(LCR));
flushinput(LCR);

% Flush the data in the input buffer.
flushinput(SCALE);
% Communicating with instrument object, SCALE.
MASS = fscanf(SCALE);MASS = strsplit(MASS);MASS = str2double(MASS{2});
% Communicating with instrument object, LCR :Resitance

%Filename save
filename = input('\nDatafile name: ','s');

%Pause for testing
input('\nPress Enter to Start Test: ');
clc
fprintf('Testing Commenced:\n');

figure(1)
clf
hold on

subplot(2,1,1)
xlabel('t (sec)')
ylabel('Mass Flow Rate (g/sec)')
subplot(2,1,2)
xlabel('t (sec)')
ylabel('Resistance (Ohms)')

flushinput(LCR);
fprintf(LCR, 'FORM:DATA ASCii');
fprintf(LCR, ':FETCh:IMPedance:FORMatted? ');
R = strsplit(fscanf(LCR),',');R = R{2}; R = str2double(R);
data(1,:) = [0, MASS,R];
fprintf(formatString,1,0, MASS,R)
tic
pause(1)
R0 = R;
M0 = MASS;


for i = 2:N
    % Flush the data in the input buffer.
    flushinput(SCALE);
    
    % Communicating with instrument object, SCALE.
    MASS = strsplit(fscanf(SCALE));MASS = str2double(MASS{2});
    % Communicating with instrument object, LCR :Resitance
    flushinput(LCR);
    fprintf(LCR, 'FORM:DATA ASCii');
    fprintf(LCR, ':FETCh:IMPedance:FORMatted? ');
    R = strsplit(fscanf(LCR),','); R = str2double(R{2});
    
    data(i,:) = [toc, MASS,R];
    clc
    fprintf('Testing Commenced:\n');
    fprintf(formatString,i,toc, MASS,R)
    subplot(2,1,1)
    if (i > 2)
    plot(data(3:i,1),(data(3:i,2)*3/2-2*data(2:i-1,2)+data(1:i-2,2)/2)./...
        (data(3:i,1)-data(2:i-1,1)),'-o')
    end
    ylabel('Mass Flow Rate (g/sec)')
    subplot(2,1,2)
    plot(data(1:i,1),data(1:i,3),'-o')
    ylabel('R (Ohm)')
    xlabel('Time (sec)')
    drawnow
    
end

R = data(:,3); M = data(:,2); t = data(:,1);
mfr = fit(t,M,'poly1');
fprintf('\nVol. Flow Rate = %0.2f g/hr \n',mfr.p1*3600./density)
fprintf('Medain Resistance = %0.3f Ohms\n',median(R))
fprintf('Mean Resistance = %0.3f Ohms\n',mean(R))
fprintf('Std. Resistance = %0.3f Ohms\n',std(R))
%std((data(:,3)))
% Disconnect from instrument onj2, LCR.
fclose(LCR);
fclose(SCALE);
clear formatString i MASS LCR SCALE N
save(filename,'IAC','R','M','t')
%fprintf(E4980AL, sprintf(':CURRent:LEVel %g', 0.0001));
%fprintf(E4980AL, sprintf(':FREQuency:CW %g', 300000.0));
