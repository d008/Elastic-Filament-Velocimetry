%clear all
clf
clc
%try
%    load('tp.mat')
%catch
testParam = ([]);  
testParam.sampleTime = 10;
testParam.pauseTime = 1;
testParam.zeroTime = 5;
testParam.Rate = 20000;
testParam.density = 1.1;

% testParam.RtopTNSTAP = 51000;
% testParam.RbalTNSTAP = 74;
% testParam.GainTNSTAP = 4000;
% 
% testParam.RtopEFV = 51000;
% testParam.RbalEFV = 107;
% testParam.GainEFV = 4080;

testParam.Vsupply = 10;
testParam.InputGain = 40;
testParam.InputGain = 3;
testParam.Channel = 2;

%end
%%%Rx = @(V0,Vs,G,R1,R2) (-V0./(Vs.*G)+R2./(R1+R2))*R1./(1-(-V0./(Vs.*G)+R2./(R1+R2)));
%%%alpha = @(V0,Vs,G,R1,R2) Vs*G.*Rx(V0,Vs,G,R1,R2)*R1./(Rx(V0,Vs,G,R1,R2)+R1)^2;
% I = testParam.Vsupply/testParam.RtopEFV
% POWER = I^2*testParam.RbalEFV
alpha_approx =@(V0,Vs,G,R1,R2) Vs*G*R2./R1;



N = testParam.sampleTime;    %  Number of points/seconds of data collection
N0 = testParam.zeroTime;    %  Number of points/seconds of data collection

masses = zeros(N,2);
masses0 = zeros(N0,2);
formatString = ' Point  : %d\t\n Time(s): %0.1f\t\n Mass(g): %0.3f\t\n R(Ohms): %0.3f\n';
fprintf('\nConnecting to Instruments...\n')

SCALE = connectScale()
%s = connectDAQ(1,testParam.Rate)


initValue = s.inputSingleScan;
EFV_V0 = initValue(1);

% Rx(EFV_V0,testParam.Vsupply, testParam.GainEFV,...
%     testParam.RtopEFV,testParam.RbalEFV)
% 
% Rx(TNSTAP_V0,testParam.Vsupply, testParam.GainTNSTAP,...
%     testParam.RtopTNSTAP,testParam.RbalTNSTAP)

% alphaEFV = alpha_approx(EFV_V0,testParam.Vsupply, testParam.GainEFV,...
%     testParam.RtopEFV,testParam.RbalEFV)
% 
% alphaTNSTAP = alpha_approx(TNSTAP_V0,testParam.Vsupply, testParam.GainTNSTAP,...
%     testParam.RtopTNSTAP,testParam.RbalTNSTAP)
%% ZEROING CYCLE

%Filename save - zero raw data
filename = input('\nDatafile name: \n','s');
raw0 = strcat(filename,'0.bin');
fid1 = fopen(raw0,'w');

%Create Listener
lh = addlistener(s,'DataAvailable', ...
    @(src,event) stashData(src,event,fid1,EFV_V0));

%Pause for testing
input('\nPress Enter to Start Zero: \n');
s.startBackground()
clc
fprintf('Zeroing Commenced:\n');
tic

figure(1)
clf
subplot(2,1,1)
hold on
subplot(2,1,2)
title('ZERO')
for i = 1:N0
    masses0(i,:)= [toc,getWeight(SCALE)];
    subplot(2,1,2)
    plot(masses0(1:i,1),masses0(1:i,2),'b-o');
    drawnow
end
beep
s.stop()
delete(lh);
fclose(fid1);
%%
%ZERO DATA
fid2 = fopen(raw0,'r');
[V,count] = fread(fid2,[2,inf],'double');
fclose(fid2);

EFV_V0 = mean(V(2,:))

% alphaEFV = alpha(EFV_V0,testParam.Vsupply, testParam.GainEFV,...
%     testParam.RtopEFV,testParam.RbalEFV);
% 
% alphaTNSTAP = alpha(TNSTAP_V0,testParam.Vsupply, testParam.GainTNSTAP,...
%     testParam.RtopTNSTAP,testParam.RbalTNSTAP);
%% TESTING CYCLE
%Filename save
s = connectDAQ(1,testParam.Rate);
beep
input('\nSTART PUMP & Press Enter to Start Test: \n');

raw = strcat(filename,'.bin');
fid1 = fopen(raw,'w');

%Create Listener
lh2 = addlistener(s,'DataAvailable', ...
    @(src,event) stashData(src,event,fid1,EFV_V0));
%Pause for testing
h = waitbar(0,sprintf('Pausing for %d secs.\n', testParam.pauseTime));
for i =1:testParam.pauseTime
    waitbar(i/testParam.pauseTime,h);
    pause(1);
end
beep
close(h)

clc
fprintf('Testing Commenced:\n');
s.startBackground()

tic
figure(1)
clf
subplot(2,1,1)
legend('EFV','TNSTAP','location','bestoutside');
hold on
subplot(2,1,2)
xlabel('Time (sec)');
ylabel('Mass (g)')
title('TEST')

for i = 1:N
    masses(i,:)= [toc,getWeight(SCALE)];
    subplot(2,1,2)
    plot(masses(2:i,1),...
        diff(smooth(masses(1:i,2),11,'moving'))...
        ./diff(masses(1:i,1))*3600/testParam.density,'b-o');
    drawnow
end
beep
s.stop()
fclose(fid1)
delete(lh2);
save(filename,'masses','masses0','filename','testParam','raw','raw0')
fprintf('DONE:\n');
fclose(SCALE);
%% Collect Results
%ZERO DATA
raw0 = strcat(filename,'0.bin');
raw = strcat(filename,'.bin');

fid2 = fopen(raw0,'r');
[V,count] = fread(fid2,[2,inf],'double');
fclose(fid2);

EFV_V0 = mean(V(2,:))

%TEST DATA
fid2 = fopen(raw,'r');
[V,count] = fread(fid2,[2,inf],'double');
fclose(fid2);
time = V(1,:);
EFV_V = V(2,:);


sf = fit( time', (EFV_V'-EFV_V0),'cubicinterp');
sf2 = fit(masses(2:end,1),diff(smooth(masses(1:end,2),11,'moving'))./diff(masses(1:end,1))...
    *3600/testParam.density, 'cubicinterp');

PROCESS
% TS = 0:testParam.sampleTime;
% x1 = sf2(TS);
% y1 = sf(TS);
% figure(1)
% clf
% plot(x1,y1,'bo')
% sf3= fit(x1,y1, 'poly1')
% hold on
% plot(0:120,sf3(0:120))
% mfr = fit(masses(:,1),masses(:,2),'poly1');
% xlabel('Flow Rate ml/hr')
% ylabel('Volts')
% 
% print FLOWvVolts_ETHYLENE_GLYCOL -dpng
% figure(2)
% clf
% yyaxis left
% plot(TS,y1)
% xlabel('Time(s)')
% ylabel('AU')
% 
% yyaxis right
% plot(TS,x1)
% 
% ylabel('Flow Rate ml/hr')
% print TIMESERIES -dpng