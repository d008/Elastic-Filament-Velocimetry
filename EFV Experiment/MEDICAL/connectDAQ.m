function [s] = connectDAQ(N,RATE)
%Open Session
s = daq.createSession('ni');
%EFV Channel
efv_Channel = addAnalogInputChannel(s,'Dev1','ai0','Voltage');

%TNSTAP Channel
%tnstap_Channel = addAnalogInputChannel(s,'Dev3','ai7','Voltage');

%Set Default Rate
s.Rate = RATE;
%Set Default Duration of Time
%s.DurationInSeconds = N;
s.IsContinuous = true;
s.NotifyWhenDataAvailableExceeds = s.Rate*N;

fprintf('\nDAQ Connected \n');

end