% Find a VISA-USB object.
clear all
obj1 = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0957::0x0909::MY54202264::0::INSTR', 'Tag', '')

% Create the VISA-USB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = visa('AGILENT', 'USB0::0x0957::0x0909::MY54202264::0::INSTR');
else
    fclose(obj1);
    obj1 = obj1(1);
end



% Connect to instrument object, obj1.
fopen(obj1);

% Communicating with instrument object, obj1 :Resitance
fprintf(obj1, 'FORM:DATA ASCii');
fprintf(obj1, ':FETCh:IMPedance:FORMatted? ');
data1 = fscanf(obj1)
data1 = strsplit(fscanf(obj1),',')
R = data1{1}; R = str2num(R)

flushinput(obj1);
fprintf(obj1, ':FETCh:IMPedance:CORRected? ');
data2 = strsplit(fscanf(obj1),',')
R = data2{1}; R = str2num(R)
flushinput(obj1);

% % Communicating with instrument object, obj1.: Current
% fprintf(obj1, ':FETCh:SMONitor:IAC? ');
% data2 = fscanf(obj1);
% IAC = str2num(data2)
% flushinput(obj1);
%pause(2)


% Disconnect from instrument object, obj1.
fclose(obj1);