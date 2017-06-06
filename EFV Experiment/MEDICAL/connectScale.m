function [SCALE] = connectScale()
% Create the serial port object if it does not exist
% otherwise use the object that was found.
%USB PORT ON - SBI , AUTO-PRINT w/o stability
SCALE = instrfind('Type', 'serial', 'Port', 'COM3', 'Tag', '');
if isempty(SCALE)
    SCALE = serial('COM3');
else
    fclose(SCALE);
    SCALE = SCALE(1);
end
% Connect to instrument object, LCR.
fopen(SCALE);
% Flush the data in the input buffer.
flushinput(SCALE);
fprintf('\n Scale Connected\n');
end