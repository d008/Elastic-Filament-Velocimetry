
% Find a serial port object.
obj2 = instrfind('Type', 'serial', 'Port', 'COM9', 'Tag', '')

% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(obj2)
    obj2 = serial('COM9');
else
    fclose(obj2);
    obj2 = obj2(1);
end

fopen(obj2);

% Flush the data in the input buffer.
flushinput(obj2);
% Communicating with instrument object, obj2.
data10 = fscanf(obj2);
MASS = strsplit(data10);MASS = str2num(MASS{2});
MASS
fclose(obj2);