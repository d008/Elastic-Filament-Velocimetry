function [m] = getWeight(scale)
flushinput(scale);
m = fscanf(scale);
m = strsplit(m);
m = str2double(m{2});
end