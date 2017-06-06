%cd(uigetdir)
listv = dir('*.mat');
nPts = length(listv);
MASS = zeros(nPts,2);
VOLT = zeros(nPts,3);

for ind = 1:nPts
    load(listv(ind).name)
    raw0 = strcat(filename,'0.bin');
    raw = strcat(filename,'.bin');
    
    %ZERO DATA
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
    
    MASS(ind,1) = median(diff(smooth(masses(1:end,2),11,'moving'))./diff(masses(1:end,1))...
        *3600/testParam.density);
    MASS(ind,2) = std(diff(smooth(masses(1:end,2),11,'moving'))./diff(masses(1:end,1))...
        *3600/testParam.density);
    VOLT(ind,1) = mean(EFV_V0);
    VOLT(ind,2) = mean(EFV_V');
    VOLT(ind,3) = std(EFV_V'-EFV_V0);
    
end
figure(1)
clf
%errorbar(MASS(:,1)*(55*15)/(pi*0.002^2)/3600/1000/1000,VOLT(:,2)-VOLT(:,1),VOLT(:,3),'o')
plot(MASS(:,1),(VOLT(:,2)-VOLT(:,1)),'o')

%hold on
%herrorbar(MASS(:,1),VOLT(:,1),MASS(:,2),'o')
xlabel('Flow Rate ml/hr')
%xlabel('m/s')

ylabel('Volts')