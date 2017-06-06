function  stashData(src, event,fid,efv0)
subplot(2,1,1)

%plot(event.TimeStamps, event.Data)
%title((mean(event.Data)-efv0)./aefv)
errorbar(mean(event.TimeStamps), (mean(event.Data)-efv0),(std(event.Data)),'bx')
temp = [event.TimeStamps, event.Data];
fwrite(fid,temp','double');
xlabel('Time (sec)');
ylabel('V')
end
