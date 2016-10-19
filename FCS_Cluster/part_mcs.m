function [mcs lifetimedata head res1] = part_mcs(ff,package_size,tg,statusfield)
% This duplicates code from part_lifetimedata. The lifetimedata has
% to be calculated, when the mcs data was loaded and the calculation
% is cheap relative to the cost for reading and calculation mcs.
tg1 = tg.tg1;
head = ht3v2read_cl(ff, [0 0]);
count_packages = ceil(head.Records / package_size);
Timeunit=1/head.SyncRate;
use_tg1 = tg1 ./ Timeunit;

mcs = zeros(head.StopAfter, 6);
mcs(:,1) = (1:head.StopAfter)' / 1000;
lifetimedata = zeros(65536,6);
lifetimedata(:,1) = (1:65536) * head.Resolution;
overcount = 0;

for i=1:count_packages
    if exist('statusfield', 'var')
        set(statusfield, 'String', sprintf('Importing: %03i / %03i', i, count_packages));
        drawnow;
    end
    fprintf('i: %i - max: %i\n', i, count_packages);
    si = (i-1)*package_size+1;
    ei = i*package_size;
    if ei > head.Records
        ei = head.Records;
    end

    [y, tcspc, chan, markers, num, head] = ht3v2read_cl(ff, [si ei-si+1]);
    chan = chan + 1;

    old_overcount = overcount;
    z = (overcount + y) * 1000 / head.SyncRate;
    overcount = y(size(y,1),1) + old_overcount;
    res1.time(i) = y(size(y,1),1) / head.SyncRate;
   
    % head.StopAfter is in ms - so we get a histogramm with a ms resolution
    mcs(:,2) = mcs(:,2) + histc(z(chan==1), 1:head.StopAfter);
    mcs(:,3) = mcs(:,3) + histc(z(chan==2), 1:head.StopAfter);
    mcs(:,4) = mcs(:,4) + histc(z(chan==3), 1:head.StopAfter);
    mcs(:,5) = mcs(:,5) + histc(z(chan==4), 1:head.StopAfter);
    mcs(:,6) = mcs(:,6) + histc(z(chan>0 & chan<5), 1:head.StopAfter);
    usable_data = z >= use_tg1(1) & z < use_tg1(2);
    lifetimedata(:,2) = lifetimedata(:,2) + histc(tcspc(chan==1 & usable_data), 1:65536);
    lifetimedata(:,3) = lifetimedata(:,3) + histc(tcspc(chan==2 & usable_data), 1:65536);
    lifetimedata(:,4) = lifetimedata(:,4) + histc(tcspc(chan==3 & usable_data), 1:65536);
    lifetimedata(:,5) = lifetimedata(:,5) + histc(tcspc(chan==4 & usable_data), 1:65536);

    
    if isfield(tg, 'tg2')
        clear flvp;
        for k=1:4
          for j=1:4
           flvp(:,(((k-1)*4)+j)) = chan==j & tcspc>=tg.tg2(((k-1)*4)+j,1) & tcspc<=tg.tg2(((k-1)*4)+j,2);
          end
        end
        res1.rate(i,:) = sum(flvp)/res1.time(i);   
    end    

end

if isfield(tg, 'tg2')
    res1.bin = lifetimedata(:,1);
    res1.tcspcdata = lifetimedata(:,2:5);
    res1.mcs = mcs;
    for i = 1:4
       for j = 1:4
           if tg.tg2(((4*(j-1))+i),1) <= 0
              tg.tg2(((4*(j-1))+i),1) = 1;
           end
           if tg.tg2(((4*(j-1))+i),2) <= 0
              tg.tg2(((4*(j-1))+i),2) = size(lifetimedata,1);
           end
           res1.tcspc(1:size(lifetimedata(tg.tg2(((4*(j-1))+i),1):tg.tg2(((4*(j-1))+i),2),j),1),j,i) =  lifetimedata(tg.tg2(((4*(j-1))+i),1):tg.tg2(((4*(j-1))+i),2),j);
           res1.tau(1:size(res1.bin(tg.tg2(((4*(j-1))+i),1):tg.tg2(((4*(j-1))+i),2)),1),j,i) = res1.bin(tg.tg2(((4*(j-1))+i),1):tg.tg2(((4*(j-1))+i),2));
       end
    end
 
end    
    
    
    
lifetimedata(:,6) = lifetimedata(:,2) + lifetimedata(:,3) + lifetimedata(:,4) + lifetimedata(:,5);
