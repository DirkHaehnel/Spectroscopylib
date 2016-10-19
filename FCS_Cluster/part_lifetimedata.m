function [lifetimedata, head] = part_lifetimedata(ff,package_size,tg1, statusfield)
% [lifetimedata, head] = part_lifetimedata(ff,package_size,tg1, statusfield
%
% ff           => file to read
% package_size => # of photons to read per cycle
% tg1          => MCS-Timegates (globale timegate) [s]
% statusfield  => optional GUIDE textfield to print statusdata to
tg1

head = ht3v2read_cl(ff, [0 0]);
Timeunit=1/head.SyncRate;

%Timeunit = head.Sync * 1e-9;
use_tg1 = tg1 ./ Timeunit;

count_packages = ceil(head.Records / package_size);
old_overcount = 0;

lifetimedata = zeros(4096,6);
lifetimedata(:,1) = (1:4096) * head.Resolution;

for i=1:count_packages
    if exist('statusfield', 'var')
        set(statusfield, 'String', sprintf('Lifetime data: %03i / %03i', i, count_packages));
        drawnow;
    end
    fprintf('i: %i - max: %i\n', i, count_packages);
    si = (i-1)*package_size+1;
    ei = i*package_size;
    if ei > head.Records
        ei = head.Records;
    end
%    [head, z, tcspc, chan, markers, num, overcount] = pt3v2read(ff, [si ei-si+1]);
    [z, tcspc, chan, markers, num, head] = ht3v2read_cl(ff, [si ei-si+1]);

%    z = (z+old_overcount);
%    old_overcount = old_overcount + overcount;
    usable_data = z >= use_tg1(1) & z < use_tg1(2);
    lifetimedata(:,2) = lifetimedata(:,2) + histc(tcspc(chan==1 & usable_data), 1:4096);
    lifetimedata(:,3) = lifetimedata(:,3) + histc(tcspc(chan==2 & usable_data), 1:4096);
    lifetimedata(:,4) = lifetimedata(:,4) + histc(tcspc(chan==3 & usable_data), 1:4096);
    lifetimedata(:,5) = lifetimedata(:,5) + histc(tcspc(chan==4 & usable_data), 1:4096);
end

lifetimedata(:,6) = lifetimedata(:,2) + lifetimedata(:,3) + lifetimedata(:,4) + lifetimedata(:,5);
