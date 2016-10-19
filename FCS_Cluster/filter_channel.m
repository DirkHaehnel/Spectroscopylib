function selected = filter_channel(y, tcspc, chan, tg1, tg2, channel)
% [selected] = filter_channel(y, tcspc, chan, tg1, tg2, channel)
%
% Function to calculate num1/num2 vectors for tttr2xfcs_cl
%
% y         => receive time of photons (large timeframe)
%              measured as some cylce count
% tcspc     => receive time of photons (small timeframe) 
%              measured as 0-4096 slots
% chan      => detector channel
% tg1       => timegate in large timeframe ( count x 2 matrix )
% tg2       => timegate in small timeframe ( count x 2 matrix )
% channel   => detector channel to select (count x 1 matrix)
%
% selected  => given the constraints of tg1, tg2 and channel this
%              array holds the information whether this datapoint in y
%              is to be considered

% For each set of constraints in tg1, tg2, channel we generate one vector
% which is stored a column in the selected matrix
%
% So: size(tg1, 1) == size(tg2, 1) == size(channel, 1) AND
%     size(tg1, 2) == size(tg2, 2) == 2  AND
%     size(channel, 2 == 1
for i=1:size(tg2,1)
selected(:, i) = (y > tg1(i,1)) & (y < tg1(i,2)) & ...
                 (tcspc > tg2(i, 1)) & (tcspc < tg2(i, 2)) & ...
                 (chan == channel(i));
end        
