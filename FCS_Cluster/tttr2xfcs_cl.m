function [auto, autotime] = tttr2xfcs_cl(y, num1, num2, Ncasc, Nsub, ax)
% Function [auto, autotime] = tttr2xfcs_cl(y, num1, num2, Ncasc, Nsub, ax) calculates
% the crosscorrelation between two photon streams with arrival times y, where
% ind is a logical vector indicating which photon corresponds to
% which photon stream
%
% auto           => autokorrelation (3 dim. matrix)
%                   This matrix has dimensions: 
%                   size(autotime) x size(num1, 2) x size(num2, 2)
% autotime       => time of autokorrelation
%
% Basicly you retrieve the correlation data for the streams selected by
% num1(x) and num2(y) by: auto(:,x,y).
%
% y              => receive time of photons (large timeframe)
% num1           => photons to be considered channel 1 (can be a (boolean) matrix or just a vector)
% num2           => photons to be considered channel 2 (can be a (boolean) matrix or just a vector)
%                    num1 == num2: Autokorrelation
%                    num1 <> num2: Kreuzkorrelation
%
% Ncasc/Nsub     => tribute to hardware correlator times:
%                   limited resolution => es soll eine "semit-log" display
%
%                   (2007-11-23 00:40:13) Bernd Müller: zeit 1 = 1 ns, 2=2ns, 3=3ns...
%                   (2007-11-23 00:40:52) Bernd Müller: 4=6ns, 5=9ns 6=12ns
%                   (2007-11-23 00:41:26) Bernd Müller: 7=24, 8=36, 9=42.... und so weiter
%
%                   You will get Ncasc * Nsub data points
%
% ax             => optional argument, when supplied, it contains a graphics handle
%                   to plot the intermediates to
%
% Assumtions(! maybe better calles educated guesses) based on explanation by Bernd:
% 
% - to derive num1 and num2 from the data of the pt3 file, we have to
%   iterate over the whole y data, check whether this is the right
%   timeframe AND
% - check we select the correct detector channel from the input
%
% If these assumtions are correct, the cross correlation we look for FCS
% are:
%
% Predefinition:
% TG1 = Start of first laser pulse (defined as green - first volume)
% TG2 = End of first laser pulse
% TG3 = Start of second laser pulse (defines as green - second volume)
% TG4 = End of second laser pulse
% TG5-TG8 defined analogous for red laser pulse
%
% Conclusion:
%
% We have to do 12 crosscorrelations:
%
% The first crosscorrelation would be done by specifing num1(1,t) as 1 if:
% - TG1 < y(t) < TG2 and
% - channel == 1
% and num2(1,t) as 1 if:
% - TG1 < y(t) < TG2 and
% - channel == 2
%
% The other correlations should be done by switching channel to 3 and 4 and
% iteration over TG3-8
%
% The interesting question: Can we go over the full range or do we have to
% cut the matrix at the TGs???

dt = double(max(y)-min(y));
% rounding an integer to an integer????
%y = round(y(:));
if size(num1,1)<size(num1,2)
    num1 = num1';
end
if size(num2,1)<size(num2,2)
    num2 = num2';
end
% Postcondition: size(num1/2, 2) < size(num1/2, 1)
autotime = zeros(Ncasc*Nsub,1);
auto = zeros(Ncasc*Nsub,size(num1,2),size(num2,2));
shift = 0;
delta = 1;
for j=1:Ncasc
    [y, k] = unique(y);
	tmp = cumsum(num1);
    num1 = diff([zeros(1,size(num1,2)); tmp(k,:)]);
	tmp = cumsum(num2);
    num2 = diff([zeros(1,size(num2,2)); tmp(k,:)]);
    for k=1:Nsub
        shift = shift + delta;
        lag = round(shift/delta);
        [tmp,i1,i2] = intersect(y+lag,y);
        auto(k+(j-1)*Nsub,:,:) = num1(i1,:)'*num2(i2,:);
        autotime(k+(j-1)*Nsub) = shift;
    end
    y = round(0.5*y);
    num1 = 0.5*num1;
    delta = 2*delta;
end

for j=1:size(auto,1)
    auto(j,:,:) = auto(j,:,:)* (dt./(dt-autotime(j)));
end

