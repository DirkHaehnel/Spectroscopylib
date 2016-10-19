function [y, tcspc, flag, markers, num, overcount, head] = tttrRead(name,cnts)

% [y, flag, tcspc, head, num, overcount, marker] = tttrRead(name,cnts);

if ~strcmp(name(end-2:end),'t3r')
    disp('Not a t3r-file!');
    return
else
    fin = fopen(name,'r');
    if (fin==-1)
        errordlg('Cannot open specified file. Please try again.');
    else        
        head.Ident = char(fread(fin, 16, 'char')');    
        head.SoftVersion = char(fread(fin, 6, 'char')');
        head.HardVersion = char(fread(fin, 6, 'char')');
        head.FileTime = char(fread(fin, 18, 'char')');
        fread(fin, 2, 'char');
        head.Comment = char(fread(fin, 256, 'char')');
        head.NChannels = fread(fin, 1, 'int32');
        head.NCurves = fread(fin, 1, 'int32');
        head.BitsPerChannel = fread(fin, 1, 'int32');
        head.RoutingChannels = fread(fin, 1, 'int32');    
        head.NumberOfBoards = fread(fin, 1, 'int32');
        head.ActiveCurve = fread(fin, 1, 'int32');    
        head.MeasMode = fread(fin, 1, 'int32');
        head.SubMode = fread(fin, 1, 'int32');
        head.RangeNo = fread(fin, 1, 'int32');
        head.Offset = fread(fin, 1, 'int32');
        head.TAcq = fread(fin, 1, 'int32');
        head.StopAt = fread(fin, 1, 'int32');
        head.StopOnOvfl = fread(fin, 1, 'int32');
        head.Restart = fread(fin, 1, 'int32');
        head.LinLog = fread(fin, 1, 'int32');
        head.MinAx = fread(fin, 1, 'int32');
        head.MaxAx = fread(fin, 1, 'int32');
        head.MinAxCnt = fread(fin, 1, 'int32');
        head.MaxAxCnt = fread(fin, 1, 'int32');
        
        head.DispCurves = fread(fin, 16, 'int32');
        head.Params = fread(fin, 9, 'int32');
        head.RepeatMode = fread(fin, 1, 'int32');
        head.RepeatsPerCurve = fread(fin, 1, 'int32');
        head.RepeatTime = fread(fin, 1, 'int32');
        head.RepeatWaitTime = fread(fin, 1, 'int32');    
        head.ScriptName = fread(fin, 20, 'char');
        head.BordSerial = fread(fin, 1, 'int32');    
        head.CFD0 = fread(fin, 1, 'int32');
        head.CFDmin = fread(fin, 1, 'int32');
        head.Sync = fread(fin, 1, 'int32');
        head.CurveOffset = fread(fin, 1, 'int32');    
        head.Resolution = fread(fin, 1, 'float');
        head.GlobClock = fread(fin, 1, 'int32');
        
        fread(fin, 6, 'int32');
        
        head.SyncRate = fread(fin, 1, 'int32');
        head.TTTRCDFRate = fread(fin, 1, 'int32');
        head.TTTRStopAfter = fread(fin, 1, 'int32');
        head.TTTRStopReason = fread(fin, 1, 'int32');

        head.Records = fread(fin, 1, 'int32');

        SpecHeaderLength = fread(fin, 1, 'int32');
        if SpecHeaderLength>0
            tmp = fread(fin, 2, 'uint32');
            head.ScanIdent = tmp(2);
            if tmp(2)==1 % PI E710 Scan Controller
                head.ScanTimePerPix = fread(fin, 1, 'int32');
                fread(fin, 1, 'int32');
                head.ScanPattern = fread(fin, 1, 'int32');
                fread(fin, 1, 'int32');
                head.ScanStartX = fread(fin, 1, 'float');
                head.ScanStartY = fread(fin, 1, 'float');
                head.ScanWidthX = fread(fin, 1, 'int32');
                head.ScanWidthY = fread(fin, 1, 'int32');
                head.ScanResolution = fread(fin, 1, 'float');    
                head.ScanTStartTo = fread(fin, 1, 'float');    
                head.ScanTStopTo = fread(fin, 1, 'float');    
                head.ScanTStartFro = fread(fin, 1, 'float');    
                head.ScanTStopFro = fread(fin, 1, 'float');    
            end
            if tmp(2)==2 % SCX 200 Scan Controller
                head.ScanTimePerPix = fread(fin, 1, 'int32');
                head.ScanPause = fread(fin, 1, 'int32');
                head.ScanPattern = fread(fin, 1, 'int32');
                fread(fin, 1, 'int32');
                head.ScanStartX = fread(fin, 1, 'float');
                head.ScanStartY = fread(fin, 1, 'float');
                head.ScanWidthX = fread(fin, 1, 'int32');
                head.ScanWidthY = fread(fin, 1, 'int32');
                fread(fin, 1, 'int32');
            end
        end
        
       if nargin>1
           fseek(fin, 4*(cnts(1)-1), 0);
           [y num] = fread(fin, cnts(2), 'uint32');
           valid = floor(y/2^30);
           flag = floor(y/2^28 - valid*2^2);
           valid = mod(valid,2);
           y = y - floor(y/2^28)*2^28;
           tcspc = floor(y/2^16);
           overflow = (1-valid).*floor(tcspc/2^11);
           markers = (1-valid).*(tcspc - floor(tcspc/2^3)*2^3);

           % overflow = overflow + (marker>0); % what is marker==7???

           y = y - tcspc*2^16;
           tcspc = 2^12 + 1 - tcspc;
           y = y + 2^16*cumsum(overflow);
           valid = overflow==1 | y==0;
           y(valid)=[];
           tcspc(valid)=[];
           flag(valid)=[];
           markers(valid)=[];
           overcount = 2^16*sum(overflow);
           fclose(fin);
       else
           y = head;
       end
    end
end
