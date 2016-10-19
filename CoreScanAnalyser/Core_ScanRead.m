function [head, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name)
% Core_ScanRead reads the .ht3 file of an image scan
% im_tcspc is a vector containing the time channel of each photon
% im_chan is a vector where the values represent the detector number of
% each photon
% im_line is a vector containing the line number of all the photons 
% im_col is a vector containing the column number of all the photons 

    if strcmp(name(end-2:end),'ht3')

        head = HT3_Read(name);

        if ~isempty(head)

            nx         = head.ImgHdr.PixX;                                              % x-pixels
            ny         = head.ImgHdr.PixY;                                              % y-pixels
%             dt         = 1e-9*nx*head.ImgHdr.PixelTime*head.SyncRate;                   % number of pulses in one line
            xcord      = head.ImgHdr.X0+(1:nx)*head.ImgHdr.PixelSize;
            ycord      = head.ImgHdr.Y0+(1:ny)*head.ImgHdr.PixelSize;
            
            head.ImgHdr.Xcord = xcord;
            head.ImgHdr.Ycord = ycord;
            
            anzch      = 32;                                                    % max number of detectors
            Resolution = max([head.Resolution]);                                % Resolution of measurement
            chDiv      = Resolution/head.Resolution;                            % 1 in this case
            Ngate      = ceil(1e9/head.SyncRate./Resolution)+1;                 % TCPSC channels

            y        = [];
            tmpx     = [];
            chan     = [];
            markers  = [];

            im_tcspc = [];
            im_chan  = [];
            im_line  = [];
            im_col   = [];
            Turns    = [];

            cnt      = 0;
            tend     = 0; 
            line     = 1;

            h.waitbar = waitbar(0,'Please wait ...');

            fprintf('\n\n');

            if head.ImgHdr.Pattern == 0 % Unidirectional scan

                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);
                                                                                % Reading 5 million counts at a time
                while (num>0)

                    cnt = cnt + num;                
                    if ~isempty(y)
                        tmpy = tmpy+tend;                    
                    end;

                    ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));

                    y       = [y; tmpy(ind)];                         %#ok<AGROW>
                    tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                    chan    = [chan; tmpchan(ind)+1];                 %#ok<AGROW>
                    markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>

                    Turns   = [Turns; y(markers==1 & (chan==2|chan==3|chan==5))]; %#ok<AGROW>     % Sync for the change of line 

                    tend  = y(end)+loc;

                    if numel(Turns)>1
                        for j=1:numel(Turns)-1

                            dT = (Turns(2)-Turns(1));                           % Syncs in each line
                            t1 = Turns(1)+head.ImgHdr.TStartTo*dT;              % sync number for the begin of a to- scan
                            t2 = Turns(1)+head.ImgHdr.TStopTo*dT;               % sync number for the end of a to- scan

                            ind = (markers~=0)|(y<t1);
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];

                            ind = (y>=t1)&(y<=t2);

                            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW> % sync number of photons 
                            im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW> % detector number of photons
                            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW> % line number of photons
                            im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW> % separation of photons into pixels based on sync numbers

                            line = line +1;
                            waitbar(line/ny);
                            drawnow

                            Turns(1) = [];

                        end                    
                    end

                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);

                end;

                t1 = Turns(end)+head.ImgHdr.TStartTo*dT;
                t2 = Turns(end)+head.ImgHdr.TStopTo*dT;

                ind          = (y<t1);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];

                ind = (y>=t1)&(y<=t2);

                im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              
                im_chan   = [im_chan; uint8(chan(ind))];                
                im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  
                im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  

                line = line +1;
                h.waitbar = waitbar(line/ny);
                drawnow

            else  % bidirectional scan

                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);
                % tmpy            = sync value of each photon detected
                % tmptcspc        = tcspc channel of each photon detected
                % tmpchan         = detector channel for each detected photon
                % tmpmarkers      = marker indicating a virtual photon or
                %                   event. In this case, when the piezo shifts to the next line it is 1
                % num             = number of photons read actually.

                % tmpchan also has a different value when there is a line shift
                while (num>0)

                    cnt = cnt + num;
                    if ~isempty(y)
                        tmpy = tmpy+tend;                    
                    end;

                    ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<=Ngate*chDiv)); % indicator for line shift

                    y       = [y; tmpy(ind)];                         %#ok<AGROW> 
                    tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                    chan    = [chan; tmpchan(ind)+1];                 %#ok<AGROW>
                    markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>

                    Turns = [Turns; y(markers==1 & (chan==2|chan==3|chan==5))]; %#ok<AGROW>

                    tend = y(end)+loc;

                    if numel(Turns)>2
                        for j=1:2:2*floor(numel(Turns)/2-1)

                            dT = (Turns(2)-Turns(1));                           % Syncs in each line
                            t1 = Turns(1)+head.ImgHdr.TStartTo*dT;              % sync number for the begin of a to- scan
                            t2 = Turns(1)+head.ImgHdr.TStopTo*dT;               % sync number for the end of a to- scan

                            ind = (markers~=0)|(y<t1);
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];

                            ind = (y>=t1)&(y<=t2);

                            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW> % sync number of photons 
                            im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW> % detector number of photons
                            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW> % line number of photons
                            im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW> % separation of photons into pixels based on sync numbers


                            line = line +1;

                            t1 = Turns(1)+head.ImgHdr.TStartFro*dT;             % begin of fro- direction scan
                            t2 = Turns(1)+head.ImgHdr.TStopFro*dT;              % end of the fro- direction scan

                            ind = (y<t1);
                            y(ind)       = [];
                            tmpx(ind)    = [];
                            chan(ind)    = [];
                            markers(ind) = [];

                            ind = (y>=t1)&(y<=t2);

                            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                            im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                            im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>

                            line = line +1;
                            waitbar(line/ny);
                            drawnow

                            Turns(1:2) = [];
                        end
                    end

                    [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);

                end;

                if ~isempty(Turns) % for the last one or two lines
                    t1 = Turns(end-1)+head.ImgHdr.TStartTo*dT;
                    t2 = Turns(end-1)+head.ImgHdr.TStopTo*dT;

                    ind = (y<t1);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];

                    ind = (y>=t1)&(y<=t2);

                    im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                    im_chan   = [im_chan; uint8(chan(ind))];
                    im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                    im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];

                    line = line +1;

                    t1 = Turns(end-1)+head.ImgHdr.TStartFro*dT;
                    t2 = Turns(end-1)+head.ImgHdr.TStopFro*dT;

                    ind = (y<t1);
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];

                    ind = (y>=t1)&(y<=t2);

                    im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                    im_chan   = [im_chan; uint8(chan(ind))];
                    im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                    im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];

                    line = line +1;
                    h.waitbar = waitbar(line/ny);
                    drawnow
                end
            end;

            delete(h.waitbar);

        end
    end
 
end
