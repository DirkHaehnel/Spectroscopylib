function [head, tag, tcspc, tcspc_x, tcspc_y, DC, AP] = MT_ScanRead(name, num_PIE, bg, timegate)

if (nargin<4)||isempty(timegate)
    timegate = [];
end;

if (nargin<3)||isempty(bg)
    BG_AP = 0;
else
    BG_AP = bg;
end;

if (nargin<2)||isempty(num_PIE)
    num_PIE = 1;
end;

fname =  [name(1:end-4) '_FLIM.dat'];
tname =  [name(1:end-4) '_INT.dat'];


if strcmp(name(end-2:end),'ht3')

    [head, im_tcspc, im_chan, im_line, im_col] = HT3_ScanRead(name);
    
elseif strcmp(name(end-2:end),'pt3')

    [head, im_tcspc, im_chan, im_line, im_col] = PT3_ScanRead(name);
    head.ImgHdr.PixX = head.ImgHdr.ScanWidthX;
    head.ImgHdr.PixY = head.ImgHdr.ScanWidthY;

end


if ~isempty(head)

    nx = head.ImgHdr.PixX;
    ny = head.ImgHdr.PixY;

    maxch = max(im_chan);
    
    maxres     = max([head.Resolution]);
    Resolution = max([maxres 0.032]);
    chDiv      = Resolution/maxres;
    
    im_tcspc = ceil(im_tcspc./chDiv);
    Ngate    = double(max(im_tcspc));

    tcspc   = zeros(maxch, Ngate);    
    tcspc_x = zeros(nx, maxch, Ngate);
    tcspc_y = zeros(ny, maxch, Ngate);
        
    for ch = 1:maxch
        tcspc(ch,:) = mHist(double(im_tcspc(im_chan == ch)),1:Ngate);
    end

    for y = 1:ny
        ind = (im_line == y);
        tmp1 = im_tcspc(ind);
        tmp2 = im_chan(ind);
        for ch = 1:maxch
            tcspc_y(y,ch,:) = mHist(double(tmp1(tmp2 == ch)),1:Ngate);
        end
    end

    for x = 1:nx
        ind = (im_col == x);
        tmp1 = im_tcspc(ind);
        tmp2 = im_chan(ind);
        for ch = 1:maxch
            tcspc_x(x,ch,:) = mHist(double(tmp1(tmp2 == ch)),1:Ngate);
        end
    end
                    
    nch                = 1:maxch;
    chind              = sum(tcspc,2)<=100;
    nch(chind)         = [];
    tcspc(chind,:)     = [];
    tcspc_x(:,chind,:) = [];
    tcspc_y(:,chind,:) = [];
    
%     ind              = sum(tcspc,1)<=10;
%     tcspc(:,ind)     = [];
%     tcspc_x(:,:,ind) = [];
%     tcspc_y(:,:,ind) = [];
    
    tau   = ((1:size(tcspc,2))-0.5).*Resolution;
    nbin  = numel(tau);    
    anzch = numel(nch);
    
    % Estimate background in measurement
    
    c_t   = [tcspc_x./ny; tcspc_y./nx];
    c_ts  = sum(c_t,3);
    nc    = size(c_ts,1);
    ind   = (1:floor(nc/100):nc);
    AP    = zeros(anzch,1);
    DC    = zeros(anzch,1);
    
    for a = 1:anzch
        [tmp, ord] = sort(c_ts(:,a),'descend');
        cnt = tmp(ind);
        tmp = squeeze(c_t(ord(ind),a,:));
        tmp = sort(tmp,2,'ascend');
        b_g = mean(tmp(:,2:floor(0.25.*nbin)),2);
        
        if BG_AP == 0
            V = [ones(size(cnt)) cnt];
            p  = lsqnonneg(V, b_g);
            DC(a) = p(1);
            AP(a) = p(2);
        else
            AP(a) = BG_AP;
            DC(a) = mean(b_g-cnt.*BG_AP);
        end
    end
    DC(DC<0) = 0;
    AP(AP<0) = 0;
    
    % Determine time-gates for PIE
    
    if isempty(timegate)
        [timegate, Ngate] = DetectTimeGates(tcspc', num_PIE, Resolution);
    else
        Ngate = 1 + timegate(1,2)+timegate(1,4)-timegate(1,1);
    end
            
    head.timegate = timegate;
    
    tcspc   = zeros(anzch, Ngate, num_PIE);
    tcspc_x = zeros(nx, anzch, Ngate, num_PIE);
    tcspc_y = zeros(ny, anzch, Ngate, num_PIE);
    tag     = zeros(ny, nx,anzch, num_PIE);
        
    fid = fopen(fname,'w');
    tid = fopen(tname,'w');
    
    for a = 1:ny
        
        tmp = zeros(nx, anzch, Ngate, num_PIE);

        ind = (im_line == a);
        tmp1 = im_tcspc(ind);
        tmp2 = im_chan(ind);
        tmp3 = im_col(ind);
        
        for x = 1:nx          
            ind = (tmp3 == x);
            tmp4 = double(tmp1(ind));
            tmp5 = tmp2(ind);
            
            n = 0;
            for pulse = 1:num_PIE
                for ch = 1:anzch
                    h2 = [];
                    n  = n+1;
                    h1 = mHist(tmp4(tmp5==nch(ch)),timegate(n,1):timegate(n,2));
                    if isempty(h1)
                        h1 = zeros(1+timegate(n,2)-timegate(n,1),1);
                    end
                    if (timegate(n,3)~=0)                    
                         h2 = mHist(tmp4(tmp5==nch(ch)),timegate(n,3):timegate(n,4));
                         if isempty(h2)
                             h2 = zeros(1+timegate(n,4)-timegate(n,3),1);
                         end
                    end
                    tmp(x,ch,:,pulse) = [h1; h2];
                end
            end
        end
                
        tcspc   = tcspc + shiftdim(sum(tmp,1),1);
        tcspc_x = tcspc_x + tmp;
        tcspc_y(a,:,:,:) = shiftdim(sum(tmp,1),1);
        tag(a,:,:,:) = permute(sum(tmp,3),[1 2 4 3]);
        
        fwrite(fid, uint16(tmp), 'uint16' );
        fwrite(tid, tag(a,:,:,:), 'double');
    end
    
    fclose(tid);
    fclose(fid);
               
end