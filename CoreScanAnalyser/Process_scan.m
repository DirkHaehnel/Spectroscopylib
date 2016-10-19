function [head, data, tag, life_imm]=Process_scan(name, task, num_PIE, timegate, flagint, cutoff, pic)
% [head, data, tag, life_imm]=Process_scan(name, task, num_PIE, timegate, flagint,cutoff, pic)
% This program reads scan images (.ht3 files) and can be used to give out
% the intensity image or even the pixel averaged lifetime image.
%
% Input Parameters:
%
% 'name'    : is the name of the .ht3 file, with its complete address.
% 'task'    : is the flag indicating the task the program should perform.
%             task can be 'intensity', to output the intensity image, or 
%             'lifetime', so that the program calculates the lifetime image.
% 'num_PIE' : is the number of laser pulses within one sync cycle. Usually 1
% 'timegate': numbers indicating the starting and ending of each pulse in
% 'cutoff'  : This determines the delay after the excitation pulse after
%             which the counting for average lifetime needs to be done.
%             each detector. For details refer to the help of DetectTimeGates.m
% 'flagint' : If flagint ==1, the intensity image will represent the
%             photons after the cutoff.
% 'pic'     : If pic is 1, the intensity and(or) the lifetime image are
%             displayed.
%
% Output parameters:
% 
% 'head'    : This is the header of the .ht3 file, contains all the
%             information regarding this file.
% 'data'    : Is a struct containing the calculated information. Such as
%             the tcspc histogram, the time axis, the calculated intensity
%             and lifetime images, the positions of the tcspc peak in
%             different detectors and the shifts calculated based on the
%             cutoffs.
% 'tag'     : is the intensity image.
% 'life_imm': is the average pixel photon arrival image (lifetime image).

% (c) Narain Karedla, 2014.


    % name='U:\Narain\140329\Image_012.ht3';
    
    if nargin<7 || isempty(pic)
        pic=0;
    end
    
    if nargin<6 || isempty(cutoff)
        cutoff=0.5; % ns
    end
    if nargin<5 || isempty(flagint)
        flagint=0; % if flagint ==1, The intensity image will represent photons after the cutoff
    else
        flagint=1;
    end
    
    if nargin<4 || isempty(timegate)
        timegate=[];
    end
    
    if nargin<3 || isempty(num_PIE)
        num_PIE=1;
    end
    if nargin<2 || isempty(task)
        task='lifetime';
    end
    if exist([name(1:end-4),'_Core_Scan.mat'],'file')
        load([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
    else
        [head, im_tcspc, im_chan, im_line, im_col] = Core_ScanRead(name);
        save([name(1:end-4),'_Core_Scan.mat'],'head','im_tcspc','im_chan','im_line','im_col');
    end

    dind       = unique(im_chan);
    maxch      = numel(dind);
    maxres     = max([head.Resolution]);
    Resolution = max([maxres 0.032]);
    chDiv      = Resolution/maxres;
    im_tcspc   = ceil(im_tcspc./chDiv);
    Ngate      = double(max(im_tcspc));
    pixel      = head.ImgHdr.PixelSize;
    tcspcdata  = zeros(maxch, Ngate);
    nx         = head.ImgHdr.PixX;
    ny         = head.ImgHdr.PixY;
    t          = (1:Ngate).*Resolution;
    x0         = head.ImgHdr.X0;
    y0         = head.ImgHdr.Y0;
    
    for ch = 1:maxch
        tcspcdata(ch,:) = mHist(double(im_tcspc(im_chan == dind(ch))),1:Ngate);
        indel(ch) = sum(tcspcdata(ch,:))/nx/ny<5; %#ok<*AGROW> % deleting image planes containing less than 25 photons average per pixel
    end
    tcspcdata(indel,:) = [];
    dind(indel)    = [];
    maxch_n        = numel(dind);
    tcspcdata=tcspcdata.';
    
    data.t=t;
    data.tcspcdata=tcspcdata;
    
    if strcmpi(task,'intensity')
       im_pixel=double(im_line)+(double(im_col)-1)*nx + nx*ny*(double(im_chan)-1);
       tag=zeros(ny*nx*maxch_n,1); 
       life_imm=NaN;
       tag=reshape(mHist(im_pixel,1:numel(tag)),[nx,ny,maxch_n]);
       data.tag=tag;
       data.life_imm=life_imm;
              
    elseif strcmpi(task,'lifetime')
        if isempty(timegate)
            [timegate, Ngate] = DetectTimeGates(tcspcdata, num_PIE, Resolution);
        else
            Ngate = 1 + timegate(1,2)+timegate(1,4)-timegate(1,1);
        end
        %                 timegate(  1   ,:): time-window for pulse  1   in detection channel 1
        %                 timegate(  2   ,:): time-window for pulse  1   in detection channel 2
        %            ...  timegate(  n   ,:): time-window for pulse  1   in detection channel n
        %                 timegate( n+1  ,:): time-window for pulse  2   in detection channel 1
        %                 timegate( n+2  ,:): time-window for pulse  2   in detection channel 2
        %            ...  timegate( 2*n  ,:): time-window for pulse  2   in detection channel n
        %            ...  timegate(num_PIE*n,:): time-window for pulse num_PIE in detection channel n
        %
        %                 timegate(:, 1)    : begin of time-window
        %                 timegate(:, 2)    : end of time-window
        %                 timegate(:, 3)    : if > 0: continuation of time-window
        %                 timegate(:, 4)    : if > 0: end of time-window
        
        tcspc=zeros(Ngate,maxch_n,num_PIE); % Ngate timechannels, maxch_n detectors, num_PIE pulses in each.
        ind=timegate(:,3)==0;
        microtime=zeros(numel(im_tcspc),maxch_n,num_PIE);
        im_pulse=zeros(numel(im_tcspc),1);
        for j=1:num_PIE
            for i=1:maxch_n
                if ind((j-1)*maxch_n+i)
                    tcspc(:,i,j)=tcspcdata(timegate((j-1)*maxch_n+i,1):timegate((j-1)*maxch_n+i,2),i);
                    chan_ind=im_chan==i;
                    microtime(chan_ind,i,j)=im_tcspc(chan_ind)-timegate((j-1)*maxch_n+i,1);
                    im_pulse(chan_ind)=j;
                else
                    tcspc(:,i,j)=[tcspcdata(timegate((j-1)*maxch_n+i,1):timegate((j-1)*maxch_n+i,2),i);tcspcdata(timegate((j-1)*maxch_n+i,3):timegate((j-1)*maxch_n+i,4),i)];
                    chan_ind=im_chan==i;
                    timefirst=im_tcspc>timegate((j-1)*maxch_n+i,1)-1;
                    timelast=im_tcspc<timegate((j-1)*maxch_n+i,end)+1;
                    microtime(chan_ind & timefirst,i,j)=im_tcspc(chan_ind & timefirst)-timegate((j-1)*maxch_n+i,1)+1;
                    microtime(chan_ind & timelast,i,j)=im_tcspc(chan_ind & timelast)+timegate((j-1)*maxch_n+i,2)-timegate((j-1)*maxch_n+i,1)+1;
                    im_pulse((chan_ind & timefirst)|(chan_ind & timelast))=j;
                end
            end
        end
        data.timegate=timegate;
        data.tcspc=tcspc;
        
        c=0; pos=zeros(maxch_n*num_PIE,1);
        for i=1:size(tcspc,2)
            for j=1:size(tcspc,3)
                c=c+1;
                [~,pos(c)]=max(tcspc(:,i,j));
            end
        end
        
        shift=ceil(cutoff/Resolution);
        %     tag=zeros(ny,nx,maxch_n,num_PIE); % y pixels, x pixels, detectors and pulses
        %     life_imm=tag;
        
        im_pixel=double(im_line)+(double(im_col)-1)*nx + nx*ny*(double(im_chan)-1)+ nx*ny*maxch_n*(double(im_pulse)-1);
        
        tag=zeros(ny*nx*maxch_n*num_PIE,1); life_imm=ones(numel(tag),1)*10;
        mic_ind=zeros(numel(im_tcspc),1);
        for j=1:num_PIE
            for ch=1:maxch_n
                mic_ind=mic_ind+double(microtime(:,ch,j)>pos((j-1)*maxch_n+ch)+shift+1);
            end
        end
        mic_ind=not(logical(mic_ind));
        if num_PIE>1
            microtime=sum(sum(microtime));
        else
            microtime=sum(microtime,2);
        end
        if maxch~=0
            if flagint
                %         im_tcspc(mic_ind)=[];
                %         im_col(mic_ind)=[];
                %         im_line(mic_ind)=[];
                %         im_chan(mic_ind)=[];
                im_pixel(mic_ind)=[];
                %         im_pulse(mic_ind)=[];
                microtime(mic_ind)=[];
                tag=mHist(im_pixel,1:numel(tag));
                indpix=unique((tag>25).*(1:numel(tag)).');
                if indpix(1)==0; indpix(1)=[]; end
                t0_pixel=zeros(numel(indpix),1);
                for j=1:num_PIE
                    for ch=1:maxch_n
                        ind=(nx*ny*(ch-1)+nx*ny*ch*(j-1))-1<indpix & indpix<(nx*ny*ch+nx*ny*ch*(j-1))+1;
                        t0_pixel(ind)=pos((j-1)*maxch_n+ch)+shift;
                    end
                end
                for i=1:numel(indpix)
                    life_imm(indpix(i))=(mean((microtime(im_pixel==indpix(i))))-t0_pixel(i))*Resolution;
                end
                tag=reshape(tag,[nx,ny,maxch_n,num_PIE]);
                life_imm=reshape(tag,[nx,ny,maxch_n,num_PIE]);
            else
                tag=mHist(im_pixel,1:numel(tag));
                indpix=unique((tag>25).*(1:numel(tag)).');
                if indpix(1)==0; indpix(1)=[]; end
                t0_pixel=zeros(numel(indpix),1);
                for j=1:num_PIE
                    for ch=1:maxch_n
                        ind=(nx*ny*(ch-1)+nx*ny*ch*(j-1))-1<indpix & indpix<(nx*ny*ch+nx*ny*ch*(j-1))+1;
                        t0_pixel(ind)=pos((j-1)*maxch_n+ch)+shift;
                    end
                end
                %         im_tcspc(mic_ind)=[];
                %         im_col(mic_ind)=[];
                %         im_line(mic_ind)=[];
                %         im_chan(mic_ind)=[];
                im_pixel(mic_ind)=[];
                %         im_pulse(mic_ind)=[];
                microtime(mic_ind)=[];
                for i=1:numel(indpix)
                    life_imm(indpix(i))=(mean((microtime(im_pixel==indpix(i))))-t0_pixel(i))*Resolution;
                end
                
                tag=reshape(tag,[nx,ny,maxch_n,num_PIE]);
                life_imm=reshape(life_imm,[nx,ny,maxch_n,num_PIE]);
            end
        end
        
        data.tag=tag;
        data.life_imm=life_imm;
        data.pos=pos;
        data.shift=shift;
        
        
        if pic
            if maxch_n>0
                for j=1:num_PIE
                    close all
                    figure;
                    set (gcf,'name','Intensity Image','NumberTitle','off')
                    for ch=1:maxch_n
                        subplot(1,maxch_n,ch)
                        imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), data.tag(:,:,ch,j))
                        colorbar
                        axis equal
                        axis tight
                        colormap('hot')
                    end
                end
                for j=1:num_PIE
                    close all
                    figure;
                    set (gcf,'name','Lifetime Image','NumberTitle','off')
                    for ch=1:maxch_n
                        subplot(1,maxch_n,ch)
                        imagesc(x0+pixel.*(1:nx),y0+pixel.*(1:ny), data.life_imm(:,:,ch,j)) % check this!!
                        colorbar
                        axis equal
                        axis tight
                        colormap('hot')
                    end
                end
            end
        end
    end
    
end 
    
