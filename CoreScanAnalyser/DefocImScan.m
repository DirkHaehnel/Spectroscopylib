function [data, imstack] = DefocImScan(dirname,campixel, pic, printim)
% [data, imstack] = DefocImScan(dirname,campixel, pic, printim)
% This program is for analyzing defocused images taken at points located by
% the search molecule LabView software.The scan images which have defocused
% images before the next scan image are read out and the intensity plot is
% shown in the end together with the defocused images taken on this
% image. The data is saved in a superstructure 'data' and each substructure
% inside corresponds to each scan image & defocused image family.'imstack'
% is the stack of all the display images containing the defocused images.
% 
% Input parameters: 
%
% 'dirname' is the name of the directory. An example is given below
% 'campixel' is the effective pixel dimension (actual pixel size/Magnification)
% 'pic' is a flag to display the intensity images along with the defocused
% images.
% 'printim' is a flag to save intensity images and the defocused images in
% the same directory.
% 
% Function Dependencies:
%
% Process_scan.m, Core_ScanRead.m, HT3_Read.m, DetectTimeGates.m, mHist.m,
% mCluster.m and fuzzylocate.m

% (c) Narain Karedla, 2014

% dirname='U:\Narain\140513\';
if nargin<4 || isempty(printim)
    printim=0;
end

if nargin<3 || isempty(pic)
    pic=0;
end

if nargin<2 ||isempty(campixel)
    campixel=0.06;
end
dispimsize=1e3; % size of the display image
ht3names=dir([dirname '*.ht3']);
idxht3= false(numel(ht3names),1);
% idxauto = idxht3;
for i=1:numel(ht3names)
    idxht3(i)=~isempty(strfind(ht3names(i).name,'Image'));
    % dateht3(i)=ht3names(i).datenum;
end
%
% for i=1:numel(ht3names)
%     idxauto(i)=~isempty(strfind(ht3names(i).name,'AutoPoint'));
% end
defocnames = dir([dirname '*.mat']);
ind=zeros(numel(defocnames,1));
for i=1:numel(defocnames)
    if ~isempty(strfind(defocnames(i).name,'Core')) % Although this doesn't matter much, deleting the list of Core_Scan.mat files
        ind(i)=1;
    end
end
defocnames(logical(ind))=[];
datedefoc=zeros(numel(defocnames),1);
for i=1:numel(defocnames)
    datedefoc(i)=defocnames(i).datenum;
end
idxht3=find(idxht3);
dateht3=zeros(numel(idxht3),1);
for j=1:numel(idxht3)
    dateht3(j)=ht3names(idxht3(j)).datenum;
end
% dateauto=zeros(sum(idxauto),1);
% idxauto=find(idxauto);
% for j=1:numel(idxauto)
%     dateauto(j)=ht3names(idxauto(j)).datenum;
% end
[dateht3,idx]=sort(dateht3);
Imagenames=ht3names(idxht3(idx));% sorted according to the real time.
% [dateauto,idx]=sort(dateauto);
% Autonames=ht3names(idxauto(idx));

dateht3=[dateht3;Inf];

for i=2:numel(Imagenames)+1
    timind=datedefoc<dateht3(i);
    %     timindauto=dateauto<dateht3(i);
    %     tmpAutonames=Autonames(timindauto);
    tmpImnames=defocnames(timind);
    tmpdates=datedefoc(timind);
    [~,ordert]=sort(tmpdates);
    tmpImnames=tmpImnames(ordert);
    if sum(timind)>0 %&& sum(timindauto)>0
        [head,~,tag,~]=Process_scan([dirname Imagenames(i-1).name],'intensity');
        tag=sum(tag,3); % Intensity image
        x=zeros(numel(tmpImnames),1); y=x; defocim=cell(numel(tmpImnames),1);
        %         for k=1:numel(timindauto);
        %             head=HT3_Read([dirname tmpAutonames(k).name]);
        %         end
        
        for j=1:numel(tmpImnames)
            load([dirname tmpImnames(j).name],'results')
            pixel=sqrt(results.Info.CcdSize);
            if pixel-round(pixel)==0 % Assuming that the CCD image is always a square number
                imm=reshape(results.ImagePlane(end:-1:1),pixel,pixel)'; % depending on the direction of the camera
                x(j)=results.Info.PositionX;
                y(j)=results.Info.PositionY;
                bim=imm-min(imm(:))>0.25*(max(imm(:)));
                if sum(bim(:))>0.9*(pixel^2)
                    defocim{j}=[];
                    
                else
                    [bim,~]=mCluster(bim,20);  % defocused image should be atleast 20 pixels
                    xim=bim.*reshape(1:numel(imm),pixel,pixel); % for estimating x coordinates of the image
                    xmax=ceil(max(xim(:))/pixel);
                    xmin=floor(min(xim(xim>0))/pixel);
                    yim=bim.*reshape(1:numel(imm),pixel,pixel).'; % for estimating y coordinates of the image
                    ymax=ceil(max(yim(:))/pixel);
                    ymin=floor(min(yim(xim>0))/pixel);
                    window=max((xmax-xmin),(ymax-ymin));
                    center=[round((xmax-xmin)/2)+xmin,round((ymax-ymin)/2)+ymin];
                    if window/2-round(window/2)==0
                        window=window+2;
                    else
                        window=window+3;
                    end
                    
                    if ~isempty(window) && ~isempty(center) && window/2+center(1)<pixel && window/2+center(2)<pixel
                        defocim{j}=imm(max(center(2)-window/2,1):center(2)+window/2,max(center(1)-window/2,1):center(1)+window/2);
                        defocim{j}=defocim{j}-min(defocim{j}(:));
                    end
                end
            end
        end
        data(i).head=head;
        data(i).tag=tag;
        data(i).ImageInfo=Imagenames(i-1);
        data(i).pixel=pixel;
        data(i).x=x; % x values sorted in time
        data(i).y=y; % y values sorted in time
        data(i).defocim=defocim;
        data(i).TotalImages=sum(timind);
    end
    datedefoc(timind)=[];
    defocnames(timind)=[];
    %  dateauto(timindauto)=[];
    %  Autonames(timindauto)=[];
end

datind=false(numel(data),1); % I dont't know a smarter way of doing this.
for i=1:numel(data)
    if isempty(data(i).TotalImages)
        datind(i)=true;
    end
end
data(datind)=[];

% Here begins image stacking
imstack=cell(numel(data),1);
for i=1:numel(data)
    x0=data(i).head.ImgHdr.X0;
    y0=data(i).head.ImgHdr.Y0;
    nx=data(i).head.ImgHdr.PixX;
    ny=data(i).head.ImgHdr.PixY;
    scanpixel=data(i).head.ImgHdr.PixelSize;
    xcoords=x0:nx*scanpixel/dispimsize:(x0+nx*scanpixel);
    ycoords=y0:ny*scanpixel/dispimsize:(y0+ny*scanpixel);
    xc=data(i).x; % sorted in time
    yc=data(i).y; % sorted in time
    %     layers = (repmat(1:sum(xc==xc(1) & yc==yc(1)),[1,numel(xc)/sum(xc==xc(1) & yc==yc(1))])).';
    nlay = sum(xc==xc(1) & yc==yc(1)); % The first point coordinates will be the first recorded again.
    npos=numel(unique(xc.*yc));
    if numel(xc)<nlay*npos % Incomplete layers
        layers = (repmat(1:nlay-1,[npos,1]));
        layers = layers(:);
        layers = [layers; ones(numel(xc)-npos*(nlay-1),1)*nlay];
    else
        layers = (repmat(1:nlay,[floor(numel(xc)/nlay),1]));
        layers = layers(:); % since the images are sorted by the time of their acquisition, all images concerning to one layer should occur before the next layer.
    end
    dispim=zeros(dispimsize,dispimsize,numel(unique(layers)));
    [~,idx] = fuzzylocate([xc,yc],[xcoords; ycoords]);
    c=0;
    for j=1:numel(unique(layers))
        ind=find(layers==j);
        for k=1:numel(ind)
            window=size(data(i).defocim{ind(k)});
            tmplocx=idx(ind(k),1)-round(window(1)/2)+1:idx(ind(k),1)+round((window(1)-1)/2);
            tmplocy= idx(ind(k),2)-round(window(2)/2)+1:idx(ind(k),2)+round((window(2)-1)/2);
            shift= [1-min([tmplocx(tmplocx<=0).';0]),1-min([tmplocy(tmplocy<=0).';0])];
            dispim(tmplocx+shift(1),tmplocy+shift(2),j)=data(i).defocim{ind(k)};
            if ~isempty(data(i).defocim{ind(k)})
                c=c+1;
            end
        end
    end
    imstack{i}=dispim;
    data(i).imtot=c;                                                          %#ok<*AGROW>
    data(i).xcoords=xcoords;
    data(i).ycoords=ycoords;
end

if pic
    close all
    for i=1:numel(data)
        if data(i).imtot>0
            layers=size(imstack{i},3);
            x0=data(i).head.ImgHdr.X0;
            y0=data(i).head.ImgHdr.Y0;
            nx=data(i).head.ImgHdr.PixX;
            ny=data(i).head.ImgHdr.PixY;
            scanpixel=data(i).head.ImgHdr.PixelSize;
            figure
            imagesc(y0:y0+ny*scanpixel,x0:x0+nx*scanpixel,data(i).tag);
            title(['Scan Image ',num2str(i)])
            axis equal
            axis tight
            colormap('hot')
            colorbar
            xlabel('\mum')
            ylabel('\mum')
            for j=1:layers
                if sum(sum(imstack{i}(:,:,j)))>0;
                    figure
                    imagesc(data(i).ycoords,data(i).xcoords,imstack{i}(:,:,j))
                    title(['Defocused Images for Scan Image ',num2str(i),' layer ',num2str(j)])
                    axis equal
                    axis tight
                    colormap('gray')
                    colorbar
                    xlabel('\mum')
                    ylabel('\mum')
                    
                end
            end
        end
    end
end

if printim
    close all
    for i=1:numel(data)
        x0=data(i).head.ImgHdr.X0;
        y0=data(i).head.ImgHdr.Y0;
        nx=data(i).head.ImgHdr.PixX;
        ny=data(i).head.ImgHdr.PixY;
        scanpixel=data(i).head.ImgHdr.PixelSize;
        figure
        imagesc(y0:y0+ny*scanpixel,x0:x0+nx*scanpixel,data(i).tag);
        title(['Scan Image ',num2str(i)])
        axis equal
        axis tight
        colormap('hot')
        colorbar
        xlabel('\mum')
        ylabel('\mum')
        print('-dpng','-r300',[dirname,'Intensity ',data(i).ImageInfo.name(1:end-4)])
        c=0;
        for j=1:numel(data(i).defocim)
            if ~isempty(data(i).defocim{j})
                clf
                c=c+1;
                imagesc((1:size(data(i).defocim{j},1))*campixel,(1:size(data(i).defocim{j},2))*campixel,data(i).defocim{j})
                axis equal
                axis tight
                colormap('gray')
                colorbar
                xlabel('\mum')
                ylabel('\mum')
                print('-dpng','-r300',[dirname,'defocusedimage ',data(i).ImageInfo.name(1:end-4),'_',num2str(c)])
            end
        end
    end
end
end


