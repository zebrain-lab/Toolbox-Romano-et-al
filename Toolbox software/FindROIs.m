%%
%%%%%%%%%%%%%%%% STEP 1: Read video %%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;clc;
%%%%%%%%%%%%%%% PARAMETERS 1 %%%%%%%%%%%%%%%%%

[filename,pathname] = uigetfile({'*.tif';'*.TIFF'},'Open TIFF file with imaging video', 'MultiSelect', 'off');
fileName=fullfile(pathname,filename);

prompt = {'\mum per pixel in x direction','\mum per pixel in y direction'};
dlg_title = 'Input data';
num_lines = 1;
def = {'1','1'};
opts.Interpreter='tex';
answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
pixelLengthX= str2num(answer{1});
pixelLengthY= str2num(answer{2});

ansMethod = questdlg('Select method to define ROIs','ROI definition method', 'Automatically detect single-neuron ROIs', 'Use hexagonal grid of ROIs','Import or manually draw all ROIs','Automatically detect single-neuron ROIs');

if strcmp(ansMethod,'Use hexagonal grid of ROIs')
    prompt = {'Diameter of hexagonal ROIs (in \mum)'};
    dlg_title = 'Size of ROIs';
    num_lines = 1;
    def = {'7'};
    opts.Interpreter='tex';
    answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
    sizeROIs= str2num(answer{1});

end

if strcmp(ansMethod,'Import or manually draw all ROIs')
    ansMethod = questdlg('Select method to define ROIs','ROI definition method', 'Manually draw all ROIs','Import ROIs', 'Manually draw all ROIs');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v = ver;
%parOn=any(strcmp('Parallel Computing Toolbox', {v.Name}));
parOn=0;
%parOn=1; %  UNCOMMENT THIS FOR PARALLELIZATION
if parOn
    numCores=feature('numcores');
    
    
    isOpen = matlabpool('size') > 0;
    if ~isOpen
        matlabpool('open','local', round(numCores*0.8))
    end
    
end

% open video
tiffInfo = imfinfo(fileName);
numFrames = numel(tiffInfo);
video=zeros(numFrames,tiffInfo(1).Height,tiffInfo(1).Width, 'single');
disp('Reading movie file...')

if parOn
    parfor i=1:numFrames
        video(i,:,:)=im2single(imread(fileName,i));
    end
else
    for i=1:numFrames
        
        video(i,:,:)=im2single(imread(fileName,i));
    end
end



%%
% graphic propieties
text_si=18;
scrsz = get(0, 'ScreenSize');
set(0,'DefaultAxesFontSize',text_si,'DefaultFigureColor','w', 'DefaultAxesTickDir', 'in','DefaultFigureWindowStyle','normal',...
    'DefaultFigurePosition', [1 1 scrsz(3) scrsz(4)])
if strcmp(ansMethod,'Import ROIs')
    step='ImportROIs';
else
    step='ROI_Mask';
end
end_pro=0;
clear first_autoROI first_okROI


while end_pro==0
    switch step
        case 'ImportROIs'
            [file,pth] = uigetfile({'*.mat';'*.MAT'},'Select file with ROIs to import', 'MultiSelect', 'off');
            fileROIs=fullfile(pth,file);
            load(fileROIs);
            
            CC2 = bwconncomp(importedROIs,4);
            cell_per=cell(CC2.NumObjects,1);
            cells=cell(CC2.NumObjects,1);
            cell_number=CC2.NumObjects;
            f=squeeze(mean(video,1));
            f=f/max(max(f));
            avg=f;
            
            bkg=importedROIs;
            
            figure
            imagesc(f);
            shading flat; colormap gray; set(gca,'xticklabel','','yticklabel','');axis image;
          
            for i=1:CC2.NumObjects
                imgROI=logical(zeros(size(importedROIs)));
                imgROI(CC2.PixelIdxList{i})=1;
                pointsTemp=bwboundaries(imgROI);
                cell_per{i}=fliplr(pointsTemp{1});
                cells{i}=CC2.PixelIdxList{i};
                lp=line(cell_per{i}(:,1),cell_per{i}(:,2)); set(lp,'color',[1 0 0],'LineWidth',1);
         
            end
          
            
               
            
            step='DF_save';

        case 'ROI_Mask'
            f=squeeze(mean(video,1));
            fOld=f;
            f=f/max(max(f));
            avg=f;
            h=figure('name','Draw a mask that contains the region with ROIs');
            imagesc(f);
            shading flat; colormap gray; set(gca,'xticklabel','','yticklabel','');axis image;
            drawnow
            bw = roipoly;
            mask=bw;
            close(h)
            
            if strcmp(ansMethod,'Automatically detect single-neuron ROIs')
                %f(~mask)=0;
                step='Adjust_Gamma';
            elseif strcmp(ansMethod,'Manually draw all ROIs')
                f(~mask)=0;
                CC2=[];
                temp=zeros(size(f));
                points2=[];
                avg=f;
                step='Adjust_Gamma';
            elseif strcmp(ansMethod,'Use hexagonal grid of ROIs')
                avg=f;
                [cells,bkg]=HexaSegment(f,mask,sizeROIs);
                cell_number=length(cells);
                
                cell_per=cell(cell_number,1);
                
                for i=1:cell_number
                    imgROI=zeros(size(avg));
                    imgROI(cells{i})=1;
                    tempBoundary=bwboundaries(imgROI,4);
                    cell_per{i}=fliplr(tempBoundary{1});
               
                end
                
                figure; 
                imagesc(f);
                shading flat; colormap gray; set(gca,'xticklabel','','yticklabel','');axis image;
                hold on
                for i=1:cell_number
                    lp=line(cell_per{i}(:,1),cell_per{i}(:,2));
                    set(lp,'color',[1 0 0],'LineWidth',1);
                end
                
                step='DF_save';
                            
            end
            
        case 'Adjust_Gamma'
            %%%%%%%%%%%%%%% STEP 3: Adjust Gamma %%%%%%%%%%%%%%%%%%%%%%%%
            valGam=1;
            valContr=1;
            locContrastOK=0;
            myColors=[linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr));linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr));linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr))]';
            gammaOK=0;
            myColorsOrig=myColors;
            
            hf=figure(200); clf; hold on
            set(gcf,'Name','Adjust Gamma & Local Contrast', 'Toolbar','figure');
            subplot(1,2,1)
            imagesc(avg); axis image; axis ij; hold on
            set(gca,'xticklabel','','yticklabel','')
            posIm=get(gca,'Position');
            width=.1; height=.1;
            x=posIm(1); y=min(1-height-.01,posIm(2)+posIm(4)+.03);
            uicontrol('Style','pushbutton','String','Create masks again','CallBack','step=''ROI_Mask'';gammaOK=1;locContrastOK=1;','Units','normalized','Position',[x y width height]);
            
            x=max(.01,posIm(1)-width-.01); y=posIm(2)+.1; 
            uicontrol('Style','pushbutton','String','Continue','CallBack','gammaOK=1;step=''Regional_minima'';','Units','normalized','position',[x y width height]);
            
            width=.1; height=.05;
            x=posIm(1); y=max(posIm(2)-2*height-0.01,0.01);  yGamma=y; 
            hGam1=uicontrol('Style','slider','String','Set Gamma','Units','normalized','Position',[x y width height],'Min', 0.005,'Max',2,...
                'value', valGam,'SliderStep',[.01 .05],'Callback','[myColors valGam]=changeGamma(hGam1, valContr,myColors);');
            
            x=x+width+0.01; y=max(posIm(2)-2*height-0.01,0.01); yCont=y;
            hGam2=uicontrol('Style','slider','String','Set Contrast','Units','normalized','position',[x y width height],'Min',1,'Max',300,...
                'value', valContr,'SliderStep',[.1 .5],'Callback','[myColors valContr]=changeContrast(hGam2, valGam,myColors);');
            
           
            
            while ~gammaOK
                colormap(myColors);
                x=posIm(1); y=yGamma+height;
                uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Gamma' valGam})
                x=x+width+0.01; y=yCont+height;
                uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Contrast' valContr})
                drawnow
            end
            freezeColors
            
            if strcmp(ansMethod,'Manually draw all ROIs')
                step='Curate_ROIs';
                continue
            end
           
            subplot(1,2,2)
            
            im=f;
            scaleLocCont=30; 
            minvalue=floor(scaleLocCont^2*0.1);
            maxvalue=floor(scaleLocCont^2*0.9);
            im2=f;
            B = ordfilt2(im2,minvalue,true(scaleLocCont),'symmetric');
            C = ordfilt2(im2,maxvalue,true(scaleLocCont),'symmetric');
            imReZero=im2-B;
            imReZero(imReZero<=0)=0;
            imRange=(C-B);
            imRange(imRange<=0)=eps;
            imagenorm=rescalegd(imReZero./imRange);
            imagenorm(imagenorm<0)=0; imagenorm(imagenorm>1)=1;
            
            sizeIm=min(size(f,1),size(f,2));
                        
            if 0.005*(round(sizeIm/2)-1)<1
               stepSmall=1/(round(sizeIm/2)-1);
            else
                stepSmall=0.005;
            end
            
            
            hold on
            set(gcf,'Name','Adjust Gamma', 'Toolbar','figure');
            imagesc(imagenorm); axis image; axis ij; colormap(gray)
            set(gca,'xticklabel','','yticklabel','')
            posIm2=get(gca,'Position');
            
            width=.1; height=.05;
            x=posIm2(1); y=max(posIm2(2)-2*height-0.01,0.01); yLoc=y;
            hLocCont=uicontrol('Style','slider','String','Local contrast scale','Units','normalized','position',[x y width height],'Min', 1,'Max',round(sizeIm/2),'value', scaleLocCont,'SliderStep',[stepSmall*5 .025*5],'Callback','[scaleLocCont scaleChanged]=changeScale(hLocCont);');
            
             x=posIm2(1); y=yLoc+height;
            uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Local contrast scale (pixels)' scaleLocCont})
           
            width=.1; height=.1;
            x=min(1-width+.01,posIm2(1)+posIm2(3)+.01); y=posIm2(2)+.1;
            uicontrol('Style','pushbutton','String','Done','CallBack','locContrastOK=1;step=''Regional_minima'';','Units','normalized','position',[x y width height]);
            
            scaleChanged=0;
            drawnow
            
            while ~locContrastOK
                
                if scaleChanged
                    im=f;
                    minvalue=floor(scaleLocCont^2*0.1);
                    maxvalue=floor(scaleLocCont^2*0.9);
                    im2=f;
                    width=.1; height=.05;
                    x=posIm2(1); y=yLoc+height;
                    uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Local contrast scale (pixels)' scaleLocCont})
                    
                    
                    drawnow
                    B = ordfilt2(im2,minvalue,true(scaleLocCont),'symmetric');
                    C = ordfilt2(im2,maxvalue,true(scaleLocCont),'symmetric');
                    imReZero=im2-B;
                    imReZero(imReZero<=0)=0;
                    imRange=(C-B);
                    imRange(imRange<=0)=eps;
                    imagenorm=rescalegd(imReZero./imRange);
                    imagenorm(imagenorm<0)=0; imagenorm(imagenorm>1)=1;
                    
                    hold off
                    imagesc(imagenorm);  axis image; axis ij; colormap(gray)
                    set(gca,'xticklabel','','yticklabel','')
                    x=posIm2(1); y=yLoc+height;
                    uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Local contrast scale (pixels)' scaleLocCont})
                    
                    drawnow
                    scaleChanged=0;
                else
                    
                    drawnow
                    
                end
                
            end
            
            close(hf)
            
        case 'Regional_minima'
             
            %%%%%%%%%%%%%%%% STEP 4: Regional minima %%%%%%%%%%%%%%%%%%%%%%%
            
            valBorder=.15; 
            valCenter=.1; 
            threshOK=0;
            labNuc=1;
            replot=0;
            plotROIs=0;
            bigger=1;
            hf=figure(200); clf; hold on
            set(gcf,'Name','Regional minima', 'Toolbar','figure');
            subplot(1,2,2)
            avg2=avg-min(min(avg));
            avg2=avg2/max(max(avg2));
            imshow(avg2); colormap(myColors); hold on
            posIm2=get(gca,'Position');
           
            width=.1; height=.1;
            x=posIm(1); y=min(1-height-.01,posIm(2)+posIm(4)+.03);
            uicontrol('Style','pushbutton','String','Create masks again','CallBack','step=''ROI_Mask'';threshOK=2;','Units','normalized','position',[x y width height]); 
            
            x=posIm(1)+width+.01; y=min(1-height-.01,posIm(2)+posIm(4)+.03);
            uicontrol('Style','pushbutton','String','Adjust gamma again','CallBack','step=''Adjust_Gamma'';threshOK=2;','Units','normalized','position',[x y width height]);
            width=.1; height=.05;
            x=posIm(1); y=max(posIm(2)-2*height-0.01,0.01);  yBorder=y; 
            hBorder=uicontrol('Style','slider','String','Border thresh.','Units','normalized','position',[x y width height],'Min', 0,'Max',1,...
                'value', valBorder,'sliderStep',[0.05 0.2],'Callback','valBorder=get(hBorder,''Value''); replot=1;');
            x=x+width+0.01; y=max(posIm(2)-2*height-0.01,0.01); yCenter=y;
            hCenter=uicontrol('Style','slider','String','Center thresh.','Units','normalized','position',[x y width height],'Min', 0,'Max',0.5,...
                'value', valCenter,'sliderStep',[0.01 0.05],'Callback','valCenter=get(hCenter,''Value''); replot=1;');
            
            
            x=posIm(1); y=yBorder+height;
            uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Border thresh.' valBorder},'ForegroundColor','g','BackgroundColor','w')

            x=x+width+0.01; y=yCenter+height;
            uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Cell center thresh.' valCenter},'ForegroundColor','r','BackgroundColor','w')
            
            width=.1; height=.1;
            x=max(.01,posIm(1)-width-.01); y=posIm(2)+.1;
            uicontrol('Style','pushbutton','String','Find ROIs','CallBack','replot=1; plotROIs=1;','Units','normalized','position',[x y width height]);
            x=max(.01,posIm(1)-width-.01); y=y+height+.01; yLab=y;
            uicontrol('Style','pushbutton','String','Labeled Nuclei','CallBack','labNuc=1; replot=1;','Units','normalized','position',[x y width height],'ForegroundColor','r');
            x=max(.01,posIm(1)-width-.01); y=y+height+.01; yUnlab=y;
            uicontrol('Style','pushbutton','String','Unlabeled Nuclei','CallBack','labNuc=0; replot=1;','Units','normalized','position',[x y width height]);
            
          
            
            subplot(1,2,1)
          
            
            if labNuc==1
                I=imagenorm;
                I(~mask)=0;
            else
                I=imcomplement(imagenorm);
                I(~mask)=0;
            end
            
            I(isnan(I))=0;
            mybw=im2bw(I,valBorder);
            
            CC = bwconncomp(mybw,4);
            
            test=zeros([size(f,1),size(f,2),CC.NumObjects],'uint8');
            
            for i=1:CC.NumObjects
                temp=logical(zeros(size(f)));
                temp(CC.PixelIdxList{i})=1;
                temp=imclose(temp, strel('disk',1));
                temp=bwareaopen(temp, 7);
                temp=imdilate(temp,strel('disk',1));
                
                test(:,:,i)=temp;
            end
            mybw4=sum(test,3);
            mybw4_perim = bwperim(mybw4);
            
            outsideCells=~mybw4;
            
            INuclei=I;
            INuclei(outsideCells)=0;
            

            mymask_em = imextendedmax(INuclei, valCenter);
     
      
            imshow(avg2); colormap(myColors); hold on
            drawOnTop=zeros(size(avg2,1),size(avg2,2),3);
            drawOnTop(:,:,2)=mybw4_perim;
            drawOnTop(:,:,1)=mymask_em & ~mybw4_perim;
            drawOnTop(:,:,3)=outsideCells;
            h=imshow(drawOnTop); hold off
            set(h,'AlphaData',0.3*logical(sum(drawOnTop,3))) % colored pixels on top (with transparency)          
            title(['Set parameters: {\color{red}Cell center ','\color{green}Border of considered region} \color{blue}Non-considered region'],'interpreter','tex')
           
            drawnow
             
               
            while ~threshOK
                if replot
                    if labNuc==1
                        width=.1;height=.1;
                        x=max(.01,posIm(1)-width-.01); y=yLab;
                        uicontrol('Style','pushbutton','String','Labeled Nuclei','CallBack','labNuc=1; replot=1;','Units','normalized','position',[x y width height],'ForegroundColor','r');
                        x=max(.01,posIm(1)-width-.01); y=yUnlab;
                        uicontrol('Style','pushbutton','String','Unlabeled Nuclei','CallBack','labNuc=0; replot=1;','Units','normalized','position',[x y width height]);
                        
                    else
                        width=.1;height=.1;
                        x=max(.01,posIm(1)-width-.01); y=yLab;
                        uicontrol('Style','pushbutton','String','Labeled Nuclei','CallBack','labNuc=1; replot=1;','Units','normalized','position',[x y width height]);
                        x=max(.01,posIm(1)-width-.01); y=yUnlab;
                        uicontrol('Style','pushbutton','String','Unlabeled Nuclei','CallBack','labNuc=0; replot=1;','Units','normalized','position',[x y width height],'ForegroundColor','r');
                        
                    end
                    height=.05;
                    x=posIm(1); y=yBorder+height;
                    uicontrol('Style','text','Units','normalized','position',[x y width height],'String','Please wait','ForegroundColor','k','BackgroundColor','w')
                    x=x+width+0.01; y=yCenter+height;
                    uicontrol('Style','text','Units','normalized','position',[x y width height],'String','Please wait','ForegroundColor','k','BackgroundColor','w')
                    
                    drawnow
                    subplot(1,2,1)
                    if labNuc==1
                        I=imagenorm;
                        I(~mask)=0;
                    else
                        I=imcomplement(imagenorm);
                        I(~mask)=0;
                    end
                    
                    I(isnan(I))=0;
                    mybw=im2bw(I,valBorder);
                    
                    CC = bwconncomp(mybw,4);
                    
                    test=zeros([size(f,1),size(f,2),CC.NumObjects],'uint8');
                    
                    for i=1:CC.NumObjects
                        temp=logical(zeros(size(f)));
                        temp(CC.PixelIdxList{i})=1;
                        temp=imclose(temp, strel('disk',1));
                        temp=bwareaopen(temp, 7);
                        temp=imdilate(temp,strel('disk',1));
                        
                        test(:,:,i)=temp;
                    end
                    mybw4=sum(test,3);
                    mybw4_perim = bwperim(mybw4);
                    
                    outsideCells=~mybw4;
                    
                    INuclei=I;
                    INuclei(outsideCells)=0;
                    
                  
                    mymask_em = imextendedmax(INuclei, valCenter);
                    
                    cla
                    
                    imshow(avg2); colormap(myColors); hold on
                    drawOnTop=zeros(size(avg2,1),size(avg2,2),3);
                    drawOnTop(:,:,2)=mybw4_perim;
                    drawOnTop(:,:,1)=mymask_em & ~mybw4_perim ;
                    drawOnTop(:,:,3)=outsideCells;
                    h=imshow(drawOnTop); hold off
                    set(h,'AlphaData',0.3*logical(sum(drawOnTop,3))) % colored pixels on top (with transparency)
                    title(['Set parameters: {\color{red}Cell center ','\color{green}Border of considered region} \color{blue}Non-considered region'],'interpreter','tex') 
                    x=posIm(1); y=yBorder+height;
                    uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Border thresh.' valBorder},'ForegroundColor','g','BackgroundColor','w')
                    x=x+width+0.01; y=yCenter+height;
                    uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Cell center thresh.' valCenter},'ForegroundColor','r','BackgroundColor','w')
                    
                    drawnow
                    %
                    if plotROIs
                        subplot(1,2,2)
                        title('Please wait...','Color','k')
                        drawnow
                        width=.08; height=.1;
                        x=min(1-width,posIm2(1)+posIm2(3)+.01); y=posIm2(2)+.1; yDone=y;
                        uicontrol('Style','pushbutton','String','Done','CallBack','threshOK=1;step=''Find_automatic_ROIs'';','Units','normalized','position',[x y width height]);
                        if bigger==1
                            y=yDone+height+.01;
                            uicontrol('Style','pushbutton','String','Bigger ROIs','CallBack','bigger=1; replot=1; plotROIs=1;','Units','normalized','position',[x y width height],'ForegroundColor','r');
                            y=yDone+2*(height+.01);
                            uicontrol('Style','pushbutton','String','Smaller ROIs','CallBack','bigger=0; replot=1; plotROIs=1;','Units','normalized','position',[x y width height]);
                        else
                            y=yDone+height+.01;
                            uicontrol('Style','pushbutton','String','Bigger ROIs','CallBack','bigger=1; replot=1; plotROIs=1;','Units','normalized','position',[x y width height]);
                            y=yDone+2*(height+.01);
                            uicontrol('Style','pushbutton','String','Smaller ROIs','CallBack','bigger=0; replot=1; plotROIs=1;','Units','normalized','position',[x y width height],'ForegroundColor','r');
                        end
                        
                        INuclei_c = imcomplement(INuclei);
                        
                                
                        a=(~mybw4_perim & ~bwperim(imerode(mybw4,strel('disk',1))));
                     
                        
                        L = bwlabel(~a);
                        props=regionprops(~a,'FilledArea','Area');
                        idx=find([props.Area]== [props.FilledArea]);
                        bw2 = ismember(L, idx);
                        
                        I_mod = imimposemin(INuclei_c, outsideCells | mymask_em & ~(~(a+bw2) | mybw4_perim));
                        
                        L = watershed(I_mod);
                        em=L==0;
                        CC = bwconncomp(~em,4);
                        
                        test=zeros(size(f));
                        count=1;
                        clear points CCFilt
                        
                        
                        
                        imagesc(avg2); hold on; axis image; axis ij
                        set(gca,'xticklabel','','yticklabel','')
                         
                        steps = CC.NumObjects;
                      
                        for i=1:CC.NumObjects
                            temp=logical(zeros(size(f)));
                         
                            temp(CC.PixelIdxList{i})=1;
                            if length(CC.PixelIdxList{i})<5  | any(any(temp & outsideCells))
                               continue 
                            end
                            if bigger
                                temp=bwmorph(temp,'thicken',1); 
                            end
                            test = test | temp;
                            
                            CCFilt.PixelIdxList{count}=find(temp);
                            CCFilt.NumObjects=count;
                            points{count}=bwboundaries(temp);
                            lp=line(points{count}{1}(:,2),points{count}{1}(:,1));
                            set(lp,'color',[1 0 0],'LineWidth',1);
                            count=count+1;
                        end
                        title('Found ROIs','Color','r')
                       plotROIs=0;
                    end
                    
                   
                    drawnow
                    replot=0;
                else
                    drawnow
                end
            end
            if threshOK==1
                CCFilt.Connectivity=CC.Connectivity; CCFilt.ImageSize=CC.ImageSize;
                stats=regionprops(CCFilt,'Area','Perimeter','Solidity');
                circularity=([stats(:).Perimeter].^2)./(4*pi*[stats(:).Area]);
                areas=[stats.Area];
                solid=[stats.Solidity];
                intensity=zeros(1,CCFilt.NumObjects);
                for i=1:CCFilt.NumObjects
                	intensity(i)=mean(avg2(CCFilt.PixelIdxList{i}));
                
                end
            
            end
            
            
        case 'Find_automatic_ROIs'
            %%%%%%%%%%%%%%% STEP 5: Find automatic ROIs %%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%% PARAMETERS 2 %%%%%%%%%%%%%%%%%
            
   
            done=0;
            replot=0;
           
            valAreasLow=5;
            valAreasHigh=95;
            valIntensityLow=5;
            valIntensityHigh=95;
            valCircularityHigh=95;
            valCircularityLow=5;
            hf=figure(200); clf; hold on
            set(gcf,'Name','Find automatic ROIs', 'Toolbar','figure');
            imagesc(avg2); hold on; colormap(myColors); axis image; axis ij
            set(gca,'xticklabel','','yticklabel','')
            posIm=get(gca,'Position');
            width=.1; height=.1;
            x=max(.01,posIm(1)-width-.01); y=posIm(2)+.1; 
            uicontrol('Style','pushbutton','String','Create masks again','CallBack','step=''ROI_Mask'';done=1; clear first_autoROI ','Units','normalized','position',[x y width height]);
            y=y+height+.01;
            uicontrol('Style','pushbutton','String','Adjust gamma again','CallBack','step=''Adjust_Gamma'';done=1;clear first_autoROI','Units','normalized','position',[x y width height]);
            
            if ~strcmp(ansMethod,'Manually draw all ROIs')
                width=.15;
                y=y+height+.01;
                uicontrol('Style','pushbutton','String','Find Regional minima again','CallBack','step=''Regional_minima'';done=1;clear first_autoROI','Units','normalized','position',[x y width height]);
            end
            
            width=.1; height=.1;
            x=min(1-width,posIm(1)+posIm(3)+.005); y=posIm(2);
            uicontrol('Style','pushbutton','String','Done','CallBack','done=1;step=''Curate_ROIs'';','Units','normalized','position',[x y width height]);           
            y=y+height+0.01; yCircH=y;
            height=.05;
            hthCircularityHigh=uicontrol('Style','slider','BackgroundColor','g','Units','normalized','position',[x y width height],'Min', 75,'Max',100,...
                'value', valCircularityHigh,'sliderStep',[0.04 0.16],'Callback','valCircularityHigh=get(hthCircularityHigh,''Value''); replot=1;');
            y=y+2*height+0.01; yCircL=y;
            hthCircularityLow=uicontrol('Style','slider','BackgroundColor','y','Units','normalized','position',[x y width height],'Min', 0,'Max',25,...
                'value', valCircularityLow,'sliderStep',[0.04 0.16],'Callback','valCircularityLow=get(hthCircularityLow,''Value''); replot=1;');
            y=y+2*height+0.01; yIntH=y;
            hthIntensityHigh=uicontrol('Style','slider','BackgroundColor','b','Units','normalized','position',[x y width height],'Min', 75,'Max',100,...
                'value', valIntensityHigh,'sliderStep',[0.04 0.16],'Callback','valIntensityHigh=get(hthIntensityHigh,''Value''); replot=1;');
            y=y+2*height+0.01; yIntL=y;
            hthIntensityLow=uicontrol('Style','slider','BackgroundColor','c','Units','normalized','position',[x y width height],'Min', 0,'Max',25,...
                'value', valIntensityLow,'sliderStep',[0.04 0.16],'Callback','valIntensityLow=get(hthIntensityLow,''Value''); replot=1;');
            y=y+2*height+0.01; yAreaH=y;
            hthAreasHigh=uicontrol('Style','slider','BackgroundColor','r','Units','normalized','position',[x y width height],'Min', 75,'Max',100,...
                'value', valAreasHigh,'sliderStep',[0.04 0.16],'Callback','valAreasHigh=get(hthAreasHigh,''Value''); replot=1;');
            y=y+2*height+0.01; yAreaL=y;
            hthAreasLow=uicontrol('Style','slider','BackgroundColor','w','Units','normalized','position',[x y width height],'Min', 0,'Max',25,...
                'value', valAreasLow,'sliderStep',[0.04 0.16],'Callback','valAreasLow=get(hthAreasLow,''Value''); replot=1;');
            
            height=0.1;
            y=y+height+0.01;
            uicontrol('Style','pushbutton','String','Run','CallBack','replot=2;','Units','normalized','position',[x y width height]);
            
            
          
            threshAreasLow=prctile(areas,valAreasLow);
            threshAreasHigh=prctile(areas,valAreasHigh);
            threshIntensityLow=prctile(intensity,valIntensityLow);
            threshIntensityHigh=prctile(intensity,valIntensityHigh);
            threshCircularityHigh=prctile(circularity,valCircularityHigh);
            threshCircularityLow=prctile(circularity,valCircularityLow);
           
            
            for i=1:CCFilt.NumObjects
                
                lp=line(points{i}{1}(:,2),points{i}{1}(:,1));
                set(lp,'color',[1 0 0],'LineWidth',1);
                
            end
            title('Found ROIs (no morphological criteria imposed)','Color','r')
            CCFilt2=CCFilt;
            pointsFilt=points;
            bkgTemp=test;
            
            height=0.05;
            y=yAreaL+height;
            uicontrol('Style','text','BackgroundColor','w','Units','normalized','position',[x y width height],'String',{'Minimal area' round(threshAreasLow)})
            y=yAreaH+height;
            uicontrol('Style','text','BackgroundColor','r','Units','normalized','position',[x y width height],'String',{'Maximal area' round(threshAreasHigh)})
            y=yIntL+height;
            uicontrol('Style','text','BackgroundColor','c','Units','normalized','position',[x y width height],'String',{'Minimal Intensity' threshIntensityLow})
            y=yIntH+height;
            uicontrol('Style','text','BackgroundColor','b','Units','normalized','position',[x y width height],'String',{'Maximal Intensity' threshIntensityHigh})
            y=yCircL+height;
            uicontrol('Style','text','BackgroundColor','y','Units','normalized','position',[x y width height],'String',{'Minimal Circ.' threshCircularityLow})
            y=yCircH+height;
            uicontrol('Style','text','BackgroundColor','g','Units','normalized','position',[x y width height],'String',{'Max. Circ.' threshCircularityHigh})
            drawnow
           
            
            while ~done
                
                if replot
                    threshAreasLow=prctile(areas,valAreasLow);
                    
                    threshAreasHigh=prctile(areas,valAreasHigh);
                    threshIntensityLow=prctile(intensity,valIntensityLow);
                    threshIntensityHigh=prctile(intensity,valIntensityHigh);
                    threshCircularityHigh=prctile(circularity,valCircularityHigh);
                    threshCircularityLow=prctile(circularity,valCircularityLow);
                    y=yAreaL+height;
                    uicontrol('Style','text','BackgroundColor','w','Units','normalized','position',[x y width height],'String',{'Minimal area' round(threshAreasLow)})
                    y=yAreaH+height;
                    uicontrol('Style','text','BackgroundColor','r','Units','normalized','position',[x y width height],'String',{'Maximal area' round(threshAreasHigh)})
                    y=yIntL+height;
                    uicontrol('Style','text','BackgroundColor','c','Units','normalized','position',[x y width height],'String',{'Minimal Intensity' threshIntensityLow})
                    y=yIntH+height;
                    uicontrol('Style','text','BackgroundColor','b','Units','normalized','position',[x y width height],'String',{'Maximal Intensity' threshIntensityHigh})
                    y=yCircL+height;
                    uicontrol('Style','text','BackgroundColor','y','Units','normalized','position',[x y width height],'String',{'Minimal Circ.' threshCircularityLow})
                    y=yCircH+height;
                    uicontrol('Style','text','BackgroundColor','g','Units','normalized','position',[x y width height],'String',{'Max. Circ.' threshCircularityHigh})
                    drawnow
                    if replot==2
                        title('Please wait...','Color','k'); drawnow
                        first_autoROI=1;
                        threshAreasLow=prctile(areas,valAreasLow);
                        threshAreasHigh=prctile(areas,valAreasHigh);
                        threshIntensityLow=prctile(intensity,valIntensityLow);
                        threshIntensityHigh=prctile(intensity,valIntensityHigh);
                        threshCircularityHigh=prctile(circularity,valCircularityHigh);
                        threshCircularityLow=prctile(circularity,valCircularityLow);
                        
                        
                        idx= find((threshAreasLow<= areas) & (areas <= threshAreasHigh) & (threshIntensityLow<= intensity) & (intensity <= threshIntensityHigh)  & (threshCircularityLow<= circularity) & (circularity <= threshCircularityHigh));
                        
                        bkgTemp=zeros(size(avg2));
                        count=1;
                        clear pointsFilt CCFilt2
                        
                        imagesc(avg2); hold on; colormap(myColors); axis image; axis ij
                        set(gca,'xticklabel','','yticklabel','')
                        
                        
                        for k=1:length(idx)
                            i=idx(k);
                            temp=logical(zeros(size(avg2)));
                            if length(CCFilt.PixelIdxList{i})==1
                                continue
                            end
                            temp(CCFilt.PixelIdxList{i})=1;
                            
                            bkgTemp = bkgTemp | temp;
                            
                            CCFilt2.PixelIdxList{count}=find(temp);
                            CCFilt2.NumObjects=count;
                            pointsFilt{count}=bwboundaries(temp);
                            lp=line(pointsFilt{count}{1}(:,2),pointsFilt{count}{1}(:,1));
                            set(lp,'color',[1 0 0],'LineWidth',1);
                            count=count+1;
                        end
                        title('ROIs after imposing filtering criteria','Color','k')
                        
                        
                    end
                    replot=0;
                else
                    drawnow
                end
                
            end
            CC2=CCFilt2;
            temp=bkgTemp;
            points2=pointsFilt;
        case 'Curate_ROIs'
            %%%%%%%%%%%%%%% STEP 6: Curate ROIs %%%%%%%%%%%%%%%%%%%%%%%%
            
            
            deleteAreaROI=0;
            deleteROI=0;
            restoreROI=0;
            addROI=0;
            hideROIs=0;
            
            okROIs=0;
            cont=0;
            
            hf=figure(200); clf; hold on
            set(gcf,'Name','Curate ROIs', 'Toolbar','figure');
            
             imagesc(avg); 
            
            
            
            colormap(myColors);  axis image; axis ij; hold on
            set(gca,'xticklabel','','yticklabel','')
            posIm=get(gca,'Position');
            for i=1:length(points2)
                lp=line(points2{i}{1}(:,2),points2{i}{1}(:,1));
                set(lp,'color',[1 0 0],'LineWidth',1);
                
            end
            
            width=.1; height=.1;
            x=max(.01,posIm(1)-width-.01); y=posIm(2);
            uicontrol('Style','pushbutton','String','Continue','CallBack','uiresume(gcbf)','Units','normalized','position',[x y width height]);
            y=y+height+.01;
            uicontrol('Style','pushbutton','String','Pause','CallBack','uiwait(gcf)','Units','normalized','position',[x y width height]);
            width=.15;
            if ~strcmp(ansMethod,'Manually draw all ROIs')
                y=y+height+.01;
                uicontrol('Style','pushbutton','String','Find automatic ROIs again','CallBack','step=''Find_automatic_ROIs'';okROIs=2;clear first_autoROI first_okROI',...
                    'Units','normalized','position',[x y width height]);
                y=y+height+.01;
                uicontrol('Style','pushbutton','String','Find Regional minima again','CallBack','step=''Regional_minima'';okROIs=2;clear first_autoROI first_okROI',...
                    'Units','normalized','position',[x y width height]);            
            end
            y=y+height+.01;
            uicontrol('Style','pushbutton','String','Adjust gamma again','CallBack','step=''Adjust_Gamma'';okROIs=2;clear first_autoROI first_okROI',...
                'Units','normalized','position',[x y width height]);
            y=y+height+.01;
            uicontrol('Style','pushbutton','String','Create masks again','CallBack','step=''ROI_Mask'';okROIs=2; clear first_autoROI first_okROI',...
                'Units','normalized','position',[x y width height]);
            width=.1; height=.05;
            y=y+2*height+0.01; yCont=y;
            hGam2=uicontrol('Style','slider','String','Set Contrast','Units','normalized','position',[x y width height],'Min',1,'Max',900,...
                'value', valContr,'SliderStep',[.1 .5],'Callback','[myColors valContr]=changeContrast(hGam2, valGam,myColors);');
            y=y+2*height+0.01; yGam=y;
            hGam1=uicontrol('Style','slider','String','Set Gamma','Units','normalized','position',[x y width height],'Min', 0.005,'Max',2,...
                'value', valGam,'SliderStep',[.01 .05],'Callback','[myColors valGam]=changeGamma(hGam1, valContr,myColors);');
            
            width=.1; height=.1;
            x=min(1-width,posIm(1)+posIm(3)+.01); y=posIm(2);
            uicontrol('Style','pushbutton','String','Done','CallBack','okROIs=1;step=''DF_save'';','Units','normalized','position',[x y width height]);
            y=y+height+0.01;
            uicontrol('Style','pushbutton','String','Hide ROIs','CallBack','hideROIs=1;','Units','normalized','position',[x y width height]);
            y=y+1.3*height+0.01;
            uicontrol('Style','pushbutton','String','Stop Draw ROI','CallBack','addROI=2;','Units','normalized','position',[x y width height]);
            y=y+height+0.01;
            uicontrol('Style','pushbutton','String','Draw ROIs','CallBack','addROI=1;','Units','normalized','position',[x y width height]);
            y=y+1.3*height+0.01;
            uicontrol('Style','pushbutton','String','Stop Delete ROIs','CallBack','deleteROI=2','Units','normalized','position',[x y width height]);
            y=y+height+0.01;
            uicontrol('Style','pushbutton','String','Delete ROIs','CallBack','deleteROI=1;','Units','normalized','position',[x y width height]);
            
            y=y+1.3*height+0.01;
            uicontrol('Style','pushbutton','String','Delete ROIs Area','CallBack','deleteAreaROI=1;','Units','normalized','position',[x y width height]);
            
           
            while okROIs==0
          
                figure(200)
                title('Waiting input')
                colormap(myColors);
                width=.1; height=.05;
                x=max(.01,posIm(1)-width-.01); y=yGam+height;
                uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Gamma' valGam})
                y=yCont+height;
                uicontrol('Style','text','Units','normalized','position',[x y width height],'String',{'Contrast' valContr})
                
                drawnow
                
                if cont==1
                    imagesc(avg);
                    colormap(myColors);  axis image; axis ij; hold on
                    set(gca,'xticklabel','','yticklabel','')
                    for i=1:length(points2)
                        lp=line(points2{i}{1}(:,2),points2{i}{1}(:,1));
                        set(lp,'color',[1 0 0],'LineWidth',1);
                        
                    end
                    hfig=imcontrast;
                    uiwait(hfig)
                    cont=0;
                end
                
                if hideROIs
                    figure(200);
                    
                    imagesc(avg);
                    colormap(myColors);  axis image; axis ij; hold on
                    set(gca,'xticklabel','','yticklabel','')
                    title('ROIs temporally hidden')
                    drawnow
                    drawnow
                    for i=1:length(points2)
                        lp=line(points2{i}{1}(:,2),points2{i}{1}(:,1));
                        set(lp,'color',[1 0 0],'LineWidth',1);
                        
                    end
                    hideROIs=0;
                end
                
                if deleteAreaROI==1 
                    if isempty(CC2)
                        title('No ROIs to delete')
                        deleteAreaROI=0;
                        pause(2);
                    else
                        title('Draw a mask that encloses the ROIs with left-clicks, then right-click inside it and select Create mask')
                        wpoly=roipoly;
                        indToRemove=[];
                        for i=1:CC2.NumObjects
                            if sum(ismember(find(temp&wpoly),CC2.PixelIdxList{i}))>0.5*length(CC2.PixelIdxList{i})
                                indToRemove=[indToRemove,i];
                            end
                        end
                        points2(indToRemove)=[];
                        for i=1:length(indToRemove)
                            temp(CC2.PixelIdxList{indToRemove(i)})=0;
                        end
                        CC2.PixelIdxList(indToRemove)=[];
                        CC2.NumObjects=CC2.NumObjects-length(indToRemove);
                        
                        figure(200);
                        
                        imagesc(avg);
                        colormap(myColors);  axis image; axis ij; hold on
                        set(gca,'xticklabel','','yticklabel','')
                        for i=1:length(points2)
                            lp=line(points2{i}{1}(:,2),points2{i}{1}(:,1));
                            set(lp,'color',[1 0 0],'LineWidth',1);
                        end
                        deleteAreaROI=0;
                    end
                end
                    
                
                
                if deleteROI ==1
                    if isempty(CC2)
                        title('No ROIs to delete')
                        deleteROI=0;
                        pause(2);
                    else
                        title('Mark ROIs for deletion with right-clicks')
                        
                        [colToDel, rowToDel] = getpts(gcf);
                        line([round(colToDel) round(colToDel)],[round(rowToDel)-2 round(rowToDel)+2],'color','y','linewidth',2)
                        line([round(colToDel)-2 round(colToDel)+2],[round(rowToDel) round(rowToDel)],'color','y','linewidth',2)
                        if deleteROI ==1
                            indxROIToDel=sub2ind(size(f),round(rowToDel),round(colToDel));
                            if ismember(indxROIToDel, find(temp));
                                for j=1:CC2.NumObjects
                                    if ismember(indxROIToDel,CC2.PixelIdxList{j})
                                        temp(CC2.PixelIdxList{j})=0;
                                        points2(j)=[];
                                        CC2.PixelIdxList(j)=[];
                                        CC2.NumObjects=CC2.NumObjects-1;
                                        
                                        break
                                    end
                                end
                            else
                                [colToDel, rowToDel] = getpts(hf);
                                
                                indxROIToDel=sub2ind(size(f),round(rowToDel),round(colToDel));
                            end
                            
                        elseif deleteROI==2
                            
                            deleteROI=0;
                            figure(200);
                            imagesc(avg);
                            colormap(myColors);  axis image; axis ij; hold on
                            set(gca,'xticklabel','','yticklabel','')
                            for i=1:length(points2)
                                lp=line(points2{i}{1}(:,2),points2{i}{1}(:,1));
                                set(lp,'color',[1 0 0],'LineWidth',1);
                                
                            end
                        end
                    end
                  
                end
              
                if addROI==1
                    title('Draw ROI perimeters with left-clicks and close each perimeter with a final right-click')
                    addedROI=roipoly_adv;
                    if addROI==1
                        pixsROI=find(addedROI);
                        pixsToDel=[];
                        for i=1:length(pixsROI)
                            clear rectFilt
                            [rowToCheck,colToCheck]=ind2sub(size(addedROI),pixsROI(i));
                            rect=[rowToCheck-1:rowToCheck+1;colToCheck-1:colToCheck+1];
                            rectFilt(1,:)=rect(1,rect(1,:)>0 & rect(1,:)<=size(temp,1));
                            rectFilt(2,:)=rect(2,rect(2,:)>0 & rect(2,:)<=size(temp,2));
                            checkROI=zeros(3,3);
                            checkROI(rectFilt(1,:)-rowToCheck+2,rectFilt(2,:)-colToCheck+2)=temp(rectFilt(1,:),rectFilt(2,:));
                            neigbours=find(checkROI);
                            if any(ismember(neigbours,[2 4 5 6 8]))
                                pixsToDel=[pixsToDel i];
                            end
                        end
                        addedROI(pixsROI(pixsToDel))=0;
                        
                        temp=temp|addedROI;
                        points2{end+1}=bwboundaries(addedROI);
                        
                        lp=line(points2{end}{1}(:,2),points2{end}{1}(:,1)); set(lp,'color',[1 0 0],'LineWidth',1);
                        
                        if isempty(CC2)
                            CC2.PixelIdxList{1}=find(addedROI);
                            CC2.NumObjects=1;
                        else
                            
                            CC2.PixelIdxList{end+1}=find(addedROI);
                            CC2.NumObjects=CC2.NumObjects+1;
                        end
                        drawnow
                        
                    elseif addROI==2
                        addROI=0;
                    end
                end
                
                
            end
            
           
            avg=fOld;
            
            cell_per=cell(CC2.NumObjects,1);
            cells=cell(CC2.NumObjects,1);
           
            temp2=zeros(size(avg));
            for i=1:CC2.NumObjects
                imgROI=zeros(size(avg));
                imgROI(CC2.PixelIdxList{i})=1;
                cell_per{i}=fliplr(points2{i}{1});
                cells{i}=CC2.PixelIdxList{i};
                temp2=temp2|imgROI;
            end
            cell_number=CC2.NumObjects;
           
            bkg=temp;
            close(hf)
            
        case 'DF_save'
            %%%%%%%%%%%%% STEP 7: Obtain fluorescence traces of cells and save %%%%%%%%%%%%%%%%%%
            
            cells_mean=zeros(numFrames,cell_number);
           
           
            for i=1:cell_number
                pixs=video(:,cells{i});
                cells_mean(:,i)=mean(pixs,2);
            end
            
            end_pro=1;
            
            
    end
end

for i=1:cell_number
    img=zeros(size(avg));
    img(cells{i})=1;
    junk=regionprops(img,'Centroid');
    centers{i}=[junk.Centroid(1) junk.Centroid(2)];
end
distances=zeros(cell_number);
for i=1:cell_number
    for j=1:cell_number
        distances(i,j)=sqrt(((centers{i}(1)-centers{j}(1))*pixelLengthX)^2+((centers{i}(2)-centers{j}(2))*pixelLengthY)^2);
    end
end

wdthNpil=20;

[columnsInImage rowsInImage] = meshgrid((1:size(avg,2))/pixelLengthX, (1:size(avg,1))/pixelLengthY);

npil_mean=zeros(numFrames,length(cells));
bkgDilated=imdilate(logical(sum(bkg,3)),strel('disk',2));
for i=1:length(cells)
    
    circlePixels = (rowsInImage - centers{i}(2)/pixelLengthY).^2 + (columnsInImage - centers{i}(1)/pixelLengthX).^2 <= wdthNpil.^2;
    
    npil=zeros(size(avg));
    npil(circlePixels)=1;
    npil(bkgDilated)=0;
    pixs=video(:,logical(npil));
    npil_mean(:,i)=mean(pixs,2);
end



cutName=strfind(fileName,'.tif');
outputFile=[fileName(1:cutName-1) '_ALL_CELLS.mat'];

save(outputFile,'avg','bkg','cell_number','cells_mean','npil_mean','cell_per','cells','distances','pixelLengthX','pixelLengthY');

