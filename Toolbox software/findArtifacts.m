function findArtifacts()


[filename,pathname] = uigetfile({'*.tif';'*.TIFF'},'Open file', 'MultiSelect', 'off');
fileName=fullfile(pathname,filename);

scrsz=get(0,'ScreenSize'); text_si=18;

set(0,'DefaultAxesFontSize',text_si,'DefaultFigureColor','w', 'DefaultAxesTickDir', 'in','DefaultFigureWindowStyle','normal',...
    'DefaultFigurePosition', [1 1 scrsz(3) scrsz(4)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v = ver;
% parOn=any(strcmp('Parallel Computing Toolbox', {v.Name}));
parOn=0;
%parOn=0; % UNCOMMENT THIS FOR PARALLELIZATION
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


disp('Reading movie file... done');
%mean,  filter and variable definition
f=squeeze(mean(video,1));
% cut video
hf=figure;
set(hf,'Name','Draw rectangular mask to remove borders');
imagesc(f);
set(gca,'xticklabel','','ytickLabel','')
axis image
colormap(gray)
drawnow
h=imrect(gca);
rect_coord=getPosition(h);
close(hf);
rowsRect=round(rect_coord(2)):round(rect_coord(4)+rect_coord(2));
colsRect=round(rect_coord(1)):round(rect_coord(3)+rect_coord(1));
f=f(rowsRect,colsRect);
video_cut=video(:,rowsRect,colsRect);

clear video

% mask
hf=figure;
set(hf,'Name','Draw a mask');
imagesc(f)
set(gca,'xticklabel','','ytickLabel','')
axis image
colormap(gray)
drawnow
mask=roipoly;mask_backup=mask;


mean_cut=mean(video_cut);
mean_cut(~mask)=0;
min_video=min(min(mean_cut));
max_video=max(max(mean_cut));



close(hf);
% template
hf=figure;
set(hf,'Name','Draw rectangular mask to define template');
imagesc(f)
set(gca,'xticklabel','','ytickLabel','')
axis image
colormap(gray)
drawnow
h=imrect;
rect_template=getPosition(h);
close(hf);
rowsTemplate=round(rect_template(2)):round(rect_template(4)+rect_template(2));
colsTemplate=round(rect_template(1)):round(rect_template(3)+rect_template(1));

hfilter1=fspecial('gaussian',[20 20],1);
fFilt1=imfilter(f,hfilter1);
fFilt1(~mask)=0;

peak=zeros(numFrames,1);
peakAv=zeros(numFrames,1);

%% set threshold by maximizing correlation with template
thres=0.0025;
vect_step=-5:95;
test_frames=[round(0.1*numFrames) round(0.5*numFrames) round(0.9*numFrames)];
peak_adjust=zeros(length(vect_step),length(test_frames));
for FR=1:length(test_frames)
    for n=1:length(vect_step)
        thres_adjust=thres+vect_step(n)*0.01*(0.05-0.0001);
        if thres_adjust<0
            continue
        else 
            test=im2single(im2bw(squeeze(video_cut(test_frames(FR),:,:)),thres_adjust));
            test=test(rowsTemplate,colsTemplate);
            template1=imfilter(test,hfilter1);
            template1(~mask(rowsTemplate,colsTemplate))=0;
           
            if all(reshape(template1,1,[])==template1(1,1))
                C1=0;
            else
                
                C1 = normxcorr2(template1, fFilt1);
            end
            peak_adjust(n,FR)=C1(round(size(C1,1)/2), round(size(C1,2)/2));
        end
    end
end

[pks, locs]=max(mean(peak_adjust'));
thres=thres+vect_step(locs)*0.01*(0.05-0.0001);


test=im2single(im2bw(squeeze(video_cut(round(0.1*numFrames),:,:)),thres));
test=test(rowsTemplate,colsTemplate);
template1=imfilter(test,hfilter1);


%% perform calculation
mask=mask(rowsTemplate,colsTemplate);
ind=[];

if parOn
    parfor i=1:numFrames
        test=im2single(im2bw(squeeze(video_cut(i,:,:)),thres));
        test=test(rowsTemplate,colsTemplate);
        template1=imfilter(test,hfilter1);
        template1(~mask)=0;
        if ~any(any(template1))
            ind=[ind i];
            continue
        end
        C1 = normxcorr2(template1, fFilt1);
        peak(i)=C1(round(size(C1,1)/2), round(size(C1,2)/2));
        peakAv(i)=sum(sum(C1(round(size(C1,1)/2-5:size(C1,1)/2+5), round(size(C1,2)/2-5:size(C1,2)/2+5))));
    end
else
    for i=1:numFrames
        test=im2single(im2bw(squeeze(video_cut(i,:,:)),thres));
        test=test(rowsTemplate,colsTemplate);
        template1=imfilter(test,hfilter1);
        template1(~mask)=0;
        if ~any(any(template1))
            ind=[ind i];
            continue
        end
        C1 = normxcorr2(template1, fFilt1);
        peak(i)=C1(round(size(C1,1)/2), round(size(C1,2)/2));
        peakAv(i)=sum(sum(C1(round(size(C1,1)/2-5:size(C1,1)/2+5), round(size(C1,2)/2-5:size(C1,2)/2+5))));
    end
    
end

for i=1:length(ind)
    peak(ind(i))=peak(ind(i)-1);
    peakAv(ind(i))=peakAv(ind(i)-1);
end

%% find artifacts (movements)
movements=zeros(numFrames,1);
step_runline=50;
thres_zscore=-3;
z_score_peak=zscore(peak-runline(peak,step_runline,1));
z_score_peakAv=zscore(peakAv-runline(peakAv,step_runline,1));
movements(z_score_peak<thres_zscore)=1;
movements(ind)=1;
hf=figure;
set(gcf,'Name', ['threshold value: ' num2str(thres_zscore) ', artifacts found: ' num2str(sum(movements)) '; press enter to continue'])
plot(z_score_peakAv)
line([1 numFrames], [thres_zscore thres_zscore],'color','g')
ylim([-7 7])
title('peakAv')
pause
ok_thresh_zscore=0;
while ~ ok_thresh_zscore
    choice = questdlg([ 'Artifacts found: ' num2str(sum(movements)) ', threshold value, ' num2str(thres_zscore) ', Are you satisfied with threshold value?'], ...
        'Threshold zscore', 'Yes','No','Yes');
    switch choice
        case 'Yes'
            ok_thresh_zscore=1;
            %             close(hf)
        case 'No'
            inputTitle = char('threshold zscore?');
            prompt = {'new threshold'};
            answer = inputdlg(prompt, inputTitle, 1);
            thres_zscore= str2num(answer{1});
            clf(hf); figure(hf)
            movements=zeros(numFrames,1);
            movements(z_score_peak<thres_zscore)=1;
            movements(ind)=1;
            set(gcf,'Name', ['threshold value: ' num2str(thres_zscore) ', artifacts found: ' num2str(sum(movements)) '; press enter to continue'])
            plot(z_score_peakAv)
            line([1 numFrames], [thres_zscore thres_zscore],'color','g')
            ylim([-7 7])
            title('peakAv')
            pause
    end
end
close all


%% clean artifacts manually

ind_movements=find(movements');
DIS_MIN_GR=5;
clear group_mov ind_group_mov
group_mov=1;
ind_group_mov{group_mov}=ind_movements(1);
for i=2:length(ind_movements)
    if ind_movements(i)-ind_movements(i-1)<DIS_MIN_GR
        ind_group_mov{group_mov}=[ind_group_mov{group_mov} ind_movements(i)];
    else
        group_mov=group_mov+1;
        ind_group_mov{group_mov}=ind_movements(i);
    end
end
clean_mov=cell(1,length(ind_group_mov));

dir=0;
windows_frame=5;
hf=figure;
allIn=0;
for quan_mov=1:group_mov
    clf
    subplot(1,2,1)
    posIm=get(gca,'Position');
    width=.1; height=.1;
    x=posIm(1); y=max(posIm(2)-2*height-0.01,0.01);
    uicontrol('Style','pushbutton','String','Pause','CallBack','uiwait(gcf)','Units','normalized','Position',[x y width height]);
    x=x+width+0.01;
    uicontrol('Style','pushbutton','String','Continue','CallBack','uiresume(gcbf)','Units','normalized','Position',[x y width height]);
    x=x+1.5*width;
    uicontrol('Style','pushbutton','String','Accept all','CallBack',{@break_function},'Units','normalized','Position',[x y width height]);
   
    
    x=posIm(1)+posIm(3)+.01; y=posIm(2)+.1;
    uicontrol('Style','pushbutton','String','<< Previous frame','CallBack',{@dir_back_function},'Units','normalized','Position',[x y width height]);
    x=x+width+0.01;
    uicontrol('Style','pushbutton','String','Next frame >>','CallBack',{@dir_foward_function},'Units','normalized','Position',[x y width height]);
    x=posIm(1)+posIm(3)+.01; y=y+height+0.01;
    uicontrol('Style','pushbutton','String','Yes','CallBack',{@mov_ok_function},'Units','normalized','Position',[x y width height]);
    x=x+width+0.01;
    uicontrol('Style','pushbutton','String','No','CallBack',{@mov_no_ok_function},'Units','normalized','Position',[x y width height]);
 
   
    if ind_group_mov{quan_mov}(1)>windows_frame
        START=ind_group_mov{quan_mov}(1)-windows_frame;
    else
        START=1;
    end
    if ind_group_mov{quan_mov}(end)+windows_frame<numFrames
        END=ind_group_mov{quan_mov}(end)+windows_frame;
    else
        END=numFrames;
    end
    
    clear video
    video=zeros(END-START,tiffInfo(1).Height,tiffInfo(1).Width, 'single');
    for i=START:END
        video(i-START+1,:,:)=im2single(imread(fileName,i));
    end
    video=video(:,rowsRect,colsRect);
    video(:,~mask_backup)=0;
    width=.2;height=.1;
    xTxt=posIm(1)+posIm(3)+.01;yTxt=y+height+.01;
    height=.05;
    yTxt2=yTxt+height+.01;
    while length(clean_mov{quan_mov})~= length(ind_group_mov{quan_mov})
        if allIn
            break
        end
        for j=1:length(ind_group_mov{quan_mov})
            mov_ok=[];
            FRA_TO_SHOW=ind_group_mov{quan_mov}(j)-START+1;
            dir=0;
            if allIn
                break
            end
            while isempty(mov_ok)
                if allIn
                   break 
                end
                if (dir==-1 & FRA_TO_SHOW > 1) | (dir==1 & FRA_TO_SHOW < size(video,1))
                    FRA_TO_SHOW=FRA_TO_SHOW+dir;
                end
                
                dir=0;
               subplot(1,2,1)
                imagesc(squeeze(video(FRA_TO_SHOW,:,:)))
                set(gca,'xticklabel','','ytickLabel','')
                axis image
                colormap(gray)
                caxis([min_video max_video*1.2])
                if ismember(START-1+FRA_TO_SHOW, ind_group_mov{quan_mov})
                    col_frame='r';
                else
                    col_frame='k';
                end
                title(['Frame number ' num2str(START-1+FRA_TO_SHOW)],'color',col_frame)
               
                uicontrol('Style','text','Units','normalized','position',[xTxt yTxt width height],...
                    'String',{'Does frame ' num2str(ind_group_mov{quan_mov}(j)) ' contain an artifact?'})
             
                uicontrol('Style','text','Units','normalized','position',[xTxt yTxt2 width height],...
                    'String',{['Checking frame ' num2str(j) ' (out of ' num2str(length(ind_group_mov{quan_mov})) ')'] ...
                    ['from event ' num2str(quan_mov) ' (out of ' num2str(group_mov) ')']})
               
                
                drawnow
                clean_mov{quan_mov}=[clean_mov{quan_mov} mov_ok];
            end
        end
    end
end
close(hf)

if allIn
   movements=ones(numFrames,1);
else
    movements_original=movements;
    for quan_mov=1:group_mov
        movements(ind_group_mov{quan_mov})=clean_mov{quan_mov};
    end
end

inputTitle = char('Inform any other artifact');
prompt = {['Are there any other frames that you want to label as artifacts? If so, which ones? (numbers separeted by spaces):']};
answer = inputdlg(prompt, inputTitle, 1);
mov_to_add= str2num(answer{1});

if mov_to_add
    movements(mov_to_add)=1;
end
%% save all relevant variables in *_ARTIFACTS.mat file
cut=strfind(fileName,'.tif');
outputFile=[fileName(1:cut-1) '_ARTIFACTS.mat'];
save(outputFile, 'movements','thres_zscore')


    function get_thres(hthres,~)
        thres=get(hthres,'Value');
    end
    function draw_now_function(hObject, ~, ~)
        draw_now=get(hObject,'Value');
    end
    function thresOK_function (hObject, ~, ~)
        thresOK=get(hObject,'Value');
    end
    function mov_ok_function(~,~,~)
        mov_ok=1;
    end
    function mov_no_ok_function(~,~,~)
        mov_ok=0;
    end
    function dir_back_function(~,~,~)
        dir=-1;
    end
    function dir_foward_function(~,~,~)
        dir=1;
    end
    function break_function(~,~,~)
        allIn=1;
    end

end

