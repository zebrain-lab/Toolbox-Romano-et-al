% %  smooth the registration videos form template_matching (ImageJ toolbox)
% %  1) the user has to selecto the registered videos for smoothing
% %  2) if the data given by template_matching (Results) has the same name
% %     that the registered video the program automatically find it if not
% %     the user has to choose for each case...
% %  3) if the registered video has the same name than the original one  +
% %     _reg*.tif automatically find the original one, if not the user has to
% %     choose for each case...
% %  4) save the smooth video in a file called  (NAME_ORIGINAL)_reg_smooth.tif
% %
% %
% % each user should change the default pathname



clear all; close all
set(0,'DefaultAxesFontSize',14,'DefaultFigureColor','w', 'DefaultAxesTickDir', 'in','DefaultFigureWindowStyle','docked')



% choose files;
[filename_video, pathname_video] = uigetfile('*.tif','Open the original video file', 'MultiSelect', 'on');

if isequal(filename_video,0)
    disp('User selected Cancel')
elseif ischar(filename_video)  % if is only one video
    disp (['Video selected: ' fullfile(pathname_video, filename_video)])
    filename_data=fullfile(pathname_video,[filename_video(1:regexp(filename_video,'.tif')-1) '_reg.txt']); % search the data file
    if exist(filename_data,'file') % if exits open the data file
        data=importdata(filename_data);
    else % if not, ask to the user to choose the data file
        h=msgbox(['The name of the data file is not the same as the video file (' filename_video '), please choose the file'],'error','error','modal');
        uiwait(h)
        [filename_data_temp, pathname_data_temp] = uigetfile('*.txt',['choose the data file for video ' filename_video]...
            , 'MultiSelect', 'off',pathname_video);
        data=importdata(fullfile(pathname_data_temp, filename_data_temp));
        clear filename_data_temp pathname_data_temp h
    end
    % looks for the original video
    file_video_or=fullfile(pathname_video, filename_video);
    if ~ exist(file_video_or,'file')
        h=msgbox(['The original video is not found for video ( ' filename_video ' ), please choose the file'],'error','error','modal');
        uiwait(h)
        [file_video_or_temp, pathname_video_or_temp] = uigetfile('*.tif',['choose the data file for video ' filename_video]...
            , 'MultiSelect', 'off',pathname_video);
        file_video_or=fullfile(pathname_video_or_temp, file_video_or_temp);
        clear h pathname_video_or_temp file_video_or_temp
    end
    % convert to cell for compatibility with the case than more than one video is choosen
    filename_data=cellstr(filename_data);
    filename_video=cellstr(filename_video);
    file_video_or=cellstr(file_video_or);
    
else  % the same but for more than one video
    filename_data=cell(length(filename_video),1);
    file_video_or=cell(length(filename_video),1);
    for i=1:length(filename_video)
        disp(['video ' num2str(i) ': ' fullfile(pathname_video, filename_video{i})])
        filename_data_temp=fullfile(pathname_video,[filename_video{i}(1:regexp(filename_video{i},'.tif')-1) '_reg.txt']);
        if exist(filename_data_temp,'file')
            data_temp=importdata(filename_data_temp);
            filename_data{i}=filename_data_temp;
        else
            h=msgbox(['The name of the data file is not the same as the video file ( ' filename_video{i} ' ), please choose the file'],'error','error','modal');
            uiwait(h)
            [filename_data_temp, pathname_data_temp] = uigetfile('*.txt',['choose the data file for video ' filename_video{i}]...
                , 'MultiSelect', 'off',pathname_video);
            data_temp=importdata(fullfile(pathname_data_temp, filename_data_temp));
            filename_data{i}=fullfile(pathname_data_temp, filename_data_temp);
            clear filename_data_temp pathname_data_temp h
        end
        if exist('data','var')
            data=[data; data_temp];
            whos data
        else
            data=data_temp;
        end
        
        file_video_or_temp2=fullfile(pathname_video, filename_video{i});
        if exist(file_video_or_temp2,'file')
            file_video_or{i}=file_video_or_temp2;
        elseif ~exist(file_video_or_temp2,'file')
            h=msgbox(['The original video is not found for video ( ' filename_video ' ), please choose the file'],'error','error','modal');
            uiwait(h)
            [file_video_or_temp, pathname_video_or_temp] = uigetfile('*.tif',['choose the data file for video ' filename_video]...
                , 'MultiSelect', 'off',pathname_video);
            file_video_or{i}=fullfile(pathname_video_or_temp, file_video_or_temp);
            clear h pathname_video_or_temp file_video_or_temp
        end
    end
end
clear data_temp i file_video_or_temp2
%%
tic
for VIDEO=1:length(filename_video)
    a=tic;
    [~,ind_frames_ord]=sort(data(VIDEO).data(:,2)); % because the toolbox not necessarily starts in frame 1.
    % -1: because the first frame is the reference one, and has no deplacement.
    
    
    if data(VIDEO).data(1,2)==2
        X_des=[0; data(VIDEO).data(ind_frames_ord,3)];
        Y_des=[0; data(VIDEO).data(ind_frames_ord,4)];
    else
        X_des=[data(VIDEO).data(ind_frames_ord(1:data(VIDEO).data(1,2)),3); 0; data(VIDEO).data(ind_frames_ord(data(VIDEO).data(1,2)+1:end),3)];
        Y_des=[data(VIDEO).data(ind_frames_ord(1:data(VIDEO).data(1,2)),4); 0; data(VIDEO).data(ind_frames_ord(data(VIDEO).data(1,2)+1:end),4)];
        
    end
    
%     %for "jumpy" registrations:
    cut=1;

    xIn=zscore(diff(X_des))<cut & zscore(diff(X_des))>-cut;
    yIn=zscore(diff(Y_des))<cut & zscore(diff(Y_des))>-cut;

    [X_line]=runline(interp1(find(xIn),X_des(xIn),1:length(X_des),'pchip'),10,1);
    [Y_line]=runline(interp1(find(yIn),Y_des(yIn),1:length(Y_des),'pchip'),10,1);

    
    
    file_info=imfinfo(file_video_or{VIDEO});
    numFrames=length(file_info);
    
    
    FILE_VIDEO=file_video_or{VIDEO}
    for frames=1: numFrames
        xform=[1 0 0; 0 1 0; X_line(frames) Y_line(frames) 1];
        tform_translate = maketform('affine',xform);
        frame_temp=imread(FILE_VIDEO,'Index',frames);%,'PixelRegion', {[1 (info(1).Height-offset2)] [1 (info(1).Width-offset2)]} );
        
        [cb_trans xdata ydata]= imtransform(frame_temp, tform_translate, 'bicubic','Xdata',[1 file_info(1).Height],'Ydata', [1 file_info(1).Width],'size', size(frame_temp), 'fill', 0);
        
        %%%% save the smooth video in a file called: NAME_ORIGINAL_reg_smooth.tif
        if(frames==1)% if first image, erase previous file if it exists
            imwrite(cb_trans,[FILE_VIDEO(1:regexp(FILE_VIDEO,'.tif')-1) '_reg_smooth.tif'],'tif','WriteMode','overwrite');
        else
            imwrite(cb_trans,[FILE_VIDEO(1:regexp(FILE_VIDEO,'.tif')-1) '_reg_smooth.tif'],'tif','WriteMode','append');
        end
        
        
    end
    disp(['End video ' num2str(VIDEO) ', elapsed time ' num2str(toc(a))])
end

toc

