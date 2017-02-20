function fileName=getVideoFile()


default_pathname='';
% choose files;
[filename_video, pathname_video] = uigetfile('*.tif','Open the video file', 'MultiSelect', 'off',default_pathname);

if isequal(filename_video,0)
    disp('User selected Cancel')
else ischar(filename_video)  % if is only one video
    disp (['Video selected: ' fullfile(pathname_video, filename_video)])
    fileName=fullfile(pathname_video,filename_video);
    
end
