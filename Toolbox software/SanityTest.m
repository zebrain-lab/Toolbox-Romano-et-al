function [fluoTraces,F0,smoothBaseline,deletedCells]=SanityTest(filename, params)

%v = ver;
%parOn=any(strcmp('Parallel Computing Toolbox', {v.Name}));
parOn=0; 
%parOn=1; % UNCOMMENT FOR PARALLELIZATION
if parOn
    numCores=feature('numcores');
    
    
    isOpen = matlabpool('size') > 0;
    if ~isOpen
        matlabpool('open','local', round(numCores*0.8))
    end
    
end


data=load(filename);

if strcmp(params.neuropileSubtraction,'Yes')
    fluoTraces=data.cells_mean-params.alpha*data.npil_mean;
    fluoTraces=bsxfun(@minus,fluoTraces,min(fluoTraces));
else
    fluoTraces=data.cells_mean;
end

% Calculation the baseline in a running time window whose length is the maximum of 15 s or 40 time decays

twdw=max(15,40*params.tauDecay); 
wdw=round(params.fps*twdw);
numCells=data.cell_number;
numFrames=size(fluoTraces,1);
smoothBaseline=zeros(size(fluoTraces));

disp('1.1 Calculating fluorescence baseline of ROIs...')
if parOn
    parfor j=1:numCells
        dataSlice=fluoTraces(:,j);
        temp=zeros(numFrames-2*wdw,1);
        for i=wdw+1:numFrames-wdw
            temp(i-wdw)=prctile(dataSlice(i-wdw:i+wdw),8);
        end
        smoothBaseline(:,j)=[temp(1)*ones(wdw,1) ; temp; temp(end)*ones(wdw,1)];
        smoothBaseline(:,j)=runline(smoothBaseline(:,j),wdw,1);
        
    end
else
    for j=1:numCells
        dataSlice=fluoTraces(:,j);
        temp=zeros(numFrames-2*wdw,1);
        for i=wdw+1:numFrames-wdw
            temp(i-wdw)=prctile(dataSlice(i-wdw:i+wdw),8);
        end
        smoothBaseline(:,j)=[temp(1)*ones(wdw,1) ; temp; temp(end)*ones(wdw,1)];
        smoothBaseline(:,j)=runline(smoothBaseline(:,j),wdw,1);
    end
end

if strcmp(params.fluoBaselineCalculation.Method,'Smooth slow dynamics')
    params.fluoBaselineCalculation.timeWindow=wdw;
    F0=smoothBaseline;
else  strcmp(params.fluoBaselineCalculation.Method,'Average fluorescence on time window')  
    indStart=max(round(params.fps*params.fluoBaselineCalculation.t0),1);
    indEnd=min(round(params.fps*params.fluoBaselineCalculation.t1),numFrames);
    F0=repmat(mean(fluoTraces(indStart:indEnd,:)),size(fluoTraces,1),1);
end
% We look for cells that are not appropiate and remove them.


[deletedCells]=checkCells(fluoTraces,smoothBaseline,data, params);
toKeep=setdiff(1:numCells,deletedCells);

if ~isempty(deletedCells)
    avg=data.avg; bkg=data.bkg; bkg=data.bkg;
    cell_number=data.cell_number; cells_mean=data.cells_mean; npil_mean=data.npil_mean; cell_per=data.cell_per; cells=data.cells;
    distances=data.distances; pixelLengthX=data.pixelLengthX; pixelLengthY=data.pixelLengthY;
    
    cut=strfind(filename,'.mat');
    filenameBackUp=[filename(1:cut-1) '_Original.mat'];
    save(filenameBackUp,'avg','bkg','cell_number','cells_mean','npil_mean','cell_per','cells','distances','pixelLengthX','pixelLengthY');
    deleted.bkg=zeros(size(bkg));
    for i=1:length(deletedCells)
        bkg(data.cells{deletedCells(i)})=0;
        deleted.bkg(data.cells{deletedCells(i)})=1;
    end
    distances=data.distances(toKeep,toKeep);
    cell_number=length(toKeep); cells_mean=data.cells_mean(:,toKeep); npil_mean=data.npil_mean(:,toKeep); cell_per=data.cell_per(toKeep); cells=data.cells(toKeep);
    deleted.cell_number=length(deletedCells); deleted.cells_mean=data.cells_mean(:,deletedCells); deleted.cell_per=data.cell_per(deletedCells); deleted.cells=data.cells(deletedCells);
    save(filename,'avg','bkg','cell_number','cells_mean','npil_mean','cell_per','cells','distances','pixelLengthX','pixelLengthY','deleted');

end


fluoTraces=fluoTraces(:,toKeep);
smoothBaseline=smoothBaseline(:,toKeep);
F0=F0(:,toKeep);