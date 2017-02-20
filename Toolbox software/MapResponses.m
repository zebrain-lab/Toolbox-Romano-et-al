function [roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,mapData,roiHSVRescaled,peakResponseLabels,peakParameterStdLabels,cellper,pixROIs,mx,cutUp,cutDown,offset,traces,onsets,offsets,responses,responsesStd,prevFrame,eventDuration,params,plotRows,transpa,mapParams,fnameOut,hIm]=MapResponses();
[filename,pathname] = uigetfile({'*_RASTER.mat';'*_RASTER.MAT'},'Open file with Raster data', 'MultiSelect', 'off');
filenameRASTER=fullfile(pathname,filename);
load(filenameRASTER);

if ~any(any(raster==0))
   raster=zeros(size(raster)); 
end

[filename,pathname] = uigetfile({'*.mat';'*.MAT'},'Open file with timing data', 'MultiSelect', 'off');
load(fullfile(pathname,filename));

cutName=strfind(filenameRASTER,'_RASTER.mat');
filenameALL_CELLS=[filenameRASTER(1:cutName-1) '_ALL_CELLS.mat'];
if exist(filenameALL_CELLS, 'file') == 2
    dataAllCells=load(filenameALL_CELLS);
end

fnameOut=[filenameRASTER(1:cutName-1) '_RESPONSE_MAP.mat'];

prompt = {'Event duration (s)','Duration of time window after event for ROI response calculation (s)'};
dlg_title = 'Parameters for color remapping';
num_lines = 1;
def = {num2str(1/params.fps),'3'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
eventDuration=round(str2num(answer{1})*params.fps);
postFrame= round(str2num(answer{2})*params.fps);
prevFrame=round(1.5*params.fps);
    
wdwFrame=postFrame;


numFrames=size(deltaFoF,1);
mapParams=nan(length(mapData),1);
mapFrames=cell(length(mapData),1);
mapLabels=cell(length(mapData),1);
numStim=length(mapData);
for i=1:numStim
    mapParams(i)=mapData{i}.value;
    mapFrames{i}=round(mapData{i}.onsetTime*params.fps);
    mapLabels{i}=str2num(mapData{i}.label);
end
[junk,indSort]=sort(mapParams);
mapParams= mapParams(indSort);
for i=1:numStim
    mapFrames{i}= mapFrames{indSort(i)};
end

x=((mapParams-mapParams(1)))/(mapParams(end)-mapParams(1));

totalCells=size(deltaFoF,2);

plotRows=1;
transpa=1;
%%
peakParameter=zeros(totalCells,1);
peakResponse=zeros(totalCells,1);
peakParameterStd=zeros(totalCells,1);
responses=nan(numStim,totalCells);
responsesStd=nan(numStim,totalCells);
roiHSV=zeros(totalCells,3);
traces=cell(totalCells,numStim);
colTraces=cell(totalCells,numStim);
onsets=cell(totalCells,numStim);
offsets=cell(totalCells,numStim);
for numCell=1:totalCells
   
    for i=1:numStim
        
        numTrials=length(mapFrames{i});
        responsesTrials=nan(numTrials,1);
        
        traces{numCell,i}=nan(numTrials,prevFrame+postFrame+1);
        colTraces{numCell,i}=nan(numTrials,prevFrame+postFrame+1);
        
        for j=1:numTrials
            tempTrace=nan(prevFrame+postFrame+1,1);
            indStart=max([mapFrames{i}(j)-prevFrame 1]);
            indEnd=min([mapFrames{i}(j)+postFrame numFrames]);
            
            indStartTemp=prevFrame-(mapFrames{i}(j)-indStart)+1;
            tempTrace(indStartTemp:indStartTemp+length(indStart:indEnd)-1)=deltaFoF(indStart:indEnd,numCell);
            
            
            traces{numCell,i}(j,:)=tempTrace;
            responsesTrials(j)=nanmean(traces{numCell,i}(j,prevFrame+1:prevFrame+1+wdwFrame));
            
            tempTrace=zeros(prevFrame+postFrame+1,1);
            tempTrace(indStartTemp:indStartTemp+length(indStart:indEnd)-1)=raster(indStart:indEnd,numCell);
            colTraces{numCell,i}(j,:)=tempTrace;
            
            temp=[0 colTraces{numCell,i}(j,:) 0];
            onsets{numCell,i}{j}=find(diff(temp)==1);
            offsets{numCell,i}{j}=find(diff(temp)==-1)-1;
             
        end
        responses(i,numCell)=nanmean(responsesTrials);
        responsesStd(i,numCell)=nanstd(responsesTrials);
    end
    responsesCorrected=responses(:,numCell);
    responsesCorrected(responsesCorrected<0)=0;
    
    [peakResponse(numCell),indPos]=max(responsesCorrected);
    peakParameter(numCell)=mapParams(indPos);
    if numStim==1
        hueInd=1;
        satInd=1;
        peakParameterStd(numCell)=0;
    else
        if length(indPos:numStim)==length(1:indPos)
            newX=x';
            newResponse=responsesCorrected;
            newInds=1:numStim;
            newMapParams=mapParams';
        else
            [junk,opt]=max([length(indPos:numStim) length(1:indPos)]);
            if opt==1
                newInds=[length(x):-1:2*indPos 1:length(x)];
                tempX=x(end:-1:2*indPos);
                newX=[(2*x(indPos)-tempX)'  x'];
                newX=newX/max(newX);
                newResponse=responsesCorrected(newInds);
                
                tempX=mapParams(end:-1:2*indPos);
                newMapParams=[(2*mapParams(indPos)-tempX)'  mapParams'];
                
            else
                newInds=[1:length(x) length(x)-(2*(length(x)-indPos)+1):-1:1];
                tempX=x(length(x)-(2*(length(x)-indPos)+1):-1:1);
                newX=[x'  (2*x(indPos)-tempX)'];
                newX=newX/max(newX);
                newResponse=responsesCorrected(newInds);
                
                tempX=mapParams(length(mapParams)-(2*(length(mapParams)-indPos)+1):-1:1);
                newMapParams=[mapParams'  (2*mapParams(indPos)-tempX)'];
                
            end
        end
        
        probNorm=zeros(1,length(newX));
        probNorm(1)=0.5; probNorm(end)=0.5;
        [valMax,indMax]=max(probNorm);
        maxStd=sqrt(sum(probNorm.*((newX-indMax).^2)));
        
        prob=newResponse/sum(newResponse);
        [valMax,indMax]=max(newResponse);
        xMax=newX(indMax);
        paramMax=newMapParams(indMax);
        valStdPeak=sqrt(sum(prob'.*((newX-xMax).^2)));
        if ~any(newResponse)
            valStdPeakNorm=0;
            peakParameterStd(numCell)=NaN;
        else
            valStdPeakNorm=valStdPeak/maxStd;
            peakParameterStd(numCell)=(mapParams(end)-mapParams(1))*sqrt(sum(prob'.*((newMapParams-paramMax).^2)))/(newMapParams(end)-newMapParams(1));
        end
        [valMax,indMax]=max(responsesCorrected);
        indMax=x(indMax);
        hueInd=indMax;
        satInd=1-valStdPeakNorm;
    end
    roiHSV(numCell,:)=[hueInd satInd peakResponse(numCell)];
    
end

%%
remap=0;
cutUp{2}=1;
cutDown{2}=0;
offset{2}=0;
cutUp{3}=1;
cutDown{3}=0;
offset{3}=0;
im=im2single(dataAllCells.avg);
cellper=dataAllCells.cell_per;
pixROIs=dataAllCells.cells;
[roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,mapData,roiHSVRescaled,peakResponseLabels,peakParameterStdLabels,cellper,mx,cutUp,cutDown,offset,transpa,hIm]=remapResponseColors(roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,remap,cellper,cutUp,cutDown,offset,mapData,pixROIs,transpa);
