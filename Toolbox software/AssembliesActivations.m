[filename,pathname] = uigetfile({'*_CLUSTERS.mat';'*_CLUSTERS.MAT'},'Open file with assemblies data', 'MultiSelect', 'off');
filenameCLUSTER=fullfile(pathname,filename);

load(filenameCLUSTER);

prompt = {'Select the threshold p-value for the significance of an assembly activation'};
dlg_title = 'User input';
num_lines = 1;
def = {'0.01'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
threshSignifMatchIndex= str2num(answer{1});


numFrames=size(matchIndexTimeSeries,2);
totalClust=length(assembliesCells);

% This is to keep the complete transient around the peak
matchIndexTimeSeriesSignificant=zeros(totalClust,numFrames);
matchIndexTimeSeriesSignificantPeaks=zeros(totalClust,numFrames);

for indClust=1:totalClust
    tempSignif=(matchIndexTimeSeriesSignificance(indClust,:)<=threshSignifMatchIndex);
    matchIndexTimeSeriesSignificant(indClust,:)=zeros(size(matchIndexTimeSeries(indClust,:)));
    matchIndexTimeSeriesSignificantPeaks(indClust,tempSignif)=matchIndexTimeSeries(indClust,tempSignif);
    indsSignif=find(tempSignif);
    zerosMatch=find(matchIndexTimeSeries(indClust,:)==0);
    for i=1:length(indsSignif)
        temp=zerosMatch-indsSignif(i);
        indEnd=zerosMatch(find(temp>0,1,'first'));
        indStart=zerosMatch(find(temp<0,1,'last'));
        
        if isempty(indEnd)
            indEnd=numFrames;
        end
        if isempty(indStart)
            indStart=1;
        end
        
        matchIndexTimeSeriesSignificant(indClust,indStart:indEnd)=matchIndexTimeSeries(indClust,indStart:indEnd);
    end
end
if strcmp(clustering.method,'PCA-promax')
    save(filenameCLUSTER,'clustering','assembliesCells', 'assembliesVectors', 'PCsRot', 'confSynchBinary', 'matchIndexTimeSeries', 'matchIndexTimeSeriesSignificance','matchIndexTimeSeriesSignificant','matchIndexTimeSeriesSignificantPeaks','threshSignifMatchIndex')
else
    save(filenameCLUSTER,'clustering','assembliesCells', 'confSynchBinary', 'matchIndexTimeSeries', 'matchIndexTimeSeriesSignificance','matchIndexTimeSeriesSignificant','matchIndexTimeSeriesSignificantPeaks','threshSignifMatchIndex')
end
  

