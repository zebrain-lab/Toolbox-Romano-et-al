function [toDelete]=checkCells(fluoTraces,baseline,data, params)

%v = ver;
%parOn=any(strcmp('Parallel Computing Toolbox', {v.Name}));
parOn=0; 
%parOn=1; %UNCOMMENT THIS FOR PARALLELIZATION
if parOn
    numCores=feature('numcores');
    
    isOpen = matlabpool('size') > 0;
    if ~isOpen
        matlabpool('open','local', round(numCores*0.8))
    end
    
end


numCells=size(baseline,2);
numFrames=size(baseline,1);
baselineMeans=zeros(size(baseline));

if parOn
    parfor i=1:numCells
        
        baselineMeans(:,i)=runline(baseline(:,i),min(2000,round(numFrames/4)),1);
    end
else
    for i=1:numCells
        
        baselineMeans(:,i)=runline(baseline(:,i),min(2000,round(numFrames/4)),1);
    end
end
baselineVariations=zscore((baseline-baselineMeans)./baselineMeans);

% Check for ROIs with few pixels
thresh=params.cutOffDev;
deviations=baselineVariations<thresh;
toDeleteArtifacts=find(logical(sum(deviations,1)));
toDeleteArtifacts=unique([toDeleteArtifacts find(isnan(sum(fluoTraces,1)))]);
for i=1:numCells
    pixs(i)=length(data.cells{i});
end
toDeletePixs=find(pixs<params.cutOffPixels);


% Check for ROIs with weak baseline
toDeleteDim1=[];
for j=1:numCells
    dataSlice=fluoTraces(:,j)-baseline(:,j);
    
    [counts,x]=hist(dataSlice,min(round(length(dataSlice)/20),100));
    counts=smooth(counts);
    if any(isnan(x))
        toDeleteDim1=[toDeleteDim1 j];
        continue
    end

    
    
    [valueCenter,indCenter]=max(counts);
    
    if counts(1)>valueCenter*params.cutOffIntensity/100
        toDeleteDim1=[toDeleteDim1 j];
    end
    
end

toDeleteDim2=[];
for j=1:numCells
    if any(baseline(:,j)==0)
        toDeleteDim2=[toDeleteDim2 j];
    end
end

toDeleteDim=union(toDeleteDim1,toDeleteDim2);
toDelete=union(union(toDeleteDim,toDeleteArtifacts),toDeletePixs);
okCells=setdiff(1:numCells,toDelete);

% Plot results of sanity test
figure;imagesc(data.avg);
hold on; shading flat; colormap gray; axis image;

    lp(1)=line(data.cell_per{okCells(1)}(:,1),data.cell_per{okCells(1)}(:,2));
    set(lp(1),'color','y','linewidth',2);

if ~isempty (toDeleteDim1)
    for ind=toDeleteDim1(1)
        lp(2)=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
        set(lp(2),'color','b','linewidth',2);
        
    end
end
if ~isempty(toDeleteDim2)
    for ind=toDeleteDim2(1)
        lp(3)=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
        set(lp(3),'color','r','linewidth',2);
        
    end
end
if ~isempty(toDeleteArtifacts)
    for ind=toDeleteArtifacts(1)
        lp(4)=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
        set(lp(4),'color','g','linewidth',2);
        
    end
end
if ~isempty(toDeletePixs)
    for ind=toDeletePixs(1)
        lp(5)=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
        set(lp(5),'color','m','linewidth',2);
        
    end
end
ind_lp=find(lp~=0);
leg={'OK','Dim','baseline=0','Artifacts','Pixs'}';
legend(leg(ind_lp))

%%%%%%%%%%%%%

for i=2:length(okCells)
    ind=okCells(i);
    lp=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
    set(lp,'color','y','linewidth',2);
end


if length(toDeleteDim1)>=2
    for ind=toDeleteDim1(2:end)
        lp(1)=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
        set(lp(1),'color','b','linewidth',2);
        
    end
end
if length(toDeleteDim2)>=2
    for ind=toDeleteDim2(2:end)
        lp(2)=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
        set(lp(2),'color','r','linewidth',2);
        
    end
end
if length(toDeleteArtifacts)>=2
    for ind=toDeleteArtifacts(2:end)
        lp(3)=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
        set(lp(3),'color','g','linewidth',2);
        
    end
end
if length(toDeletePixs)>=2
    for ind=toDeletePixs(2:end)
        lp(4)=line(data.cell_per{ind}(:,1),data.cell_per{ind}(:,2));
        set(lp(4),'color','m','linewidth',2);
        
    end
end



