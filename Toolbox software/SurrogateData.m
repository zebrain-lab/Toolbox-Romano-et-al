function SurrogateData()
[filename,pathname] = uigetfile({'*_CLUSTERS.mat';'*_CLUSTERS.MAT'},'Open file with assemblies data', 'MultiSelect', 'off');
filenameCLUSTER=fullfile(pathname,filename);

dataClust=load(filenameCLUSTER);
assembliesCells=dataClust.assembliesCells;
cutName=strfind(filenameCLUSTER,'_CLUSTERS.mat');
dataAllCells=load([filenameCLUSTER(1:cutName-1) '_ALL_CELLS.mat']);

whichNull = questdlg('What kind of surrogate data do you want to generate?', 'Surrogate data', 'Random surrogate assemblies', 'Topographical surrogate assemblies', 'Both','Random surrogate data');

prompt = {['Number of surrogate assemblies per real assembly']};
def={'100'};
answer = inputdlg(prompt, 'Indicate the number of syntethic surrogate assemblies per assembly', 1,def);
numRepet=str2num(answer{1});

numAssemblies=length(assembliesCells);
numNeurons=dataAllCells.cell_number;

if strcmp(whichNull,'Random surrogate assemblies')
    assembliesSurrogateRandom=generateRandom(numRepet);
    save([filenameCLUSTER(1:cutName-1) '_SURROGATE_CLUSTERS.mat'],'assembliesSurrogateRandom')
    disp('Done');
elseif strcmp(whichNull,'Topographical surrogate assemblies')
    assembliesSurrogateTopo=generateTopo(numRepet);
    save([filenameCLUSTER(1:cutName-1) '_SURROGATE_CLUSTERS.mat'],'assembliesSurrogateTopo')
    disp('Done');
else
    assembliesSurrogateRandom=generateRandom(numRepet);
    assembliesSurrogateTopo=generateTopo(numRepet);
    save([filenameCLUSTER(1:cutName-1) '_SURROGATE_CLUSTERS.mat'],'assembliesSurrogateRandom','assembliesSurrogateTopo')
    disp('Done');
end

    function assembliesSurrogateRandom=generateRandom(numRepet,~)
        assembliesSurrogateRandom=cell(numAssemblies,numRepet);
  
        for j=1:numRepet;
            permCells=randperm(numNeurons);
            for i=1:numAssemblies
               assembliesSurrogateRandom{i,j}=permCells(assembliesCells{i});
          
            end
        end
    end

    function assembliesSurrogateTopo=generateTopo(numRepet,~)
        
        distancesOriginal=dataAllCells.distances;
        pValueShuff=nan(length(assembliesCells),numRepet);
        assembliesSurrogateTopo=cell(length(assembliesCells),numRepet);
        for numShuff=1:numRepet
            for k=1:length(assembliesCells)
                
                cellsAssembly=assembliesCells{k};
                
                temp=distancesOriginal(cellsAssembly,cellsAssembly);
                distsAssembly=temp(logical(triu(ones(size(temp)),1)));
                countTry=0;
                prob=0;
                while prob<0.05 & countTry<20
                    distsShuff=nan(size(distsAssembly));
                    shuffledAssembly=nan(size(cellsAssembly));
                    cellsToPick=setdiff(1:numNeurons,cellsAssembly);
                    distsToCopy=distsAssembly;
                    for i=1:length(cellsAssembly)
                        if i==1
                            targetCell=cellsToPick(randi(length(cellsToPick)));
                            shuffledAssembly(i)=targetCell;
                            cellsToPick=setdiff(cellsToPick,targetCell);
                        else
                            distsToPick=distancesOriginal(shuffledAssembly(~isnan(shuffledAssembly)),cellsToPick);
                            squaredErrorTemp=nan(size(distsToPick));
                            copiedFrom=nan(size(distsToPick));
                            for j=1:size(distsToPick,2) 
                                distsToCopyTemp=distsToCopy;
                                for n=1:size(distsToPick,1) 
                                    [minDist ind]=nanmin((distsToCopyTemp-distsToPick(n,j)).^2);
                                    squaredErrorTemp(n,j)=minDist;
                                    copiedFrom(n,j)=ind; 
                                    distsToCopyTemp(ind)=nan;
                                end
                            end
                            [junk, indMins]=min(mean(squaredErrorTemp,1));
                            
                            targetCell=cellsToPick(indMins);
                            shuffledAssembly(i)=targetCell;
                            cellsToPick=setdiff(cellsToPick,targetCell);
                            squaredErrors((i-1)*((i-1)-1)/2+1:i*(i-1)/2)=squaredErrorTemp(:,indMins);
                            distsShuff((i-1)*((i-1)-1)/2+1:i*(i-1)/2)=distsToPick(:,indMins);
                            distsToCopy(copiedFrom(:,indMins))=[];
                            
                        end
                        
                    end
                    
                    if length(distsAssembly)>2
                        [h,prob]=kstest2(distsAssembly,distsShuff);
                    else
                        prob=1;
                    end
                    
                    if prob>0.05
                        pValueShuff(k,numShuff)=prob;
                        assembliesSurrogateTopo{k,numShuff}=shuffledAssembly;
                    end
                    countTry=countTry+1;
                    
                    
                end
                disp(['#shuffle ' num2str(numShuff) ' of assembly ' num2str(k)])
            end
        end
    end


end