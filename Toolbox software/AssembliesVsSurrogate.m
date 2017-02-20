function AssembliesVsSurrogate()



[filename,pathname] = uigetfile({'*_CLUSTERS.mat';'*_CLUSTERS.MAT'},'Open file with assemblies data', 'MultiSelect', 'off');
filenameCLUSTER=fullfile(pathname,filename);

clust=load(filenameCLUSTER);
assembliesCells=clust.assembliesCells;
cutName=strfind(filenameCLUSTER,'_CLUSTERS.mat');
dataAllCells=load([filenameCLUSTER(1:cutName-1) '_ALL_CELLS.mat']);

nms=load([filenameCLUSTER(1:cutName-1) '_SURROGATE_CLUSTERS.mat']);

[filename,pathname] = uigetfile({'*.mat';'*.MAT'},'Open file with variable to control against surrogate assemblies', 'MultiSelect', 'off');
filenameVARIABLE=fullfile(pathname,filename);

set(0,'DefaultFigureColor','w', 'DefaultAxesTickDir', 'out','DefaultFigureWindowStyle','docked')

vars=load(filenameVARIABLE);

repeat=1;

numCells=dataAllCells.cell_number;
numVars=length(vars.variable);
wrongInp=0;
for i=1:numVars
    if isfield(vars.variable{i},'x')
        sz=size(vars.variable{i}.y);
    else
        sz=size(vars.variable{i});
    end
    if sz(1)~=numCells 
        disp(['Incorrect format of input: Variable #' num2str(i) ' has a size of ' num2str(sz) '. First dimension (' num2str(sz(1)) ') should match the number or ROIs (' num2str(dataAllCells.cell_number) ')'])
        wrongInp=1;
        break
    end
end

if wrongInp==0
    for i=1:numVars
        sz=size(vars.variable{i});
        
         if sz(2)==1 & ~isfield(vars.variable{i},'x')
            if repeat
                varAssemblies{i}.data=vars.variable{i}([assembliesCells{:}]);
            else
                varAssemblies{i}.data=vars.variable{i}(unique([assembliesCells{:}]));
            end
            if isfield(nms,'assembliesSurrogateRandom')
                if repeat
                    varSurrogate{i}.Random.data=vars.variable{i}([nms.assembliesSurrogateRandom{:,:}]);
                else
                    varSurrogate{i}.Random.data=vars.variable{i}(unique([nms.assembliesSurrogateRandom{:,:}]));
                end
                compareVars(varAssemblies{i}.data,varSurrogate{i}.Random.data,i,'Random');
            end
            if isfield(nms,'assembliesSurrogateTopo')
                if repeat
                    varSurrogate{i}.Topo.data=vars.variable{i}([nms.assembliesSurrogateTopo{:,:}]);
                else
                    varSurrogate{i}.Topo.data=vars.variable{i}(unique([nms.assembliesSurrogateTopo{:,:}]));
                end
                compareVars(varAssemblies{i}.data,varSurrogate{i}.Topo.data,i,'Topographical');
            end
            
        end
        
        if sz(2)==dataAllCells.cell_number
            
            asym=any(reshape(vars.variable{i},[],1)~=reshape(vars.variable{i}',[],1));
            
            varAssemblies{i}.data=pickValsMatrix(assembliesCells,vars.variable{i},dataAllCells.cell_number,asym,repeat);
            if isfield(nms,'assembliesSurrogateRandom')
                varSurrogate{i}.Random.data=pickValsMatrix(nms.assembliesSurrogateRandom,vars.variable{i},dataAllCells.cell_number,asym,repeat);
                compareVars(varAssemblies{i}.data,varSurrogate{i}.Random.data,i,'Random');
            end
            if isfield(nms,'assembliesSurrogateTopo')
                varSurrogate{i}.Topo.data=pickValsMatrix(nms.assembliesSurrogateTopo,vars.variable{i},dataAllCells.cell_number,asym,repeat);
                compareVars(varAssemblies{i}.data,varSurrogate{i}.Topo.data,i,'Topographical');
            end
            
            
        end
        
        if isfield(vars.variable{i},'x')
            
            if repeat
                varAssemblies{i}.data=vars.variable{i}.y([assembliesCells{:}],:);
            else
                varAssemblies{i}.data=vars.variable{i}.y(unique([assembliesCells{:}]),:);
            end
            if isfield(nms,'assembliesSurrogateRandom')
                if repeat
                    varSurrogate{i}.Random.data=vars.variable{i}.y([nms.assembliesSurrogateRandom{:,:}],:);
                else
                    varSurrogate{i}.Random.data=vars.variable{i}.y(unique([nms.aassembliesSurrogateRandom{:,:}]),:);
                end
                compareVars2(vars.variable{i}.x,varAssemblies{i}.data,varSurrogate{i}.Random.data,i,'Random');
            end
            if isfield(nms,'assembliesSurrogateTopo')
                if repeat
                    varSurrogate{i}.Topo.data=vars.variable{i}.y([nms.assembliesSurrogateTopo{:,:}],:);
                else
                    varSurrogate{i}.Topo.data=vars.variable{i}.y(unique([nms.assembliesSurrogateTopo{:,:}]),:);
                end
                compareVars2(vars.variable{i}.x,varAssemblies{i}.data,varSurrogate{i}.Topo.data,i,'Topographical');
            end
            
        end
        
    end
end

save([filenameCLUSTER(1:cutName-1)  '_ASSEMBLIES_vs_SURROGATE.mat'],'varAssemblies','varSurrogate')

    function compareVars(dataAss,dataNM,varNum,kind)
        if strcmp(kind,'Random')
            figure('Name',['Variable #' num2str(varNum) '; Assemblies vs. Random Surrogate Data'])
        elseif  strcmp(kind,'Topographical')
            figure('Name',['Variable #' num2str(varNum) '; Assemblies vs. Topographical Surrogate Data'])
        end
        mn1=min(dataAss); mx1=max(dataAss);
        mn2=min(dataNM); mx2=max(dataNM);
        x=linspace(min(mn1,mn2),max(mx1,mx2),40);
        subplot(2,1,1)
        [count]=hist(dataAss,x);
        bar(x,count/sum(count),1,'k');
        xlabel(['Variable #' num2str(varNum)]); ylabel('Frequency')
        title('Assemblies')
        subplot(2,1,2)
        [count2]=hist(dataNM,x);
        bar(x,count2/sum(count2),1,'k');
        xlabel(['Variable #' num2str(varNum)]); ylabel('Frequency')
        title(kind)
        
    end

    function compareVars2(x,dataAss,dataNM,varNum,kind)
        if strcmp(kind,'Random')
            figure('Name',['Variable #' num2str(varNum) '; Assemblies vs. Random Surrogate Data'])
        elseif  strcmp(kind,'Topographical')
            figure('Name',['Variable #' num2str(varNum) '; Assemblies vs. Topographical Surrogate Data'])
        end
        
        subplot(2,1,1)
        plot(x,dataAss,'Color',[.8 .8 .8])
        hold on;
        plot(x,mean(dataAss),'k','LineWidth',2)
        xlabel('x'); ylabel('y')
        title('Assemblies')
        subplot(2,1,2)
        plot(x,dataNM,'Color',[.8 .8 .8])
        hold on;
        plot(x,mean(dataNM),'k','LineWidth',2)
        xlabel('x'); ylabel('Variable')
        title(kind)
        
    end

    function valsPicked=pickValsMatrix(list,Mat,numCells,asym,repeat)
        if size(list,1)>1
            list=reshape(list,1,[]);
        end
        
        pairs=zeros(numCells,numCells);
        if ~asym
            for k=1:length(list)
                temp=list{k};
                for n=1:length(temp)
                    for m=n+1:length(temp)
                        pairs(temp(n),temp(m))=pairs(temp(n),temp(m))+1;
                    end
                end
            end
        else
            for k=1:length(list)
                temp=list{k};
                for n=1:length(temp)
                    for m=1:length(temp)
                        if m==n
                            continue
                        else
                            pairs(temp(n),temp(m))=pairs(temp(n),temp(m))+1;
                        end
                    end
                end
            end
        end
        pairs=reshape(pairs,[],1);
        Mat=reshape(Mat,[],1);
        inds=find(pairs);
        
        if repeat
            numRepeats=pairs(inds);
        else
            numRepeats=ones(length(inds),1);
        end
        
        valsPicked=zeros(sum(numRepeats),1);
        count=1;
        for k=1:length(inds)
            for kk=1:numRepeats(k)
                valsPicked(count)=Mat(inds(k));
                count=count+1;
            end
        end
        
        
    end



end