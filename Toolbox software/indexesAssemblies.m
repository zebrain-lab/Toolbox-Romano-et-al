function indexesAssemblies()

[filename,pathname] = uigetfile({'*.mat';'*.MAT'},'Open file with raster data', 'MultiSelect', 'off');
filenameRASTER=fullfile(pathname,filename);


[filename,pathname] = uigetfile({'*.mat';'*.MAT'},'Open file with assemblies data', 'MultiSelect', 'off');
filenameCLUSTER=fullfile(pathname,filename);


temp=load(filenameRASTER);
raster=temp.raster; movements=temp.movements; deltaFoF=temp.deltaFoF;

myClusters=load(filenameCLUSTER);
temp=raster;
temp(logical(movements),:)=0;
rasterAnalog=deltaFoF;
rasterAnalog(~temp)=0;


data=zscore(rasterAnalog);

dists = pdist(data','euclidean');
distMat=squareform(dists);


[PCs,Score,eigenvals]=princomp(data);
maxEigenValPastur=(1+sqrt(size(data,2)/size(data,1)))^2;
minEigenValPastur=(1-sqrt(size(data,2)/size(data,1)))^2;
correctionTracyWidom=size(data,2)^(-2/3);
smaller = eigenvals < maxEigenValPastur + correctionTracyWidom;
cutOffPC = find(smaller,1)-1; % The significant PCs go up to cutOffPC
dataPC=PCs(:,1:cutOffPC)';


distsPCs = pdist(dataPC','euclidean');
distMatPC=squareform(distsPCs);


hubertValue=hubert(distMat,myClusters);
hubertValuePC=hubert(distMatPC,myClusters);
DBValue=dbIndex(data,myClusters);
DBValuePC=dbIndex(dataPC,myClusters);

disp(['The Normalized Hubert Gamma index for the clusters is (the higher the better): ' num2str(hubertValue) ' in data space; ' num2str(hubertValuePC) ' in reduced-dimensionality space.']);
disp(['The Davies-Bouldin index for the clusters is (the lower the better): ' num2str(DBValue) ' in data space; ' num2str(DBValuePC) ' in reduced-dimensionality space.']);

    function hubertValue=hubert(distMat,myClusters)
        Q=ones(size(distMat));
        for n=1:length(myClusters.assembliesCells)
            for i=1:length(myClusters.assembliesCells{n})
                cellsIn=myClusters.assembliesCells{n};
                for j=1:length(cellsIn)
                    Q(cellsIn(j),cellsIn)=0;
                end
                
            end
        end
        
        muP=mean(distMat(logical(triu(ones(size(distMat)),1))));
        muQ=mean(Q(logical(triu(ones(size(distMat)),1))));
        stdP=std(distMat(logical(triu(ones(size(distMat)),1))));
        stdQ=std(Q(logical(triu(ones(size(distMat)),1))));
        multiplied=((distMat-muP).*(Q-muQ));
        
        sumMultiplied=0;
        for i=1:size(distMat,1)-1
            sumMultiplied=sumMultiplied+sum(multiplied(i,i+1:size(distMat,1)));
            
        end
        
        pairs=size(distMat,1)*(size(distMat,1)-1)/2;
        hubertValue=sumMultiplied/(stdP*stdQ*pairs);
    end

    function DBValue=dbIndex(data,myClusters)
        sizeTemp=length([myClusters.assembliesCells{:}]);
        dataTemp=nan(sizeTemp,size(data,1)); clust=nan(sizeTemp,1);
        
        cellsEns=unique(reshape([myClusters.assembliesCells{:}],[],1));
        count=1;
        for n=1:length(myClusters.assembliesCells)
            for j=1:length(myClusters.assembliesCells{n})
                targetCell=myClusters.assembliesCells{n}(j);
                dataTemp(count,:)=data(:,targetCell)';
                clust(count)=n;
                count=count+1;
                
            end
        end
        %Include cells not in assemblies
        cellsNotInEns=setdiff(1:size(myClusters.PCsRot,1),cellsEns);
        dataTemp=[dataTemp; data(:,cellsNotInEns)'];
        clust=[clust; (n+1)*ones(length(cellsNotInEns),1)];
        [DBValue,~] = db_index(dataTemp, clust);
    end
end