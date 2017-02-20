function [densityData, densityNoise, xev, yev] = NoiseModel(filename, deltaFoF, sigma, movements, plotFlag)

lambda=8;        

deltaFoF(logical(movements),:)=NaN;
transfDataMatrix=bsxfun(@rdivide,deltaFoF,sigma);
points(:,1)=reshape(transfDataMatrix(1:end-1,:),[],1);
points(:,2)=reshape(transfDataMatrix(2:end,:),[],1);
points(isnan(points(:,1)),:)=[];
points(isnan(points(:,2)),:)=[];


pointsNeg=points(points(:,1)<0 & points(:,2)<0,1:2);
Sigma=cov(pointsNeg(:,1), pointsNeg(:,2));
Mu = [0 0];
Sigma=[1 Sigma(1,2); Sigma(2,1) 1];
dataGaussianCov = mvnrnd(Mu,Sigma,4*length(pointsNeg));
clear pointsNeg

mn=floor(min(min(min(transfDataMatrix)),min(reshape(dataGaussianCov,[],1))));
mx=ceil(max(max(max(transfDataMatrix)),max(reshape(dataGaussianCov,[],1))));
nevs=1000;
[xev yev]=meshgrid(linspace(mn,mx,nevs),linspace(mn,mx,nevs));

binColIdx=reshape(xev,[],1);
binRowIdx=repmat(linspace(mn,mx,nevs),1,1000.)';

[idx,dist]=knnsearch([points(:,1) points(:,2)],[binColIdx binRowIdx],'K',100);
binSpread=reshape(mean(dist,2),nevs,nevs);

clear idx dist

[idxNoise,distNoise]=knnsearch([dataGaussianCov(:,1) dataGaussianCov(:,2)],[binColIdx binRowIdx],'K',100);
binSpreadNoise=reshape(mean(distNoise,2),nevs,nevs);

clear binColIdx binRowIdx idxNoise distNoise

sizeWinFilt=8; % in sigmas
smoothParam= 1/(sizeWinFilt*min(min(binSpread)));
smoothParamNoise= 1/(sizeWinFilt*min(min(binSpreadNoise)));


pointCol=interp1(xev(1,:),1:nevs,points(:,1),'nearest');
pointRow=interp1(xev(1,:),1:nevs,points(:,2),'nearest');

hist2dim=accumarray([pointRow pointCol],1,[nevs nevs]);

pointColNoise=interp1(xev(1,:),1:nevs,dataGaussianCov(:,1) ,'nearest');
pointRowNoise=interp1(xev(1,:),1:nevs,dataGaussianCov(:,2) ,'nearest');

hist2dimNoise=accumarray([pointRowNoise pointColNoise],1,[nevs nevs]);



densityData=zeros(size(xev));
[binRows,binCols]=find(hist2dim);
for i=1:length(binRows)
    sigmaFilt=smoothParam*binSpread(binRows(i),binCols(i));
    if mod(ceil(sizeWinFilt*sigmaFilt),2) == 0
        binFilter=fspecial('gaussian',double([ceil(sizeWinFilt*sigmaFilt)+1 ceil(sizeWinFilt*sigmaFilt)+1]),double(sigmaFilt));
    else
        binFilter=fspecial('gaussian',double([ceil(sizeWinFilt*sigmaFilt) ceil(sizeWinFilt*sigmaFilt)]),double(sigmaFilt));
    end
   
    widthRect=(size(binFilter,1)-1)/2;
    centerRect=widthRect+1;
         
    rectRows=binRows(i)-min(binRows(i)-1,widthRect):binRows(i)+min(size(densityData,1)-binRows(i),widthRect);
    rectCols=binCols(i)-min(binCols(i)-1,widthRect):binCols(i)+min(size(densityData,2)-binCols(i),widthRect);
    rectFiltRows=centerRect-min(binRows(i)-1,widthRect):centerRect+min(size(densityData,1)-binRows(i),widthRect); 
    rectFiltCols=centerRect-min(binCols(i)-1,widthRect):centerRect+min(size(densityData,2)-binCols(i),widthRect);
    
    densityData(rectRows,rectCols)=densityData(rectRows,rectCols) + hist2dim(binRows(i),binCols(i))*binFilter(rectFiltRows,rectFiltCols);
    
end

densityNoise=zeros(size(xev));
[binRowsNoise,binColsNoise]=find(hist2dimNoise);
for i=1:length(binRowsNoise)
    sigmaFilt=smoothParamNoise*binSpreadNoise(binRowsNoise(i),binColsNoise(i));
    if mod(ceil(sizeWinFilt*sigmaFilt),2) == 0
        binFilter=fspecial('gaussian',double([ceil(sizeWinFilt*sigmaFilt)+1 ceil(sizeWinFilt*sigmaFilt)+1]),double(sigmaFilt));
    else
        binFilter=fspecial('gaussian',double([ceil(sizeWinFilt*sigmaFilt) ceil(sizeWinFilt*sigmaFilt)]),double(sigmaFilt));
    end
    
    widthRect=(size(binFilter,1)-1)/2;
    centerRect=widthRect+1;
   
    rectRows=binRowsNoise(i)-min(binRowsNoise(i)-1,widthRect):binRowsNoise(i)+min(size(densityNoise,1)-binRowsNoise(i),widthRect);
    rectCols=binColsNoise(i)-min(binColsNoise(i)-1,widthRect):binColsNoise(i)+min(size(densityNoise,2)-binColsNoise(i),widthRect);
    rectFiltRows=centerRect-min(binRowsNoise(i)-1,widthRect):centerRect+min(size(densityNoise,1)-binRowsNoise(i),widthRect);
    rectFiltCols=centerRect-min(binColsNoise(i)-1,widthRect):centerRect+min(size(densityNoise,2)-binColsNoise(i),widthRect);
    
    densityNoise(rectRows,rectCols)=densityNoise(rectRows,rectCols) + hist2dimNoise(binRowsNoise(i),binColsNoise(i))*binFilter(rectFiltRows,rectFiltCols);
    
end

densityData = Smooth1D(densityData,lambda);
densityData = Smooth1D(densityData',lambda)';
densityData = densityData./(sum(sum(densityData)));
densityNoise = Smooth1D(densityNoise,lambda);
densityNoise = Smooth1D(densityNoise',lambda)';
densityNoise = densityNoise./(sum(sum(densityNoise)));

if plotFlag
       
        figure
        z = log10(densityData); 
        contour(xev,yev,z,linspace(max(max(z))-4,max(max(z)),20)); 
        hold on
        z = log10(densityNoise); 
        contour(xev,yev,z,linspace(max(max(z))-4,max(max(z)),20)); 
        axis tight
        axis square 
        ext=get(gca,'XLim');
        xlim([ext(1)-1 ext(2)+1]); ylim([ext(1)-1 ext(2)+1])
        hC=colorbar;
        tcks=get(hC,'Ytick');
        set(hC,'Ytick',unique(round(tcks)),'YTicklabel',10.^unique(round(tcks)));
        hxlab=xlabel('Baseline noise normalized dFoF @ sample i'); hylab=ylabel('Baseline noise normalized dFoF @ sample i+1'); set(gcf,'color','w');
        set(gca,'FontSize',14); set([hxlab hylab],'FontSize',14)
        
        cut=strfind(filename,'_ALL_CELLS.mat');
        outfile=[filename(1:cut-1) '_densities.png'];
        export_fig(outfile);
end


function Z = Smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;   
