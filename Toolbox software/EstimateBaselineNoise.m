function [deltaFoF, mu, sigma, params]=EstimateBaselineNoise(filename, deltaFoF, params)

ansMethod = questdlg('Choose method for estimation of noise in baseline fluorescence', 'Select method', 'Gaussian model', 'Standard deviation', 'Gaussian model');
params.BaselineNoiseMethod=ansMethod;

data=load(filename);
if strcmp(ansMethod,'Gaussian model')
    
    
    numCells=size(deltaFoF,2);
    numFrames=size(deltaFoF,1);
    
    
    % We calculate the ROI's baseline noise by fitting a
    % gaussian to the distribution (local density estimation) of negative values of deltaFoF
    
    for numcell=1:numCells
        
        dataCell=deltaFoF(:,numcell);
          
        [smoothDist,x] = ksdensity(dataCell);
        
        [valuePeak,indPeak]=max(smoothDist);
        
        xFit=x(1:indPeak);
        dataToFit=smoothDist(1:indPeak)/numFrames;
        [sigma(numcell),mu(numcell),A]=mygaussfit(xFit',dataToFit);
        
        if ~isreal(sigma(numcell))
            dev=nanstd(dataCell);
            outliers=abs(deltaFoF)>2*dev;
            
            deltaF2=dataCell;
            deltaF2(outliers)=NaN;
            sigma(numcell)=nanstd(deltaF2);
            mu(numcell)=nanmean(deltaF2);
            
        end
        
        
        distFit=A*exp(-(x-mu(numcell)).^2./(2*sigma(numcell)^2));
        
    end
    
    deltaFoF=bsxfun(@minus,deltaFoF, mu);
    
elseif  strcmp(ansMethod,'Standard deviation')
    
    dev=nanstd(deltaFoF);
    outliers=bsxfun(@gt,abs(deltaFoF),2*dev);
    % We calculate the standard deviation of each cell without outliers
    deltaF2=deltaFoF;
    deltaF2(outliers)=NaN;
    sigma=nanstd(deltaF2);
    mu=nanmean(deltaF2);
    deltaFoF=bsxfun(@minus,deltaFoF, mu);
    
    
end
end



