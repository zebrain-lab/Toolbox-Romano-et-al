function [mapOfOdds] = SignificantOdds(deltaFoF, sigma, movements, densityData, densityNoise, xev, params, plotFlag)
    
    deltaFoF(logical(movements),:)=NaN;
    transfDataMatrix=bsxfun(@rdivide,deltaFoF,sigma);
    points(:,1)=reshape(transfDataMatrix(1:end-1,:),[],1);
    points(:,2)=reshape(transfDataMatrix(2:end,:),[],1);
    points(isnan(points(:,1)),:)=[];
    points(isnan(points(:,2)),:)=[];
    
    
    pCutOff=(100-params.confCutOff)/100;
    mapOfOdds=densityNoise<=pCutOff*densityData;
    
    if plotFlag
        % We plot the result
        figure; plot(points(:,1),points(:,2),'k.'); axis equal; hold on
        h=imagesc(xev(1,:),xev(1,:), mapOfOdds); axis xy; axis tight
        alpha(h,0.6);
        hxlab=xlabel('z-transformed value @ sample i'); hylab=ylabel('z-transformed value @ sample i+1'); set(gcf,'color','w');
        set(gca,'FontSize',14); set([hxlab hylab],'FontSize',14)
        
    end
end