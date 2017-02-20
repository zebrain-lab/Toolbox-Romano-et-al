function [raster, mapOfOddsJoint]=Rasterize(deltaFoF, sigma, movements, mapOfOdds, xev, yev, params)

temp=deltaFoF;
temp(logical(movements),:)=0;
transfDataMatrix=bsxfun(@rdivide,temp,sigma);

cc=bwconncomp(~mapOfOdds);
stats=regionprops(cc,'PixelList');
indZero=find(xev(1,:)>0,1,'first');% 


for i=1:cc.NumObjects
    if ismember(indZero,stats(i).PixelList(:,1))
        break
    end
end
mapOfOddsCorrected=ones(size(mapOfOdds));
mapOfOddsCorrected(cc.PixelIdxList{i})=0;

%%
noiseBias=1.5;
factorDecay=exp(-1/(params.fps*params.tauDecay));
decayMap=ones(size(mapOfOdds));
rowsDecay=1:size(decayMap,1);
for i=1:size(decayMap,2)
    decayMap(rowsDecay(yev(:,1)<factorDecay*(xev(1,i)-noiseBias)-noiseBias),i)=0;
end

riseMap=ones(size(mapOfOdds));
riseMap(end,:)=0; 

mapOfOddsJoint = mapOfOddsCorrected & riseMap & decayMap;

raster=zeros(size(transfDataMatrix,1),size(transfDataMatrix,2));
for numNeuron=1:size(transfDataMatrix,2)
    [junk,bins] = histc(transfDataMatrix(1:end,numNeuron),xev(1,:));
    for numFrame=3:size(transfDataMatrix,1)-2
        optA= (mapOfOddsJoint(bins(numFrame+1),bins(numFrame)) & mapOfOddsJoint(bins(numFrame),bins(numFrame-1)));
        optB= (mapOfOddsJoint(bins(numFrame),bins(numFrame-1)) & mapOfOddsJoint(bins(numFrame-1),bins(numFrame-2)));
        optC= (mapOfOddsJoint(bins(numFrame+2),bins(numFrame+1)) & mapOfOddsJoint(bins(numFrame+1),bins(numFrame)));
        if optA | optB | optC
            raster(numFrame,numNeuron)=1;
          
        end
        
        
    end
end



