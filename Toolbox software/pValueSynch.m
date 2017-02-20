function pValue=pValueSynch(raster,repetitions)

    rasterShuffle=raster;
    x=0:1:size(raster,2);
    histTotal=zeros(1,size(x,2));
    for j=1:repetitions
        for i=1:size(raster,2)
            % This is for permutations of spike times, preserving only spike count
            % threshShuffle(:,i)=datos(randperm(size(datos,1)),i);

            % This is for shuffling ISIs, preserving therefore both single cell spike count
            % and ISI
            spikeTrain=raster(:,i);
            indSpike=find(spikeTrain);
            ISIs=diff(indSpike);

            spikeTrainShuffled=zeros(size(spikeTrain));
            spikeTrainShuffled(1)=spikeTrain(1);
            spikeTrainShuffled(cumsum(ISIs(randperm(size(ISIs,1))))+1)=1;
            rasterShuffle(:,i)=spikeTrainShuffled;
        end

        rate=sum(rasterShuffle(2:end,:),2);
        [histRate, x]=hist(rate,x);
        histRate=histRate./sum(histRate);
       
        histTotal=histTotal + histRate./repetitions;

    end

    pValue=1-cumsum(histTotal);
    pValue(find(pValue<0))=pValue(find(pValue<0,1,'First')-1);



