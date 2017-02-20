function AssembliesDynamicsVisualization();
[filename,pathname] = uigetfile({'*_CLUSTERS.mat';'*_CLUSTERS.MAT'},'Open file with assemblies data', 'MultiSelect', 'off');
filenameCLUSTER=fullfile(pathname,filename);

temp=load(filenameCLUSTER);
assembliesCells=temp.assembliesCells;
matchIndexTimeSeriesSignificant=temp.matchIndexTimeSeriesSignificant;
confSynchBinary=temp.confSynchBinary;

cutName=strfind(filenameCLUSTER,'_CLUSTERS.mat');
dataAllCells=load([filenameCLUSTER(1:cutName-1) '_ALL_CELLS.mat']);
temp=load([filenameCLUSTER(1:cutName-1) '_RASTER.mat']);
deltaFoF=temp.deltaFoF;
raster=temp.raster;
movements=temp.movements;

rasterAnalog=deltaFoF;
rasterClean=double(raster);
rasterClean(logical(movements),:)=0;
rasterAnalog(logical(movements),:)=0;
rasterAnalog(~logical(rasterClean))=0;

numFrames=size(rasterClean,1);

[filename,pathname] = uigetfile({'*_ORDER_TOPO.mat';'*_ORDER_TOPO.MAT'},'Open file with assemblies order', 'MultiSelect', 'off');
orderData=load(fullfile(pathname,filename));


areThereTags = questdlg('Do you want to add time tags?', 'Time tags', 'Yes', 'No', 'Yes');
if strcmp(areThereTags,'Yes')
    [filename,pathname] = uigetfile({'*.mat';'*.MAT'},'Open file with experimental time tags', 'MultiSelect', 'off');
    tagsData=load(fullfile(pathname,filename));
    
    tagColors=jet(max(tagsData.timeTags));
    
    for j=1:max(tagsData.timeTags)
        tempTags=zeros(size(tagsData.timeTags,1),1);
        tempTags(tagsData.timeTags==j)=1;
        indsSignif=find(tempTags);
        zerosMatch=find(tempTags==0);
        for i=1:length(indsSignif)
            temp=zerosMatch-indsSignif(i);
            indEnd=zerosMatch(find(temp>0,1,'first'))-1;
            indStart=zerosMatch(find(temp<0,1,'last'))+1;
            if isempty(indEnd)
                indEnd=numFrames;
            end
            if isempty(indStart)
                indStart=1;
            end
            if indEnd==indStart
                indEnd=indEnd+1;
            end
            tagsData.tagStart{j}(i)=indStart;
            tagsData.tagEnd{j}(i)=indEnd;
            
        end
        tagsData.tagStart{j}=unique(tagsData.tagStart{j});
        tagsData.tagEnd{j}=unique(tagsData.tagEnd{j});
    end
end

%%

scrsz=get(0,'ScreenSize');
figure;
set(gcf,'Position',scrsz,'Color','w')
h1=subplot(4,1,[1 2]);
imagesc(rasterAnalog(:,orderData.orderOfCells)')
colormap(1-gray)
temp=reshape(rasterAnalog,[],1); mxDFoF=prctile(temp(temp>0),95);
hold on; plot([1 numFrames],[length(orderData.orderOfCellsInEnsembles) length(orderData.orderOfCellsInEnsembles)],'r-','LineWidth',2);
caxis([0 mxDFoF])
pos=get(h1,'Position');
hBar=colorbar;
posBar=get(hBar,'Position');
set(h1,'XTickLabel',[],'LineWidth',2);
set(h1,'TickLength',[0.005 0.025],'TickDir','out')
set(hBar,'Position',[pos(1)+pos(3)+0.01 pos(2) posBar(3)/3 pos(4)])
set(hBar,'YTick',[0 mxDFoF/2 mxDFoF], 'YTickLabel',{'0', num2str(sprintf('%6.1f',mxDFoF/2)), ['>' sprintf('%6.1f',mxDFoF)]})
ylabel(hBar,'\DeltaF/F')
set(h1,'Position',pos);
ylabel('Neuron number')
h2=subplot(4,1,[3]);
set(h2,'XTickLabel',[],'YTick',[0 20],'LineWidth',2)
set(h2,'TickLength',[0.005 0.025],'TickDir','out')
ylabel('Counts')
mx=max(sum(rasterClean,2));
ylim([-mx/10 mx])
xlim([1 numFrames])
set(gca,'units','points')
p=get(gca,'position');
hold on
count=1;
if strcmp(areThereTags,'Yes')
    for i=1:size(tagColors,1)
        for j=1:length(tagsData.tagStart{i})
            start=tagsData.tagStart{i}(j);
            stop=tagsData.tagEnd{i}(j);
            if j==1
                hleg(count)=line([start stop],[-mx/20 -mx/20],'Color',tagColors(i,:),'LineWidth',p(4)/10);
                leg{count}=['tag ' num2str(count)];
                count=count+1;
            else
                line([start stop],[-mx/20 -mx/20],'Color',tagColors(i,:),'LineWidth',p(4)/10);
                
            end
        end
    end
    h = legend(hleg,leg{:});
end

bar(sum(rasterClean,2),1.1,'k')
hold on; plot([1 numFrames],[confSynchBinary confSynchBinary],'r-','LineWidth',2);
freezeColors
h3=subplot(4,1,4);
set(h3,'LineWidth',2);
set(h3,'TickLength',[0.005 0.025],'TickDir','out')
hold on;

plotActivations(1:length(assembliesCells));
xlabel('Frame number'); ylabel('Assembly activation'); set(gca,'YAxisLocation','right');
linkaxes([h1 h2 h3],'x')
posIm=get(gca,'Position');
width=.1; height=.05;
x=max(.01,posIm(1)-width-.01); y=posIm(2)+.1;
uicontrol('Style','pushbutton','String','Subset','CallBack',{@get_Subset},'Units','normalized','position',[x y width height]);
y=y+height+.01;
uicontrol('Style','pushbutton','String','All','CallBack',{@get_All},'Units','normalized','position',[x y width height]);
y=y+height+.01;
uicontrol('Style','text','Units','normalized','position',[x y width height],'String','Which assemblies?')



set(gcf,'toolbar','figure')

%%

    function get_All(~,~,~)
        set(gcf,'CurrentAxes',h3);
        cla
        plotActivations(1:length(assembliesCells));
    end

    function get_Subset(~,~,~)
        prompt = {['Comma-separated list of assemblies to be displayed (e.g.: 3 6 8)']};
        answer = inputdlg(prompt, ['Select assemblies from a total of ' num2str(length(assembliesCells))], 1);
        subset=str2num(answer{1});
        set(gcf,'CurrentAxes',h3);
        cla
        plotActivations(subset);
    end
    function plotActivations(toPlot,~)
        set(gcf,'CurrentAxes',h3);
        colorsMatch=jet(length(assembliesCells));
        for k=1:length(toPlot)
            plot(h3,1:size(matchIndexTimeSeriesSignificant,2),matchIndexTimeSeriesSignificant(toPlot(k),:),'Color',colorsMatch(toPlot(k),:))
        end
        
    end
end