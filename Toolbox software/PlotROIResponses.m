function plotRows=PlotROIResponses(numCell,traces,onsets,offsets,responses,responsesStd,mapData,im,cellper,roiHSVRescaled,prevFrame,eventDuration,params,plotRows)

mapParams=nan(length(mapData),1);
mapFrames=cell(length(mapData),1);
mapLabels=cell(length(mapData),1);
numStim=length(mapData);
for i=1:numStim
    mapParams(i)=mapData{i}.value;
   
    mapFrames{i}=round(mapData{i}.onsetTime*params.fps);
    mapLabels{i}=str2num(mapData{i}.label);
end
screensize = get( 0, 'Screensize' );

cont=1;
if isempty(numCell)
    disp('You clicked outside a ROI, please click inside')
    cont=0;
end
plotCols=ceil((numStim+1)/plotRows);

if cont
    figure('Name',['ROI # ' num2str(numCell)]);
    uicontrol('Style','pushbutton','String','Plot order','CallBack','plotRows= plotOrder(mapData); plotRows=PlotROIResponses(numCell,traces,onsets,offsets,responses,responsesStd,mapData,im,cellper,roiHSVRescaled,prevFrame,eventDuration,params,plotRows);','Position',[screensize(3)*0.05 screensize(4)*8/10 90 30]);
    set(gcf,'Color','w')
    set(gcf,'toolbar','figure');
    
    yl=zeros(numStim,2);
    for i=1:numStim
       
        subplot(plotRows*2,plotCols,i)
        hold on
        
        
        plot(traces{numCell,i}','color',[.5 .5 .5]);
        numTrials=size(traces{numCell,i},1);
        for j=1:numTrials
            if ~isempty(onsets{numCell,i})
                
                for k=1:length(onsets{numCell,i}{j})
                    plot(onsets{numCell,i}{j}(k):offsets{numCell,i}{j}(k),traces{numCell,i}(j,onsets{numCell,i}{j}(k):offsets{numCell,i}{j}(k)),'color',[1 0 0],'LineWidth',.5)
                end
            end
            
        end
        plot(mean(traces{numCell,i})','k','LineWidth',2)
        yl(i,:)=get(gca,'YLim');
        axis square
        set(gca,'Visible','off')
    end
    for i=1:numStim
       
        subplot(plotRows*2,plotCols,i)
        plot([prevFrame prevFrame+eventDuration],[min(yl(:,1)) min(yl(:,1))],'b','LineWidth',5)
        ylim([min(yl(:,1)) max(yl(:,2))])
        text(1,max(yl(:,2))*0.9,mapData{i}.label)
        
    end
    
    subplot(plotRows*2,plotCols,i+1)
    rg=max(yl(:,2))-min(yl(:,1));
    plot([2 2+params.fps],[min(yl(:,1))+.2*rg min(yl(:,1))+.2*rg],'k','LineWidth',5); hold on
    text(2,min(yl(:,1)),'1s');
    plot([2 2],[min(yl(:,1))+.5*rg min(yl(:,1))+rg],'k','LineWidth',5)
    h=text(2+5,min(yl(:,1))+.5*rg,[num2str(.5*rg) ' \DeltaF/F0'],'Interpreter','tex');
    set(h, 'rotation', 90)
    ylim([min(yl(:,1)) max(yl(:,2))])
    xlim([1 length(traces{numCell,1})])
    axis square
    set(gca,'Visible','off')
    
    subplotsqueeze(gcf, 1.1)
   
    remaining=1:plotRows*2*plotCols;
    remaining(1:length(remaining)/2)=[];
   
    remaining=reshape(remaining,plotCols,plotRows)';
  
    subplot(plotRows*2,plotCols,sort(reshape(remaining(:,1:floor(size(remaining,2)/2)),[],1)));
    imagesc(im); hold on;
    shading flat; colormap gray; axis image; set(gcf,'color','w'); set(gca,'YTick',[],'Xtick',[])
    title(['ROI #' num2str(numCell)])
    
    ind=numCell;
    verts=[cellper{ind}(:,1), cellper{ind}(:,2)];
    faces=1:1:length(verts);
    p=patch('Faces',faces,'Vertices',verts,'FaceColor',hsv2rgb(roiHSVRescaled(ind,:)),'EdgeColor',hsv2rgb(roiHSVRescaled(ind,:)),'FaceAlpha',roiHSVRescaled(ind,3),'EdgeAlpha',roiHSVRescaled(ind,3));
    
    remaining=setdiff(reshape(remaining,[],1),sort(reshape(remaining(:,1:floor(size(remaining,2)/2)),[],1)));
    
    subplot(plotRows*2,plotCols,remaining);
    if length(mapParams)==1
        errorbar(mapParams,responses(:,numCell),responsesStd(:,numCell)/sqrt(numTrials),'k')
        ylim([0 (responses(:,numCell)+responsesStd(:,numCell)/sqrt(numTrials))*1.2])
        xlim(sort([mapParams*.9 mapParams*1.1]))
    else
        boundedline(mapParams,responses(:,numCell),responsesStd(:,numCell)/sqrt(numTrials),'k');
    end
    set(gca,'LineWidth',2)
    xlabel('Mapping parameter'); ylabel('Average \DeltaF/F0','Interp','tex')
    axis square
 
end
