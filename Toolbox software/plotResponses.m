
function [transpa,hIm]=plotResponses(roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,mapData,roiHSVRescaled,peakResponseLabels,peakParameterStdLabels,cellper,pixROIs,mx,remap,transpa,hIm)
                 
figure
screensize = get( 0, 'Screensize' );
set(gcf,'Position',[screensize(1:4)])
hIm=subplot(1,6,1:4);
posIm=get(gca,'Position');
width=.1; height=.1;
x=max(.01,posIm(1)-width-.01);x=.01; y=posIm(2)+.1; 
uicontrol('Style','pushbutton','String','Save','CallBack','saveMap(fnameOut,roiHSV,peakParameter,peakResponse,peakParameterStd,roiHSVRescaled,peakResponseLabels,peakParameterStdLabels,cutUp,cutDown,offset,mapParams,traces,onsets,responsesStd); disp(''Mapping results saved'')','Units','normalized','position',[x y width height]);
y=y+height+.01;
uicontrol('Style','pushbutton','String','Select ROI','CallBack','numCell=lookForROI(im,totalCells,pixROIs,hIm); plotRows=PlotROIResponses(numCell,traces,onsets,offsets,responses,responsesStd,mapData,im,cellper,roiHSVRescaled,prevFrame,eventDuration,params,plotRows);','Units','normalized','position',[x y width height]);
y=y+height+.01;
uicontrol('Style','pushbutton','String','Remap colors','CallBack','remap=1; [roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,mapData,roiHSVRescaled,peakResponseLabels,peakParameterStdLabels,cellper,mx,cutUp,cutDown,offset,transpa,hIm]=remapResponseColors(roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,remap,cellper,cutUp,cutDown,offset,mapData,pixROIs,transpa);','Units','normalized','position',[x y width height]);
y=y+height+.01;
uicontrol('Style','pushbutton','String','Transp. on/off','CallBack','transpa=setdiff([1 0],transpa); remap=1; [transpa,hIm]=plotResponses(roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,mapData,roiHSVRescaled,peakResponseLabels,peakParameterStdLabels,cellper,pixROIs,mx,remap,transpa);','Units','normalized','position',[x y width height]);

set(gcf,'toolbar','figure');

imagesc(im); hold on;
shading flat; colormap gray; axis image; set(gcf,'color','w'); set(gca,'YTick',[],'Xtick',[])

for ind=1:totalCells
    
    verts=[cellper{ind}(:,1), cellper{ind}(:,2)];
    faces=1:1:length(verts);
    if transpa
        patch('Faces',faces,'Vertices',verts,'FaceColor',hsv2rgb(roiHSVRescaled(ind,:)),'EdgeColor',hsv2rgb(roiHSVRescaled(ind,:)),'FaceAlpha',roiHSVRescaled(ind,3),'EdgeAlpha',roiHSVRescaled(ind,3));
    else
        patch('Faces',faces,'Vertices',verts,'FaceColor',hsv2rgb(roiHSVRescaled(ind,:)),'EdgeColor',hsv2rgb(roiHSVRescaled(ind,:)));
    end
    
    
end
freezeColors;
hues=unique(roiHSVRescaled(:,1));

labs=cell(mapData);
for i=1:length(labs)
    labs{i}=mapData{i}.label;
end
cmap=zeros(length(hues),3);
for i=1:length(hues)
        cmap(i,:)=hsv2rgb([hues(i) 1 1]);
end
colormap(cmap);
hc=colorbar;
yCbar=get(hc,'Ylim');
posLabel=linspace(.5/length(labs),1-.5/length(labs),length(labs));
posLabel=yCbar(1)+(yCbar(2)-yCbar(1))*posLabel;
set(hc,'YTick',posLabel,'YTickLabel',labs,'TickDir','out')
ylabel(hc,'Peak mapping parameter');

sats=linspace(max(roiHSVRescaled(:,2)),min(roiHSVRescaled(:,2)),length(unique(roiHSVRescaled(:,2))));
satsLabels=linspace(min(peakParameterStdLabels),max(peakParameterStdLabels),length(unique(roiHSVRescaled(:,2))));
satsTicks=linspace(1,length(sats),5);
satsTicksLabelsTemp=linspace(satsLabels(1),satsLabels(end),5);
satsTicksLabels=cell(length(satsTicksLabelsTemp),1);
for i=1:length(satsTicksLabels)
    satsTicksLabels{i}=num2str(satsTicksLabelsTemp(i),'%3.1f');
    if remap
        if i==1
            satsTicksLabels{i}=['<=' satsTicksLabels{i}];
        end
        if i==length(satsTicksLabels)
            satsTicksLabels{i}=['>=' satsTicksLabels{i}];
        end
    end
end
matSat=zeros(length(sats),length(labs),3);
for i=1:length(hues)
    for j=1:length(sats)
        matSat(j,i,:)=[hues(i) sats(j) 1];
    end
end

vals=linspace(max(roiHSVRescaled(:,3)),min(roiHSVRescaled(:,3)),length(unique(roiHSVRescaled(:,3))));
valsLabels=linspace(max(peakResponseLabels),min(peakResponseLabels),length(unique(roiHSVRescaled(:,3))));
valsTicks=linspace(1,length(vals),5);
valsTicksLabelsTemp=linspace(valsLabels(1),valsLabels(end),5);
valsTicksLabels=cell(length(valsTicksLabelsTemp),1);
for i=1:length(valsTicksLabels)
    valsTicksLabels{i}=num2str(valsTicksLabelsTemp(i),'%3.2f');
    if remap
        if i==1
            valsTicksLabels{i}=['>=' valsTicksLabels{i}];
        end
        if i==length(valsTicksLabels)
            valsTicksLabels{i}=['<=' valsTicksLabels{i}];
        end
    end
end
matVals=zeros(length(vals),length(labs),3);
for i=1:length(hues)
    for j=1:length(vals)
        matVals(j,i,:)=[hues(i) 1 vals(j)];
    end
end


if length(mapData)>1
    hSat=subplot(1,6,5);
   
    imagesc(hsv2rgb(matSat));
    pause(1)
    posCbarHue=get(hc,'Position');
    set(hSat,'XTick',[],'YTick',satsTicks,'YTickLabel',satsTicksLabels,'YAxisLocation','right','TickDir','out','Position',[(posCbarHue(1)+posCbarHue(3))+0.075 posCbarHue(2:4)])
    ylabel(hSat,'Tuning width');
    
    hVals=subplot(1,6,6);
   
    imagesc(hsv2rgb(matVals));
    pause(1)
    posCbarSat=get(hSat,'Position');
    set(hVals,'XTick',[],'YTick',valsTicks,'YTickLabel',valsTicksLabels,'YAxisLocation','right','TickDir','out','Position',[(posCbarSat(1)+posCbarSat(3))+0.075 posCbarSat(2:4)])
    ylabel(hVals,'Response strength (average \DeltaF/F0)','Interpreter','tex');
    drawnow
else
    hVals=subplot(1,6,6);
 
    imagesc(hsv2rgb(matVals));
    pause(1)
    posCbarHue=get(hc,'Position');
    set(hVals,'XTick',[],'YTick',valsTicks,'YTickLabel',valsTicksLabels,'YAxisLocation','right','TickDir','out','Position',[(posCbarHue(1)+posCbarHue(3))+0.075 posCbarHue(2:4)])
    ylabel(hVals,'Response strength (average \DeltaF/F0)','Interpreter','tex');
    drawnow
    
end

   