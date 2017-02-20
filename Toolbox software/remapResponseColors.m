function [roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,mapData,roiHSVRescaled,peakResponseLabels,peakParameterStdLabels,cellper,mx,cutUp,cutDown,offset,transpa,hIm]=remapResponseColors(roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,remap,cellper,cutUp,cutDown,offset,mapData,pixROIs,transpa)

roiHSVRescaled=roiHSV;
peakResponseLabels=peakResponse;
peakParameterStdLabels=peakParameterStd;
numStim=length(mapData);
if remap
    figure
    for i=2:3
        if i==2 & numStim==1
            continue
        end
        subplot(1,2,i-1);
        roiHSVRescaled(:,i)=roiHSVRescaled(:,i)-min(roiHSVRescaled(:,i));
        roiHSVRescaled(:,i)=roiHSVRescaled(:,i)/max(roiHSVRescaled(:,i));
        plot(sort( roiHSVRescaled(:,i)),'-ok')
        axis square
        hold on
        legend('Original')
        if i==2
            ylabel('Normalized saturation (inversely prop. to tuning width)'); xlabel('Sorted ROIs')
        else
            ylabel('Normalized value (proportional to peak response)'); xlabel('Sorted ROIs')
        end
    end
    drawnow
    if numStim>1
        prompt = {'Lower bound for saturation (0 to 1)','Upper bound for saturation (0 to 1)','Offset for saturation (0 to 1)','Lower bound for saturation (0 to 1)','Upper bound for saturation (0 to 1)','Offset for saturation(0 to 1)'};
        dlg_title = 'Parameters for color remapping';
        num_lines = 1;
        def = {'0','1','0','0','1','0'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        cutDown{2}= str2num(answer{1});
        cutUp{2}=str2num(answer{2});
        offset{2}= str2num(answer{3});
        cutDown{3}= str2num(answer{4});
        cutUp{3}= str2num(answer{5});
        offset{3}= str2num(answer{6});
    else
        prompt = {'Lower bound for saturation (0 to 1)','Upper bound for saturation (0 to 1)','Offset for saturation(0 to 1)'};
        dlg_title = 'Parameters for color remapping';
        num_lines = 1;
        def = {'0','1','0'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        cutDown{3}= str2num(answer{1});
        cutUp{3}= str2num(answer{2});
        offset{3}= str2num(answer{3});
    end
end
for i=2:3
    if i==2
       
        if numStim>1
            indsCutUp=roiHSVRescaled(:,i)>=cutUp{i};
            indsCutDown=roiHSVRescaled(:,i)<=cutDown{i};
            
            roiHSVRescaled(indsCutUp,i)=cutUp{i};
            roiHSVRescaled(indsCutDown,i)=cutDown{i};
            roiHSVRescaled(:,i)=roiHSVRescaled(:,i)-min(roiHSVRescaled(:,i));
            roiHSVRescaled(:,i)=roiHSVRescaled(:,i)/max(roiHSVRescaled(:,i));
            roiHSVRescaled(:,i)=roiHSVRescaled(:,i)*(1-offset{i});
            roiHSVRescaled(:,i)=roiHSVRescaled(:,i)+offset{i};
            peakParameterStdLabels(indsCutDown)=min(peakParameterStdLabels(indsCutDown));
            peakParameterStdLabels(indsCutUp)=max(peakParameterStdLabels(indsCutUp));
        end
        
       
        if remap

            if i==2 & numStim==1
                continue
            end
            subplot(1,2,i-1);
            plot(sort( roiHSVRescaled(:,i)),'-or')
            legend('Original','Rescaled')
        end
        
    end
    
    if i==3
       
        indsCutUp=roiHSVRescaled(:,i)>=cutUp{i};
        indsCutDown=roiHSVRescaled(:,i)<=cutDown{i};
        roiHSVRescaled(indsCutUp,i)=cutUp{i};
        roiHSVRescaled(indsCutDown,i)=cutDown{i};
        roiHSVRescaled(:,i)=roiHSVRescaled(:,i)-min(roiHSVRescaled(:,i));
        roiHSVRescaled(:,i)=roiHSVRescaled(:,i)/max(roiHSVRescaled(:,i));
        roiHSVRescaled(:,i)=roiHSVRescaled(:,i)*(1-offset{i});
        roiHSVRescaled(:,i)=roiHSVRescaled(:,i)+offset{i};
        peakResponseLabels(indsCutUp)=min(peakResponseLabels(indsCutUp));
        peakResponseLabels(indsCutDown)=max(peakResponseLabels(indsCutDown));
        if remap
            
            subplot(1,2,i-1);
            plot(sort( roiHSVRescaled(:,i)),'-or')
            legend('Original','Rescaled')
           
            drawnow
            pause(1)
        end
    end
end
mx=.8;
roiHSVRescaled(:,1)=mx*roiHSVRescaled(:,1); 


[transpa, hIm]=plotResponses(roiHSV,peakParameter,peakResponse,peakParameterStd,im,totalCells,mapData,roiHSVRescaled,peakResponseLabels,peakParameterStdLabels,cellper,pixROIs,mx,remap,transpa);
        