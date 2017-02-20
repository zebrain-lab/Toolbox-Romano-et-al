function ROItoPlot=lookForROI(im,totalCells,pixROIs,hIm)
[colROI, rowROI] = getpts(hIm);
pixROI=sub2ind(size(im),round(rowROI),round(colROI));
ROItoPlot=[];
for i=1:totalCells
    if ismember(pixROI,pixROIs{i})
        ROItoPlot=i;
    end
end
