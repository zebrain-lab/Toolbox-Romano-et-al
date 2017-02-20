
[filename,pathname] = uigetfile({'*_CLUSTERS.mat';'*_CLUSTERS.MAT'},'Open file with assemblies data', 'MultiSelect', 'off');
filenameCLUSTER=fullfile(pathname,filename);

load(filenameCLUSTER);

cutName=strfind(filenameCLUSTER,'_CLUSTERS.mat');
dataAllCells=load([filenameCLUSTER(1:cutName-1) '_ALL_CELLS.mat']);
load([filenameCLUSTER(1:cutName-1) '_RASTER.mat']);
rasterAnalog=deltaFoF;
rasterAnalog(~raster)=0;
rasterAnalog(isnan(sum(rasterAnalog,2)),:)=0;
numColors=1000;
colorstopaint=flipud(hsv(numColors));
colorstopaint=colorstopaint(numColors*.161:numColors,:);
numColors=size(colorstopaint,1);
nInterp=numColors;
%% Plotting assemblies together according to a particular topography
orderSource = questdlg('How do we order assemblies?', 'User input','Similarity of dynamics','Along topographical axis', 'Similarity of dynamics');

if strcmp(orderSource,'Similarity of dynamics')
    D = pdist(matchIndexTimeSeriesSignificant,'correlation');
    Z = linkage(matchIndexTimeSeriesSignificant,'average','correlation');
    assembliesOrdered{1} = fliplr(optimalleaforder(Z,D));
    projectionOrderedTotal=linspace(0,1,size(matchIndexTimeSeriesSignificant,1));
    
elseif strcmp(orderSource,'Along topographical axis')
    
    orderSourceTopo = questdlg('Select topographical method', 'User input', 'Along curve', 'Left-right', 'Bottom-top','Along curve');
     
    clear X Y
    
    if strcmp(orderSourceTopo,'Along curve')
        prompt = {'How many regions?'};
        dlg_title = 'Divide image in regions';
        num_lines = 1;
        def = {'1'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        numStruct=str2num(answer{1});
        pointsCurve=cell(numStruct,1);
        mask=zeros(numStruct,size(dataAllCells.avg,1),size(dataAllCells.avg,2));
        for i=1:numStruct
             h=figure('Name','Draw a mask for each region');
            imagesc(dataAllCells.avg)
            shading flat; colormap gray; axis image; set(gcf,'color','w');
            set(gcf, 'Position', get(0,'Screensize'));
            drawnow
            BW = roipoly;
            mask(i,:,:)=BW;
        end
        
         c=linspace(0,1,numStruct*nInterp);
         yy=linspace(min(c),max(c),size(colorstopaint,1));
         cm = spline(yy,colorstopaint',c);
         cm(cm>1)=1;
         cm(cm<0)=0;
        
        for i=1:numStruct
            satisfied='No';
            while strcmp(satisfied,'No')
                h=figure('Name','Draw a curve to establish an axis for the topography of that structure');
                imagesc(dataAllCells.avg);
                shading flat; colormap gray; axis image; set(gcf,'color','w');
                set(gcf, 'Position', get(0,'Screensize'));
                drawnow
                [Y{i},X{i}]=getpts(h);
                hold on;
                
                monoX=(all(diff(X{i})<=0) | all(diff(X{i})>=0));
                
                monoY=(all(diff(Y{i})<=0) | all(diff(Y{i})>=0));
                
                temp=find([monoX monoY]);
                if ~isempty(temp)
                    useThis=temp(1);
                else
                    disp('Retry drawing curve')
                    close(h)
                    continue
                end
                
                if useThis==2
                    Xb=X{i}; Yb=Y{i};
                    X{i}=Yb; Y{i}=Xb;
                end
                
                if unique(X{i})~=length(X{i})
                    [junk,indNonRepeat]=unique(X{i});
                    indToCorrect=setdiff(1:length(X{i}),indNonRepeat);
                    X{i}(indToCorrect)=X{i}(indToCorrect)+0.1;
                end
                
                x=linspace(min(X{i}),max(X{i}),nInterp);
                y=interp1(X{i},Y{i},x,'pchip');
                pointsCurve{i}=[x;y]';
                if useThis==2
                    pointsCurve{i}=fliplr(pointsCurve{i});
                end
                
                if useThis==2
                    [junk which]=min([abs(pointsCurve{i}(1,2)-X{i}(1)) abs(pointsCurve{i}(1,2)-X{i}(end))]);
                    if which==2
                        pointsCurve{i}=flipud(pointsCurve{i});
                    end
                else
                    [junk which]=min([abs(pointsCurve{i}(1,2)-Y{i}(1)) abs(pointsCurve{i}(1,2)-Y{i}(end))]);
                    if which==2
                        pointsCurve{i}=flipud(pointsCurve{i});
                    end
                end
              
                
                                
                for j=1:nInterp-1
                    indColor=j+(i-1)*(nInterp);
                    line([pointsCurve{i}(j,2) pointsCurve{i}(j+1,2)],[pointsCurve{i}(j,1) pointsCurve{i}(j+1,1)],'color',[cm(:,indColor)],'LineWidth',5);
                    
                end
                pause
                satisfied = questdlg('Satisfied with the curve?', 'Question', 'Yes', 'No', 'No');
                
            end
        end
    elseif strcmp(orderSourceTopo,'Left-right')
        numStruct=1;
        pointsCurve=cell(numStruct,1);
        mask=ones(numStruct,size(dataAllCells.avg,1),size(dataAllCells.avg,2));
        
        pointsCurve{1}(:,2)=linspace(1,size(dataAllCells.avg,2),nInterp);
        pointsCurve{1}(:,1)=round(size(dataAllCells.avg,1)/2)*ones(nInterp,1);
    else
        numStruct=1;
        pointsCurve=cell(numStruct,1);
        mask=ones(numStruct,size(dataAllCells.avg,1),size(dataAllCells.avg,2));
        pointsCurve{1}(:,2)=round(size(dataAllCells.avg,2)/2)*ones(nInterp,1);
        pointsCurve{1}(:,1)=linspace(size(dataAllCells.avg,1),1,nInterp);
    end
    
    % we check the correspondance of assemblies to the structures. the structure
    % with the biggest number of cells of an assembly "owns" the assembly
    
    inStruct=cell(numStruct,1);
    for k=1:length(assembliesCells)
        mode=k;
        cellsInStruct=zeros(numStruct,1);
        for j=1:length(assembliesCells{mode})
            numCell=assembliesCells{mode}(j);
            testCell=zeros(size(dataAllCells.avg));
            testCell(dataAllCells.cells{numCell})=1;
            pixStruct=zeros(numStruct,1);
            for i=1:numStruct
                pixStruct(i)= sum(sum(squeeze(mask(i,:,:)) & testCell));
            end
            [junk,ind]=max(pixStruct);
            cellsInStruct(ind)= cellsInStruct(ind)+1;
        end
        [junk,ind]=max(cellsInStruct);
        inStruct{ind}= [inStruct{ind} mode];
    end
    
    % We calculate the centroids of the assemblies and project them on the axis
    % of topography
    
    centerAssembly=cell(numStruct,1);
    for k=1:numStruct
        for i=1:length(inStruct{k})
            centerAssembly{k}(i,1:2)=[0 0];
            mode=inStruct{k}(i);
            for j=1:length(assembliesCells{mode})
                test=zeros(size(dataAllCells.avg));
                test(dataAllCells.cells{assembliesCells{mode}(j)})=1;
                statsCell=regionprops(test,'Centroid');
                centerAssembly{k}(i,1:2)=(statsCell.Centroid/length(assembliesCells{mode}) + centerAssembly{k}(i,1:2));
            end
        end
    end
    
    projections=cell(numStruct,1);
    order=cell(numStruct,1);
    projectionOrdered=cell(numStruct,1);
    for j=1:numStruct
        for i=1:length(inStruct{j})
            idx=knnsearch(pointsCurve{j}, fliplr(centerAssembly{j}(i,1:2)));
            projections{j}(i)=sum(sqrt(diff(pointsCurve{j}(1:idx,1)).^2+diff(pointsCurve{j}(1:idx,2)).^2))/sum(sqrt(diff(pointsCurve{j}(:,1)).^2+diff(pointsCurve{j}(:,2)).^2));
            
        end
        [projectionOrdered{j},order{j}]=sort(projections{j},'ascend');
    end
    
    
    assembliesOrdered{1}=[];
    projectionOrderedTotal=[];
    separation=.05;
    for i=1:numStruct
        assembliesOrdered{1}=[assembliesOrdered{1} inStruct{i}(order{i})];
        temp=projectionOrdered{i};
       
        if i==1
            temp=temp/numStruct;
        else
           
            temp=(1/numStruct)+temp/numStruct;
        end
        
        projectionOrderedTotal=[projectionOrderedTotal temp];
    end
end

corrs=corr(rasterAnalog);
transp=cell(size(rasterAnalog,2),1);
for i=1:length(assembliesCells)
    corrsAss=corrs(assembliesCells{i},assembliesCells{i});
    
    normColor=mean(corrsAss);
    normColor=normColor-min(normColor);
    normColor=normColor+.1;
    normColor=normColor/max(normColor);
    for j=1:length(assembliesCells{i})
        transp{assembliesCells{i}(j)}=[transp{assembliesCells{i}(j)} normColor(j)];
    end
end

%    now we plot the assemblies together

figure('Name','Topographical organization of assemblies');
totalModes=length(assembliesOrdered{1});
cellsToPlot=unique(cat(1,[assembliesCells{assembliesOrdered{1}}]));


orderOfCell=zeros(dataAllCells.cell_number,1);
if isempty(cellsToPlot)
    text(.5,.5,'No assemblies','Fontsize',14,'FontWeight','Bold');
else
    imagesc(dataAllCells.avg); hold on;
    shading flat; colormap gray; axis image; set(gcf,'color','w'); set(gca,'YTick',[],'Xtick',[])
    set(gcf, 'Position', get(0,'Screensize'));
    for j=1:length(cellsToPlot)
        ind=cellsToPlot(j);
        clear inMode
        for i=1:totalModes
            numMode=assembliesOrdered{1}(i);
            inMode(i)=ismember(cellsToPlot(j),assembliesCells{numMode});
        end
        [junk, belongsMode]=find(inMode);
        
        
        verts=[dataAllCells.cell_per{ind}(:,1), dataAllCells.cell_per{ind}(:,2)];
        faces=1:1:length(verts);
        indColor=round(sum(projectionOrderedTotal(belongsMode))/length(belongsMode)*(numColors-1)+1);
        
        p=patch('Faces',faces,'Vertices',verts,'FaceColor',colorstopaint(indColor,:),'EdgeColor',colorstopaint(indColor,:),'FaceAlpha',mean(transp{cellsToPlot(j)}),'EdgeAlpha',mean(transp{cellsToPlot(j)}));
        orderOfCell(cellsToPlot(j)) = indColor;
        
    end
end
pause(2)
if strcmp(orderSource,'Similarity of dynamics') 
    export_fig([filenameCLUSTER(1:cutName-1) '_assemblies_topo_temporal.png']);
else
    export_fig([filenameCLUSTER(1:cutName-1) '_assemblies_topo_spatial.png']);
end

[junk,orderOfCells]=sort(orderOfCell,'descend');
cutCell=find(junk<=0,1,'first');
orderOfCellsInEnsembles=orderOfCells(1:cutCell);
%%
save([filenameCLUSTER(1:cutName-1) '_ORDER_TOPO.mat'],'orderOfCells', 'orderOfCellsInEnsembles','assembliesOrdered')

