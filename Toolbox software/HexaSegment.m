function [ROI,MaskROI]=HexaSegment(img,Mask,L)

%% Hexagonal mapping of space using voronoi:

% Generate hexagonal grid
Nx=size(img,2);Ny=size(img,1);
% Generate hexagonal grid
Rad3Over2 = sqrt(3) / 2;
xP=0:L:Nx;
xI=(0:L:Nx)+0.5*L;
y=0:L*sqrt(3) / 2:Ny;

X=nan(numel(y),numel(xP));
Y=nan(numel(y),numel(xP));
for j=1:numel(y)
    if mod((j-1),2)==0
        X(j,:)=xP;
    else
        X(j,:)=xI;
    end
    Y(j,:)=(j-1)*L*sqrt(3) / 2;
end

%% Create mask using nearest neighboor

% Iterate on pixel:
Numpix=numel(img);
clear Cell
Cell=nan(size(img));
for i=1:size(img,1)
    for j=1:size(img,2)
        
        % Assign Pixel To Celler
        Dist=(j-X).^2+(i-Y).^2;
        [~,id]=min(Dist(:));
        
        Cell(i,j)=id;
    end
end

%% Color Celler:
NumCell=numel(unique(Cell));
Perm=randperm(numel(unique(Cell)));

CellCol=nan(size(Cell));
for c=1:NumCell
    id=find(Cell==c);
    CellCol(id)=Perm(c);
end

%% Intersect Mask & ROI

ROI=[];
ROIMask=zeros(size(img));
k=1;
for c=1:NumCell
    id=find(Cell==c);
    if sum(Mask(id))==numel(id)
        ROI{k}=id;
        ROIMask(id)=k;
        k=k+1;
    end
end

%% Overlay voxel and brain:
color=jet(numel(ROIMask));
color=color(randperm(numel(ROIMask)),:);

MaskROI=zeros(size(img,1),size(img,2),3);
for i=1:numel(ROI)
    [I,J] = ind2sub(size(img),ROI{i});
    for k=1:numel(I)
        MaskROI(I(k),J(k),:)=color(i,:);
    end
end

