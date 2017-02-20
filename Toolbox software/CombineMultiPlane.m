clear all
%%
%%%%%%% INPUT FILES
[filename,pathname] = uigetfile({'*.mat';'*.MAT'},'Open file with raster data', 'MultiSelect', 'on');
for i=1:length(filename)
    filenameRASTER{i}=fullfile(pathname,filename{i});
    cutName=strfind(filenameRASTER{i},'_RASTER.mat');
    filenameALL_CELLS{i}=[filenameRASTER{i}(1:cutName-1) '_ALL_CELLS.mat'];
end


dataRaster=load(filenameRASTER{1});
dataAllCells=load(filenameALL_CELLS{1});
raster=[];
deltaFoF=[];
fromPlane=[];
movements=zeros(size(dataRaster.movements,1),1);
cells=cell(1,0);
cell_per=cell(0,1);
avg=zeros(size(dataAllCells.avg,1),size(dataAllCells.avg,2),length(filename));
for i=1:length(filename)
    
    dataRaster=load(filenameRASTER{i});
    dataAllCells=load(filenameALL_CELLS{i});
    
    raster=[raster, dataRaster.raster];
    deltaFoF=[deltaFoF, dataRaster.deltaFoF];
    movements=movements+ dataRaster.movements;
    cells=cat(2,cells, dataAllCells.cells);
    cell_per=cat(1,cell_per, dataAllCells.cell_per);
    avg(:,:,i)=dataAllCells.avg;
    fromPlane=[fromPlane; i*ones(size(dataRaster.deltaFoF,2),1)];
end

save([filenameRASTER{1}(1:cutName-1) '_MultiPlane_RASTER.mat'],'raster','deltaFoF','movements');
save([filenameRASTER{1}(1:cutName-1) '_MultiPlane_ALL_CELLS.mat'],'cells','cell_per','fromPlane','avg');


