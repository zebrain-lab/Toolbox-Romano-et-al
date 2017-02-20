%%%%%%  CHOOSE INPUT FILE
clear all
[filename,pathname] = uigetfile({'*_ALL_CELLS.mat';'*_ALL_CELLS.MAT'},'Open file with fluorescence time series', 'MultiSelect', 'off');
filenameALL_CELLS=fullfile(pathname,filename);

cut=strfind(filenameALL_CELLS,'_ALL_CELLS.mat');
filenameARTIFACTS=[filenameALL_CELLS(1:cut-1) '_ARTIFACTS.mat'];
load(filenameARTIFACTS);

%%
%%%%%%  PARAMETERS
plotFlag=0;
prompt = {'Frequency of imaging (frames per second, Hz)','Fluorescence decay time constant of reporter (\tau , seconds)'};
dlg_title = 'Imaging parameters';
num_lines = 1;
def = {'4','0.38'};
opts.Interpreter='tex';% x=1:size(raster,1);
answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
params.fps= str2num(answer{1});
params.tauDecay=str2num(answer{2});

plotFlag=0;
prompt = {'Minimal number of pixels per ROI','Minimal ROI fluorescence relative to baseline level (%)','Maximal sudden decrease in ROI baseline fluorescence (z-score of baseline fluorescence variations)'};
dlg_title = 'Parameters for control of data sanity';
num_lines = 1;
def = {'5','25','-4'};

answer = inputdlg(prompt,dlg_title,num_lines,def);
params.cutOffPixels= str2num(answer{1});
params.cutOffIntensity= str2num(answer{2});
params.cutOffDev= str2num(answer{3});



%%
data=load(filenameALL_CELLS);
imageAvg=data.avg;
cut=strfind(filenameALL_CELLS,'_ALL_CELLS.mat');
outputFile=[filenameALL_CELLS(1:cut-1) '_RASTER.mat'];

ansMethod = questdlg('Subtract the fluorescence from surrounding neuropile?', 'Neuropile fluorescence contamination correction', 'Yes','No','Yes');
params.neuropileSubtraction=ansMethod;
if strcmp(params.neuropileSubtraction,'Yes')
    prompt = {'Set the coefficient \alpha (0 to 1) for neuropile subtraction (Fcorrected = F - \alpha * Fneuropile)'};
    dlg_title = 'Set neuropile correction';
    num_lines = 1;
    def = {'0.4'};
    opts.Interpreter='tex';
    answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
    params.alpha= str2num(answer{1});
    
end

ansMethod = questdlg('Choose method to estimate baseline fluorescence', 'Baseline fluorescence (F0) calculation', 'Average fluorescence on time window','Smooth slow dynamics','Average fluorescence on time window');
params.fluoBaselineCalculation.Method=ansMethod;

if strcmp(params.fluoBaselineCalculation.Method,'Average fluorescence on time window')
    prompt = {'Select t0 (in seconds)','Select t1 (in seconds)'};
    dlg_title = 'Select time window (from t0 to t1) for calculation of F0 ';
    num_lines = 1;
    def = {'0','5'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    params.fluoBaselineCalculation.t0= str2num(answer{1});
    params.fluoBaselineCalculation.t1= str2num(answer{2});
    
end


disp('Step 1: Sanity test for ROIs and calculation relative fluorescence variation.')
[fluoTraces,F0,smoothBaseline,deletedCells]=SanityTest(filenameALL_CELLS,params);

deltaFoF=(fluoTraces-F0)./F0;

disp('Step 2: Calculation of noise level in baseline fluorescence ROIs.')
[deltaFoF, mu, sigma, params]=EstimateBaselineNoise(filenameALL_CELLS, deltaFoF, params);

deltaFoF(logical(movements),:)=NaN;

ansMethod = questdlg('Choose method for detection on significant trasients', 'Select method', 'Static threshold', 'Dynamic threshold', 'Import data','Static threshold');
params.methodSignificatTransients=ansMethod;

disp('Step 3: Detection of significant fluorescence transients.')
if strcmp(ansMethod,'Static threshold')
    prompt = {'Minimal \DeltaF/F0 * 1/\sigma of significant ROI fluorescence transient'};
    dlg_title = 'Parameters for inference signinficant fluorescence transients';
    num_lines = 1;
    def = {'3'};
    opts.Interpreter='tex';
    answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
    params.deltaFoFCutOff= str2num(answer{1});
    
    disp('Step 3.1: Producing raster plot.')
    raster=double(bsxfun(@gt,deltaFoF,params.deltaFoFCutOff*sigma+mu));
    
    disp('Step 3.2: Saving and quitting.')
    save(outputFile,'raster', 'params','deltaFoF', 'deletedCells','movements', 'mu', 'sigma','imageAvg', 'F0');
    
elseif strcmp(ansMethod,'Dynamic threshold')
    
    prompt = {'Minimal cofindence that ROI fluorescence transient is not noise (%)'};
    dlg_title = 'Parameter for inference signinficant fluorescence transients';
    num_lines = 1;
   
    def = {'95'};
    opts.Interpreter='tex';
    answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
   
    params.confCutOff= str2num(answer{1});
    
    disp('Step 3.1: Estimating noise model.')
    [densityData, densityNoise, xev, yev] = NoiseModel(filenameALL_CELLS, deltaFoF, sigma, movements, plotFlag);
    [mapOfOdds] = SignificantOdds(deltaFoF, sigma, movements, densityData, densityNoise, xev, params, plotFlag);
    
    disp('Step 3.2: Producing raster plot.')
    
    [raster, mapOfOddsJoint]=Rasterize(deltaFoF, sigma, movements, mapOfOdds, xev, yev, params);
    
    disp('Step 3.1: Saving and quitting.')
    save(outputFile,'raster', 'params','deltaFoF', 'deletedCells','movements', 'mu', 'sigma', 'mapOfOdds', 'mapOfOddsJoint', 'xev', 'yev', 'densityData', 'densityNoise','imageAvg', 'F0');
    
elseif strcmp(ansMethod,'Import data')
    [filename,pathname] = uigetfile({'*.mat';'*.MAT'},'Open file with data to import', 'MultiSelect', 'off');
    dataImport=load(fullfile(pathname,filename));
    raster=dataImport.significantTransients;
    
    disp('Step 3.1: Saving and quitting.')
    save(outputFile,'raster', 'params','deltaFoF', 'deletedCells','movements', 'mu', 'sigma','imageAvg', 'F0');
    
end

ansPlot = questdlg('Plot examples of significant fluorescent traces?', 'Plot results', 'Yes','No','Yes');

if strcmp(ansPlot,'Yes')
    
    ansContinue='Yes';
    alreadyPlot=[];
    [a,b]=sort(nansum(raster),'descend');
    t=linspace(0,(size(raster,1)-1)/params.fps,size(raster,1));
    z=zeros(size(t));
    while strcmp(ansContinue,'Yes')

        prompt = {'Number of traces per plot','\DeltaF/F0 plot scale'};
        dlg_title = 'Plot parameters';
        num_lines = 1;
    
        def = {'20','1'};
        opts.Interpreter='tex';
        answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
        
        totalPlot= str2num(answer{1});
        scalePlot= str2num(answer{2});
        
        
        numPlotCells=1:totalPlot;
        y=bsxfun(@plus,deltaFoF(:,b(numPlotCells)),scalePlot*[1:size(deltaFoF(:,b(numPlotCells)),2)]);
        col = raster(:,b(numPlotCells));
        figure; hold on
        for i=1:size(col,2)
            surface([t(1:end-1);t(1:end-1)],[i;i],[y(1:end-1,i)';y(1:end-1,i)'],[col(1:end-1,i)';col(1:end-1,i)'],'facecol','no','edgecol','interp','linew',1);
        end
        posBar=t(end)*1.05;
        plot3([posBar;posBar],[1;1],[1;1+scalePlot],'k','linewidth',5);
        colormap([0 0 0; 1 0 0])
        xlim([t(1) t(end)*1.07]); xlabel('Time (s)')
        set(gca,'View', [0 40],'ycolor','w','zcolor','w','TickLength',[.01; .01])
        h=text(posBar+5,1,1,[num2str(scalePlot) ' \DeltaF/F0']);
        set(h, 'rotation', 90)
        for i=1:size(col,2)
            text(t(1)-25,i,scalePlot*i,num2str(b(i)));
        end
        
        ansContinue = questdlg('Continue with more traces?', 'Plot results', 'Yes','No','Yes');
        b(1:totalPlot)=[];

    end
    
end

