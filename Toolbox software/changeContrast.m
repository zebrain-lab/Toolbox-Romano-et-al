function [myColors, valContr]=changeContrast(hObject, valGam,myColorsOrig)
%         
        valContr = get(hObject,'Value');
       
       
        myColors=[linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr));linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr));linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr))]';
        
        
        
        
end


