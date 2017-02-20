function [myColors, valGam]=changeGamma(hObject, valContr,myColorsOrig)
%       
        valGam = get(hObject,'Value');
     
       
        myColors=[linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr));linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr));linspace(0,1,256-floor(valContr)).^ valGam,ones(1,floor(valContr))]';
        
        
        
end


