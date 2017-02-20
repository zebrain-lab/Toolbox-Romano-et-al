function [scaleLocCont scaleChanged]=changeScale(hObject)
%       
        scaleLocCont = get(hObject,'Value');
        scaleLocCont=round(scaleLocCont);
        if scaleLocCont<1
          
            scaleLocCont=1;
        end
        scaleChanged=1;
       
        
end


