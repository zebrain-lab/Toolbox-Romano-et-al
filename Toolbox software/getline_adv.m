function varargout = getline_adv(varargin)

global frm_num ;

%MODIFIED VERSION OF GETLINE, TO BE USED WITH ROI AND ROIPOLY_ADV

%GETLINE Select polyline with mouse.
%   [X,Y] = GETLINE(FIG) lets you select a polyline in the
%   current axes of figure FIG using the mouse.  Coordinates of
%   the polyline are returned in X and Y.  Use normal button
%   clicks to add points to the polyline.  A shift-, right-, or
%   double-click adds a final point and ends the polyline
%   selection.  Pressing RETURN or ENTER ends the polyline
%   selection without adding a final point.  Pressing BACKSPACE
%   or DELETE removes the previously selected point from the
%   polyline.
%
%   [X,Y] = GETLINE(AX) lets you select a polyline in the axes
%   specified by the handle AX.
%
%   [X,Y] = GETLINE is the same as [X,Y] = GETLINE(GCF).
%
%   [X,Y] = GETLINE(...,'closed') animates and returns a closed
%   polygon.
%
%   Example
%   --------
%       imshow('moon.tif')
%       [x,y] = getline 
%
%   See also GETRECT, GETPTS.

%   Callback syntaxes:
%        getline('KeyPress')
%        getline('FirstButtonDown')
%        getline('NextButtonDown')
%        getline('ButtonMotion')

%   Grandfathered syntaxes:
%   XY = GETLINE(...) returns output as M-by-2 array; first
%   column is X; second column is Y.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 5.27.4.3 $  $Date: 2004/08/10 01:39:34 $

global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_X GETLINE_Y
global GETLINE_ISCLOSED 


xlimorigmode = xlim('mode');
ylimorigmode = ylim('mode');
xlim('manual');
ylim('manual');

if ((nargin >= 1) && (ischar(varargin{end})))
    str = varargin{end};
    if (str(1) == 'c')
        % getline(..., 'closed')
        GETLINE_ISCLOSED = 1;
        varargin = varargin(1:end-1);
    end
else
    GETLINE_ISCLOSED = 0;
end

if ((length(varargin) >= 1) && ischar(varargin{1}))
    % Callback invocation: 'KeyPress', 'FirstButtonDown',
    % 'NextButtonDown', or 'ButtonMotion'.
    feval(varargin{:});
    return;
end

GETLINE_X = [];
GETLINE_Y = [];

if (length(varargin) < 1)
    GETLINE_AX = gca;
    GETLINE_FIG = ancestor(GETLINE_AX, 'figure');
else
    if (~ishandle(varargin{1}))
        CleanUp(xlimorigmode,ylimorigmode);
        eid = 'Images:getline:expectedHandle';
        error(eid, '%s', 'First argument is not a valid handle');
    end
    
    switch get(varargin{1}, 'Type')
    case 'figure'
        GETLINE_FIG = varargin{1};
        GETLINE_AX = get(GETLINE_FIG, 'CurrentAxes');
        if (isempty(GETLINE_AX))
            GETLINE_AX = axes('Parent', GETLINE_FIG);
        end

    case 'axes'
        GETLINE_AX = varargin{1};
        GETLINE_FIG = ancestor(GETLINE_AX, 'figure');

    otherwise
        CleanUp(xlimorigmode,ylimorigmode);
        eid = 'Images:getline:expectedFigureOrAxesHandle';
        error(eid, '%s', 'First argument should be a figure or axes handle');
    end
end

% Remember initial figure state
old_db = get(GETLINE_FIG, 'DoubleBuffer');
state= uisuspend(GETLINE_FIG);


% CREATING A CUSTOM POINTER
% P = ones(16)+1;
% P(1,:) = 1; P(16,:) = 1;
% P(:,1) = 1; P(:,16) = 1;
% P(1:4,8:9) = 1; P(13:16,8:9) = 1;
% P(8:9,1:4) = 1; P(8:9,13:16) = 1;
% P(5:12,5:12) = NaN; % Create a transparent region in the center

P=ones(16,16)*NaN;
P(8:10,8:10)=1;
P(9,9)=2;



% Set up initial callbacks for initial stage  CHANGE THE TYPE OF CURSOR
set(GETLINE_FIG, ...
    'Pointer','custom','PointerShapeCData',P,'PointerShapeHotSpot',[9 9],...
    'WindowButtonDownFcn', 'getline_adv(''FirstButtonDown'');',...
    'KeyPressFcn', 'getline_adv(''KeyPress'');', ...
    'DoubleBuffer', 'on');

% Bring target figure forward
%fg=figure(GETLINE_FIG);
%maximize(fg);





% Initialize the lines to be used for the drag HERE CHOOSE LINE
% CHARACTERISTICS (E.G. COLOR DASHED, ETC.)
GETLINE_H1 = line('Parent', GETLINE_AX, ...
                  'XData', GETLINE_X, ...
                  'YData', GETLINE_Y, ...
                  'Visible', 'off', ...
                  'Clipping', 'off', ...
                  'Color', 'r', ...
                  'LineStyle', '-');

GETLINE_H2 = line('Parent', GETLINE_AX, ...
                  'XData', GETLINE_X, ...
                  'YData', GETLINE_Y, ...
                  'Visible', 'off', ...
                  'Clipping', 'off', ...
                  'Color', 'y', ...
                  'LineStyle', ':');

% We're ready; wait for the user to do the drag
% Wrap the call to waitfor in try-catch so we'll
% have a chance to clean up after ourselves.
errCatch = 0;
try
    waitfor(GETLINE_H1, 'UserData', 'Completed');
catch
    errCatch = 1;
end

% After the waitfor, if GETLINE_H1 is still valid
% and its UserData is 'Completed', then the user
% completed the drag.  If not, the user interrupted
% the action somehow, perhaps by a Ctrl-C in the
% command window or by closing the figure.

if (errCatch == 1)
    errStatus = 'trap';
    
elseif (~ishandle(GETLINE_H1) || ...
            ~strcmp(get(GETLINE_H1, 'UserData'), 'Completed'))
    errStatus = 'unknown';
    
else
    errStatus = 'ok';
    x = GETLINE_X(:);
    y = GETLINE_Y(:);
    % If no points were selected, return rectangular empties.
    % This makes it easier to handle degenerate cases in
    % functions that call getline.
    if (isempty(x))
        x = zeros(0,1);
    end
    if (isempty(y))
        y = zeros(0,1);
    end
end

% Delete the animation objects
if (ishandle(GETLINE_H1))
    delete(GETLINE_H1);
end
if (ishandle(GETLINE_H2))
    delete(GETLINE_H2);
end

% Restore the figure's initial state
if (ishandle(GETLINE_FIG))
   uirestore(state);
   set(GETLINE_FIG, 'DoubleBuffer', old_db);
end

CleanUp(xlimorigmode,ylimorigmode);

% Depending on the error status, return the answer or generate
% an error message.
switch errStatus
case 'ok'
    % Return the answer
    if (nargout >= 2)
        varargout{1} = x;
        varargout{2} = y;
    else
        % Grandfathered output syntax
        varargout{1} = [x(:) y(:)];
    end
    
case 'trap'
    % An error was trapped during the waitfor
    eid = 'Images:getline:interruptedMouseSelection';
    error(eid, '%s', 'Interruption during mouse selection.');
    
case 'unknown'
    % User did something to cause the polyline selection to
    % terminate abnormally.  For example, we would get here
    % if the user closed the figure in the middle of the selection.
    eid = 'Images:getline:interruptedMouseSelection';
    error(eid, '%s', 'Interruption during mouse selection.');
end

%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress %#ok
global BORRAR FINISH_ROI
global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_PT1 
global GETLINE_ISCLOSED
global GETLINE_X GETLINE_Y
global frm_num hlv line_lim fast frames img;

key = double(get(GETLINE_FIG, 'CurrentCharacter'));%for ascii characters
if isempty(key)
    key=get(GETLINE_FIG, 'Currentkey'); %for special characters non ascii
end
    
    
switch key
case {char(27)}  % esc  to
    % remove the previously selected point
    switch length(GETLINE_X)
    case 0
        % nothing to do
    case 1
        GETLINE_X = [];
        GETLINE_Y = [];
        % remove point and start over
        
        set([GETLINE_H1 GETLINE_H2], ...
                'XData', GETLINE_X, ...
                'YData', GETLINE_Y);
        set(GETLINE_FIG, 'WindowButtonDownFcn', ...
                'getline_adv(''FirstButtonDown'');', ...
                'WindowButtonMotionFcn', '');
    otherwise
        % remove last point
       
        if (GETLINE_ISCLOSED)
            GETLINE_X(end-1) = [];
            GETLINE_Y(end-1) = [];
        else
            GETLINE_X(end) = [];
            GETLINE_Y(end) = [];
        end
        set([GETLINE_H1 GETLINE_H2], ...
                'XData', GETLINE_X, ...
                'YData', GETLINE_Y);
    end
    

        
    
 case {char(13), char(3)}   % enter and return keys
     % return control to line after waitfor
    set(GETLINE_H1, 'UserData', 'Completed');     
     
    
case 28,
  ROI_callbacks('frm_back'); 
  fast=0;
  
 
case 29,
 ROI_callbacks('frm_forward');
 fast=0;
 
 
    case 30,
        fast=1;
 ROI_callbacks('frm_forward');
 
    case 31,
        fast=1;
        ROI_callbacks('frm_back');
 
 
case 'home',
    
     set(findobj('Tag','frm_sld'),'Value',0);
     set(findobj('Tag','frm_num'),'String',sprintf('Average'));
     axes(findobj('Tag','axesavg')); 
      set(hlv,'Xdata',[1 1],'Ydata',[line_lim(1) line_lim(2)]); 
       set(hlv,'color','k');
     frm_num=0;
     
case 'end',
    
     set(findobj('Tag','frm_sld'),'Value',frames);
     set(findobj('Tag','frm_num'),'String',num2str(frames));
    axes(findobj('Tag','axesavg')); 
            
         axes(findobj('Tag','axesavg')); 
      set(hlv,'Xdata',[frames frames],'Ydata',[line_lim(1) line_lim(2)]); 
        set(hlv,'color','k');
     frm_num=frames;
     axes(findobj('Tag','axesfrm'));
    imshow(img(:,:,frm_num));
    
     


% ADD CONTROLS SO HIGH LOW AND GAMMA CONTROLS CAN BE MOVED !!!! also add
% home to go back fast to avg image

case  32,  %ends marking of rois and saves data
     set(GETLINE_FIG,'Pointer','arrow'); 
        
        
    
        
    
    otherwise, 
    %disp ('THE KEY IS'); 
    key;
    
    
    
 
end

%--------------------------------------------------
% Subfunction FirstButtonDown
%--------------------------------------------------
function FirstButtonDown %#ok

global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_ISCLOSED
global GETLINE_X GETLINE_Y
pause(0.001);
P=ones(16,16)*NaN;
P(8:10,8:10)=1;
P(9,9)=2;



% Set up initial callbacks for initial stage  CHANGE THE TYPE OF CURSOR
set(GETLINE_FIG, ...
    'Pointer','custom','PointerShapeCData',P,'PointerShapeHotSpot',[9 9],...
    'WindowButtonDownFcn', 'getline_adv(''FirstButtonDown'');',...
    'KeyPressFcn', 'getline_adv(''KeyPress'');', ...
    'DoubleBuffer', 'on');

[x,y] = getcurpt(GETLINE_AX);

% check if GETLINE_X,GETLINE_Y is inside of axis
xlim = get(GETLINE_AX,'xlim');
ylim = get(GETLINE_AX,'ylim');
if (x>=xlim(1)) && (x<=xlim(2)) && (y>=ylim(1)) && (y<=ylim(2))
    % inside axis limits
    GETLINE_X = x;
    GETLINE_Y = y;
else
    % outside axis limits, ignore this FirstButtonDown
    return
end

if (GETLINE_ISCLOSED)
    GETLINE_X = [GETLINE_X GETLINE_X];
    GETLINE_Y = [GETLINE_Y GETLINE_Y];
end

set([GETLINE_H1 GETLINE_H2], ...
        'XData', GETLINE_X, ...
        'YData', GETLINE_Y, ...
        'Visible', 'on');
%if it is not left click or double click, then ends (usually mid
if (strcmp(get(GETLINE_FIG, 'SelectionType'), 'extend'))  %normal is left mouse button, Extend is mid button 
    %Alternate, is ctrl left click, and Open is double click of any button
    % We're done!
      set(findobj('Tag','notes'),'visible','on');

    set(GETLINE_H1, 'UserData', 'Completed');
     
 
    
%if there is left click continuos with the drawing of polygon    
elseif (strcmp(get(GETLINE_FIG, 'SelectionType'), 'normal'))
     
    
    % Let the motion functions take over.
    set(GETLINE_FIG, 'WindowButtonMotionFcn', 'getline_adv(''ButtonMotion'');', ...
            'WindowButtonDownFcn', 'getline_adv(''NextButtonDown'');');
end

%--------------------------------------------------
% Subfunction NextButtonDown
%--------------------------------------------------
function NextButtonDown %#ok

global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_ISCLOSED
global GETLINE_X GETLINE_Y

selectionType = get(GETLINE_FIG, 'SelectionType');

if (~strcmp(selectionType, 'open'))
    % We don't want to add a point on the second click
    % of a double-click

    [x,y] = getcurpt(GETLINE_AX);
    if (GETLINE_ISCLOSED)
        GETLINE_X = [GETLINE_X(1:end-1) x GETLINE_X(end)];
        GETLINE_Y = [GETLINE_Y(1:end-1) y GETLINE_Y(end)];
    else
        GETLINE_X = [GETLINE_X x];
        GETLINE_Y = [GETLINE_Y y];
    end
    
    set([GETLINE_H1 GETLINE_H2], 'XData', GETLINE_X, ...
            'YData', GETLINE_Y);
    
end
% only right button closes the polygon.
if (~strcmp(get(GETLINE_FIG, 'SelectionType'), 'normal')) && (~strcmp(get(GETLINE_FIG, 'SelectionType'), 'open'))
    % We're done!
    set(GETLINE_H1, 'UserData', 'Completed');
end

%-------------------------------------------------
% Subfunction ButtonMotion
%-------------------------------------------------
function ButtonMotion %#ok (position of the pointer)

global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_ISCLOSED
global GETLINE_X GETLINE_Y



[newx, newy] = getcurpt(GETLINE_AX);
if (GETLINE_ISCLOSED && (length(GETLINE_X) >= 3))
    x = [GETLINE_X(1:end-1) newx GETLINE_X(end)];
    y = [GETLINE_Y(1:end-1) newy GETLINE_Y(end)];
else
    x = [GETLINE_X newx];
    y = [GETLINE_Y newy];
end

set([GETLINE_H1 GETLINE_H2], 'XData', x, 'YData', y);

%---------------------------------------------------
% Subfunction CleanUp
%--------------------------------------------------
function CleanUp(xlimmode,ylimmode)

xlim(xlimmode);
ylim(ylimmode);
% Clean up the global workspace
clear global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
clear global GETLINE_X GETLINE_Y
clear global GETLINE_ISCLOSED
