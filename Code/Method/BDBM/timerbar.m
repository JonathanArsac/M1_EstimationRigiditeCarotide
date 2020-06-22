function fout = timerbar(x,whichbar)
%TIMERBAR Display timer bar.
%   H = TIMERBAR(X,'title') creates and displays a wait bar of 
%   fractional length X.  The handle to the timerbar figure is
%   returned in H.  X should be between 0 and 1. It also displays 
%   the estimated time left for computation
%
%   TIMERBAR(X,H) will set the length of the bar in waitbar H
%   to the fractional length X.
%
%   TIMERBAR(X) will set the length of the bar in the most recently
%   created waitbar window to the fractional length X.
%
%   TIMERBAR is typically used inside a FOR loop that performs a 
%   lengthy computation.  A sample usage is shown below:
%
%       h = timerbar(0,'Please wait...');
%       for i=1:100,
%           % computation here %
%           timerbar(i/100,h)
%       end
%       close(h)
%
%   See also WAITBAR

%   Suresh E Joel 06-03-2002
%   Copyright (c) 2002 by Suresh Joel
%

if nargin==2 & ischar(whichbar)
   type=2; % we are initializing
   name=whichbar;
   t0=clock;
elseif nargin==2 & ishandle(whichbar)
   type=1; % we are updating, given a handle
   f=whichbar;
elseif nargin==1
   f = findobj(allchild(0),'flat','Tag','TMWTimerbar');
   
   if isempty(f)
      type=2;
      name='TimerBar';
   else
      type=1;
      f=f(1);
   end   
elseif nargin==2
   error('Two-argument syntax requires TIMERBAR(X,''title'') or TIMERBAR(X,H)')
else
   error('Input arguments not valid.');
end

x = max(0,min(100*x,100));

switch type
case 1,  % waitbar(x)    update
   p = findobj(f,'Type','patch');
   l = findobj(f,'Type','line');
   a = findobj(f,'Tag','Timebx');
   if isempty(f) | isempty(p) | isempty(l) | isempty(a),
      error('Couldn''t find timerbar handles.'); 
   end
   t0= get(f,'UserData');
   xpatch = get(p,'XData');
   %replace 0 with xpatch(2)
   xpatch = [xpatch(2) x x xpatch(2)];
   set(p,'XData',xpatch')
   xline = get(l,'XData');
   set(l,'XData',xline);
   t=etime(clock,t0);
   if(x~=0),t=(t/x)*(100-x);else return; end;
   t=round(t);
   if t>60 & t<3600,
      str = strcat('Time left :',num2str(fix(t/60)),'mins,',num2str(round(60*rem(t/60,1))),'secs');
   elseif t>3600
      str = strcat('Time left :',num2str(fix(t/3600)),'hrs,',num2str(round(60*rem((t/3600),1))),'mins');
   else
      str = strcat('Time left :',num2str(t),'secs');
   end
   set(a,'String',str);

case 2,  % waitbar(x,name)  initialize

   oldRootUnits = get(0,'Units');

   set(0, 'Units', 'points');
   screenSize = get(0,'ScreenSize');
   
   axFontSize=get(0,'FactoryAxesFontSize');
   
   pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');
      
   width = 360 * pointsPerPixel;
   height = 75 * pointsPerPixel;
   pos = [screenSize(3)/2-width/2 screenSize(4)/2-height/2 width height];

   f = figure(...
           'Units', 'points', ...
           'Position', pos, ...
           'Resize','off', ...
           'CreateFcn','', ...
           'NumberTitle','off', ...
           'IntegerHandle','off', ...
           'MenuBar', 'none', ...
           'Tag','TMWTimerbar',...
           'Visible','off',...
           'UserData',t0);
   colormap([]);
   
   axNorm=[.05 .3 .9 .2];
   axPos=axNorm.*[pos(3:4),pos(3:4)];
   
   tbox=uicontrol('Parent',f,...
      'Units','points',...
      'Position',[pos(3)/4 1 pos(3)/2 pos(4)/4],...
      'String','Time left : ...',...
      'Style','text',...
      'Tag','Timebx');
   
   h = axes('XLim',[0 100],...
      'YLim',[0 1],...
      'Box','on', ...
      'Units','Points',...
      'FontSize', axFontSize,...
      'Position',axPos,...
      'XTickMode','manual',...
      'YTickMode','manual',...
      'XTick',[],...
      'YTick',[],...
      'XTickLabelMode','manual',...
      'XTickLabel',[],...
      'YTickLabelMode','manual',...
      'YTickLabel',[]);
   
   tHandle=title(name);
   tHandle=get(h,'title');
   oldTitleUnits=get(tHandle,'Units');
   set(tHandle,...
      'Units',      'points',...
      'String',     name);
   
   tExtent=get(tHandle,'Extent');
   set(tHandle,'Units',oldTitleUnits);
   
   titleHeight=tExtent(4)+axPos(2)+axPos(4)+5;
   if titleHeight>pos(4)
      pos(4)=titleHeight;
      pos(2)=screenSize(4)/2-pos(4)/2;
      figPosDirty=logical(1);
   else
      figPosDirty=logical(0);
   end
   
   if tExtent(3)>pos(3)*1.10;
      pos(3)=min(tExtent(3)*1.10,screenSize(3));
      pos(1)=screenSize(3)/2-pos(3)/2;
      
      axPos([1,3])=axNorm([1,3])*pos(3);
      set(h,'Position',axPos);
      
      figPosDirty=logical(1);
   end
   
   if figPosDirty
      set(f,'Position',pos);
   end

   xpatch = [0 x x 0];
   ypatch = [0 0 1 1];
   xline = [100 0 0 100 100];
   yline = [0 0 1 1 0];
   
   p = patch(xpatch,ypatch,'r','EdgeColor','r','EraseMode','none');
   l = line(xline,yline,'EraseMode','none');
   set(l,'Color',get(gca,'XColor'));

   set(f,'HandleVisibility','callback','visible','on');
   
   set(0, 'Units', oldRootUnits);
end  % case
drawnow;

if nargout==1,
  fout = f;
end