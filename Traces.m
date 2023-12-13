% Traces: class to plot multi-channel data
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%


classdef Traces
  
  properties
    data = [];      % [nch,nsamp]
    time = [];      % [1,nsamp]
    label = '';     % traces label
    chlabels = {};   % {nch}
    dataunits = '';
    timeunits = 's';
    prop;
    opt;
  end
  
  methods
    
    %------------------------------------------------------------------
    function t = Traces(varargin)
      %   t = traces(data)
      %   t = traces(data,time)
      
      switch nargin
        
        case 1
          if isstruct(varargin{1}) && isfield(varargin{1},'data') && isfield(varargin{1},'time')
            tt = varargin{1};
            t.data = tt.data;
            t.time = tt.time;
            if isfield(tt,'chlabels') && iscell(tt.chlabels) && length(tt.chlabels)==size(t.data,1)
              for i=1:size(t.data,1)
                if ischar(tt.chlabels{i})
                  t.chlabels{i} = tt.chlabels{i};
                else
                  t.chlabels{i} = sprintf('ch %02i',i);
                end
              end
            end
            
          else
            % assume data matrix
            t.data = varargin{1};
            t.time = [1:size(t.data,2)];
            for i=1:size(t.data,1); t.chlabels{i} = sprintf('ch %02i',i); end
          end
          
        case 2
          % assume data and time
          data = varargin{1};
          t.data = data;
          
          time = varargin{2};
          if max(size(time))==size(t.data,2)
            t.time = time(:)';
          else
            t.time = [1:size(data,2)];
          end
          for i=1:size(t.data,1); t.chlabels{i} = sprintf('ch %02i',i); end
          
        case 3
          % assume data, time, and chlabels
          data = varargin{1};
          time = varargin{2};
          t.data = data;
          t.time = time;
          
          chlabels = varargin{3};
          if iscell(chlabels) && length(chlabels)==size(t.data,1)
            for i=1:size(t.data,1)
              if ischar(chlabels{i})
                t.chlabels{i} = chlabels{i};
              else
                t.chlabels{i} = sprintf('ch %02i',i);
              end
            end
          end
          
      end
      
      t = setdefaults(t);
      
    end
    
    %------------------------------------------------------------------
    function t = setdefaults(t)
      % set default prop and opt
      t.prop.color = 'k';
      t.prop.linewidth = 1.5;
      t.opt.autoscale = 1;
      t.opt.autotrange = 1;
      t.opt.fill = 0;
      t.opt.baseline = 1;
      t.opt.baselinestyle = ':';
      t.opt.box = 0;
      t.opt.yscale = 0;
      t.opt.yscalelabel = '';
      t.opt.ytick = 0;
      
    end
    
    %------------------------------------------------------------------
    function hh = plot(t,h_axes)
      %plot traces object
      
      h = {};
      if nargin<2 | ~ishandle(h_axes) | ~all(strcmp(get(h_axes,'type'),'axes'))
        h_axes = gca;
      end
      
      % plot parameters
      % fillpar = .2; % hsv format -> luminance change
      
      % trace parameters
      ntrace = length(t);
      n = size(t(1).data,1);
      
      if length(h_axes)==ntrace  % one trace per axes
        
        for i=1:ntrace
          axes(h_axes(i))
          taxes = getappdata(h_axes(i),'Traces');
          if isempty(taxes) % first traces being plotted in axes
            setappdata(h_axes(i),'Traces',t(i))
            if t(i).opt.autoscale
              [tmin,tmax] = range(t(i));
              tmin = min(tmin);
              tmax = max(tmax);
            else
              tmin = t(i).opt.ylim(1);
              tmax = t(i).opt.ylim(2);
            end
            dt = tmax-tmin;
            if t(i).opt.autotrange
              [trmin,trmax] = trange(t(i));
              trmin = min(trmin);
              trmax = max(trmax);
            else
              trmin = t(i).opt.xlim(1);
              trmax = t(i).opt.xlim(2);
            end
            
            for j=1:n
              hnew(j,i) = plot(t(i).time,t(i).data(j,:)+dt*(n-j),'linewidth',t(i).prop.linewidth);
              if j==1, hold on, end
              if t(i).opt.fill, hfill(j)=fill([t(i).time(1),t(i).time([1:end]),t(i).time(end)],[0,t(i).data(j,:),0]+dt*(n-j),[.5 .5 .5]); end
            end
            
            set(hnew(:,i),t(i).prop)
            if isfield(t(i).opt,'chcolor') & length(t(i).opt.chcolor)==n;
              for j=1:n
                set(hnew(j,i),'color',t(i).opt.chcolor{j});
                if t(i).opt.fill, set(hfill(j),'facecolor',t(i).opt.chcolor{j}); end
              end
            else
              fillcol = rgb2hsv(get(hnew(1,i),'color'));
              if ~isfield(t(i).opt,'fillpar')
                t(i).opt.fillpar = .2;
              end
              fillcol(3) = mod(fillcol(3)-t(i).opt.fillpar,1);
              if t(i).opt.fill, set(hfill,'facecolor',hsv2rgb(fillcol)), end
            end
            set(h_axes(i),'xlim',[trmin trmax],'ylim',[tmin,dt*(n-1)+tmax],'ytick',dt*[0:n-1],'yticklabel',flipud(t(i).chlabels(:)))
            if t(i).opt.baseline, hl=line(xlim'*ones(1,n),[1;1]*dt*[0:n-1]); set(hl,'linestyle',t(i).opt.baselinestyle,'color',get(hnew(1,i),'color')), end %QUI
            
            xlabel(sprintf('Time [%s]',t(i).timeunits))
            h_tit = title(t(i).label);
            set(h_tit,'color',t(i).prop.color)
            if ~t(i).opt.box, box off, end
            if t(i).opt.yscale  % make yaxis scale
              pos_scale = get(h_axes(i),'Position')*[1 0 0 0;0 1 0 0;1.01 0 .02 0;0 0 0 1];  % x width 5 %
              yl = get(h_axes(i),'YLim');
              h_scale = axes('Position',pos_scale,'Visible','off');
              %         line(trmax*[1 1],tmin+[0 t(i).opt.yscale],'linewidth',1.5,'color','k')
              line(0*[1 1],tmin+[0 t(i).opt.yscale],'linewidth',1.5,'color','k')
              ylim(yl)
              if ~isempty(t(i).opt.yscalelabel), text(.1,tmin+t(i).opt.yscale/2,t(i).opt.yscalelabel,'verticalalignment','middle'), end
              axes(h_axes(i))
            end
            if  t(i).opt.ytick % yticks (additional invisible axes with only ytick)
              axPos = get(h_axes(i),'Position');
              axXlim = get(h_axes(i),'Xlim');
              axYlim = get(h_axes(i),'Ylim');
              for j=1:n
                axPos_j = [axPos(1) axPos(2)+axPos(4)/n*(j-1) axPos(3) axPos(4)/n];
                axYlim_j = [tmin,tmax];
                hh = axes('Position',axPos_j,'Color','none','Visible','off','Xlim',axXlim,'Ylim',axYlim_j);
                if tmin<0
                  ytickmin = ceil((tmin+dt/20)*100)/100;
                  line(axXlim(1)*[1 1]+diff(axXlim)*[0 .02],[1 1]*ytickmin,'color','k','linewidth',1)
                  text(axXlim(1)+diff(axXlim)*.05,ytickmin,num2str(ytickmin))
                end
                if tmax>0
                  ytickmax = floor((tmax-dt/20)*100)/100;
                  line(axXlim(1)*[1 1]+diff(axXlim)*[0 .02],[1 1]*ytickmax,'color','k','linewidth',1)
                  text(axXlim(1)+diff(axXlim)*.05,ytickmax,num2str(ytickmax))
                end
              end
            end
            
          elseif isa(taxes,'Traces') % traces already in axes
            taxes = [taxes, t(i)]; % add axes
            cla(h_axes(i))
            h{i} = plot(taxes,h_axes(i));
            setappdata(h_axes(i),'Traces',taxes)
          else
            warning('wrong app data in axes')
          end
        end
        
        
      else % overlap traces (on first axes)
        
        axes(h_axes(1))
        if t(1).opt.autoscale
          [tmin,tmax] = range(t);  % over all traces
          tmin = min(tmin);
          tmax = max(tmax);
        else
          tmin = t(1).opt.ylim(1);  % use limits of first trace
          tmax = t(1).opt.ylim(2);
        end
        dt = tmax-tmin;
        if t(1).opt.autotrange
          [trmin,trmax] = trange(t);  % over all traces
          trmin = min(trmin);
          trmax = max(trmax);
        else
          trmin = t(1).opt.xlim(1);
          trmax = t(1).opt.xlim(2);
        end
        for i=1:ntrace
          for j=1:n
            hnew(j,i) = plot(t(i).time,t(i).data(j,:)+dt*(n-j),'linewidth',t(i).prop.linewidth); 
            if i==1 & j==1, hold on, end
            if t(i).opt.fill, hfill(j)=fill([t(i).time(1),t(i).time([1:end]),t(i).time(end)],[0,t(i).data(j,:),0]+dt*(n-j),[.5 .5 .5]); end
          end
          set(hnew(:,i),t(i).prop)  
          fillcol = rgb2hsv(get(hnew(1,i),'color')); 
          if ~isfield(t(i).opt,'fillpar')
            t(i).opt.fillpar = .2;
          end
          fillcol(3) = mod(fillcol(3)-t(i).opt.fillpar,1);
          if t(i).opt.fill, set(hfill,'facecolor',hsv2rgb(fillcol)), end
          if t(i).opt.baseline, hl=line(xlim'*ones(1,n),[1;1]*dt*[0:n-1]); set(hl,'linestyle',t(i).opt.baselinestyle,'color',get(hnew(1,i),'color')), end 
          set(h_axes(1),'xlim',[trmin trmax],'ylim',[tmin, dt*(n-1)+tmax])
          if ~isempty(t(i).chlabels), set(h_axes(1),'ytick',dt*[0:n-1],'yticklabel',flipud(t(i).chlabels(:))), end
        end
        xlabel(sprintf('Time [%s]',t(1).timeunits))
        titfmt = repmat('%s ',1,ntrace);
        h_tit = title(sprintf(titfmt(1:end-1),t.label));
        if ~t(1).opt.box, box off, end
        if t(1).opt.yscale  % make yaxis scale
          pos_scale = get(h_axes(1),'Position')*[1 0 0 0;0 1 0 0;1.01 0 .02 0;0 0 0 1];  % x width 5 %
          yl = get(h_axes(1),'YLim');
          h_scale = axes('Position',pos_scale,'Visible','off');
          line(0*[1 1],tmin+[0 t(1).opt.yscale],'linewidth',1.5,'color','k')
          ylim(yl)
          if ~isempty(t(1).opt.yscalelabel), text(.1,tmin+t(1).opt.yscale/2,t(1).opt.yscalelabel,'verticalalignment','middle'), end
        end
        if t(1).opt.ytick % yticks (additional invisible axes with only ytick)
          axPos = get(h_axes(1),'Position');
          axXlim = get(h_axes(1),'Xlim');
          axYlim = get(h_axes(1),'Ylim');
          for j=1:n
            axPos_j = [axPos(1) axPos(2)+axPos(4)/n*(j-1) axPos(3) axPos(4)/n];
            axYlim_j = [tmin,tmax];
            hh = axes('Position',axPos_j,'Color','none','Visible','off','Xlim',axXlim,'Ylim',axYlim_j);
            if tmin<0
              ytickmin = ceil((tmin+dt/20)*100)/100;
              line(axXlim(1)*[1 1]+diff(axXlim)*[0 .02],[1 1]*ytickmin,'color','k','linewidth',1)
              text(axXlim(1)+diff(axXlim)*.05,ytickmin,num2str(ytickmin))
            end
            if tmax>0
              ytickmax = floor((tmax-dt/20)*100)/100;
              line(axXlim(1)*[1 1]+diff(axXlim)*[0 .02],[1 1]*ytickmax,'color','k','linewidth',1)
              text(axXlim(1)+diff(axXlim)*.05,ytickmax,num2str(ytickmax))
            end
          end
        end
        
        
      end
      
      if nargout>0, hh = h; end
      
    end
    
    %------------------------------------------------------------------
    function [tmin,tmax] = range(t)
      %returns min and max over all channels
      
      ntrace = length(t);
      for i=1:ntrace
        tmin(i) = min(min(t(i).data,[],2));
        tmax(i) = max(max(t(i).data,[],2));
      end
      
    end
    
    %------------------------------------------------------------------
    function [trmin,trmax] = trange(t)
      % compute time tange of traces
      
      ntrace = length(t);
      for i=1:ntrace
        trmin(i) = min(t(i).time,[],2);
        trmax(i) = max(t(i).time,[],2);
      end
      
    end
    
    %------------------------------------------------------------------
    
  end
  
  methods (Static)
  end
  
end