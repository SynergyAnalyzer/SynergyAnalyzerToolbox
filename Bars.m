% Bars: class to plot multiple bar plots
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%

classdef Bars
  
  properties
    data = [];      % [nch,n]
    chlabels = {};   % {nch}
    labels = {};  % {n}
    dataunits = '';
    prop;
    opt;
  end
  
  methods
    
    %------------------------------------------------------------------
    function b = Bars(varargin)
      %   b = traces(data)
      %   b = traces(data,chlabels)
      %   b = traces(data,chlabels,labels)
      
      if nargin>0
          % assume data matrix
          b.data = varargin{1};
          for i=1:size(b.data,1);  b.chlabels{i} = sprintf('ch %02i',i); end
          for i=1:size(b.data,2);  b.labels{i} = sprintf('%02i',i); end
      end
      
      if nargin>1
          % assume chlabels
          chlabels = varargin{2};
          if iscell(chlabels) && length(chlabels)==size(b.data,1)
            for i=1:size(b.data,1)
              if ischar(chlabels{i})
                b.chlabels{i} = chlabels{i};
              end
            end
          end
      end
         
      if nargin>2
          % assume labels
          labels = varargin{3};
          if iscell(labels) && length(labels)==size(b.data,2)
            for i=1:size(b.data,2)
              if ischar(labels{i})
                b.labels{i} = labels{i};
              end
            end
          end
          
      end
      
      b = setdefaults(b);
      
    end
    
    %------------------------------------------------------------------
    function b = setdefaults(b)
      % set default prop and opt
      b.prop.color = 'k';
      
      b.opt.type = 'barh'; % 'bar'
      b.opt.color = {};    % column colors
      b.opt.chcolor = {};    % channel colors
      b.opt.posaxes  = [.04 .02 .94 .96];  % outer axes position

      
    end
    
    %------------------------------------------------------------------
    function h_axes = plot(b,h_axes)
      %plot bars object
      
      nset = length(b);
      
      for iset=1:nset
        
        W = b(iset).data;
        [M,N] = size(W);
        
        if nargin<2 || ~all(ishandle(h_axes(:))) || ~all(strcmp(get(h_axes,'type'),'axes')) || ~isequal(size(h_axes,1),N) || ~isequal(size(h_axes,2),nset)
          % create figure
          if iset==1, h_fig = figure; end
          for i=1:N
            % create axes
            switch b.opt.type
              case 'barh'
                h_axes(i,iset) = axes('Position',ArrayBox(b.opt.posaxes,nset,N,iset,i,[.02 .08]));
              case 'bar'
                h_axes(i,iset) = axes('Position',ArrayBox(b.opt.posaxes,N,nset,i,iset,[.02 .08]));
            end
          end
        end
        
        % limits
        maxW = max(max(W));
        minW = min(min(W));
        
        switch b(iset).opt.type
          case 'barh'
            
            for i=1:N
              axes(h_axes(i,iset))
              for k=1:M
                hb(k) = barh(k,W(k,i));
                hold on
              end
              if ~isempty(b(iset).opt.color)&& length(b(iset).opt.color)==N && ~isempty(b(iset).opt.chcolor) && length(b(iset).opt.chcolor)==M
                 for k=1:M
                  set(hb(k),'facecolor',b(iset).opt.chcolor{k}*b(iset).opt.color{i})
                end 
              elseif ~isempty(b(iset).opt.color) && iscell(b(iset).opt.color) && length(b(iset).opt.color)==N
                set(hb,'facecolor',b(iset).opt.color{i})
              elseif ~isempty(b(iset).opt.chcolor) && iscell(b(iset).opt.chcolor) && length(b(iset).opt.chcolor)==M
                for k=1:M
                  set(hb(k),'facecolor',b(iset).opt.chcolor{k})
                end
              end
              axis('ij')
              box off
              set(gca,'XLim',[min(0,minW*1.1) maxW*1.1],'YLim',[0 M+1],'YTick',[1:M],'Yticklabel','')
              if i==1
                set(gca,'YTickLabel',b(iset).chlabels)
              end
              xlabel(b(iset).labels{i})
            end
            
        end
      end
      
    end
   
    %------------------------------------------------------------------
    
  end
  
  methods (Static)
  end
  
end