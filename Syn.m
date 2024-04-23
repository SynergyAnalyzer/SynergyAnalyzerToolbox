% Syn: class to extract and plot muscle synergies
%
% s = Syn(data,type,inds,nch,chlabels)
%
% list of methods
%------------------------------------------------------------------------
% Syn                   class constructor
% find                  run extraction algorithms
% getDefFindOpt         function to define the defaults of find function
% reconstruct           recostruct EMG data with synergies
% getDefReconstructOpt  function to define the defaults of reconstruct function
% plot                  plot function
% getDefPlotOpt         function to define the defaults of plot function
% indmaxR               find set and repetition with max R^2 for N synergies
% num                   return number of extracted synergies
% numsel                return number of synergies selected with criterion
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%


classdef Syn
    
    properties
        type
        data
        inds
        nch
        W
        C
        R
        Riter
        chlabels
        ind_nonneg
        opt
    end
    
    methods
        %------------------------------------------------------------------
        function obj = Syn(data,type,datapar)
            % create object
            
            if nargin>0
                obj.data = data;
            end
            
            if nargin>1
                obj.type = type;
            else
                % assume spatial type
                obj.type = 'spatial';
            end
            
            if nargin>2
                obj.inds = datapar.inds;
                
                obj.nch = datapar.nch;
                
                if isfield(datapar,'nonnegch')
                    obj.ind_nonneg = datapar.nonnegch;
                end
                
                obj.chlabels = datapar.chlabels;
            end
            
            
        end
        
        %------------------------------------------------------------------
        function obj = find(obj,opt)
            % run extraction algorithm
            
            % set options or use defaults
            defopt = obj.getDefFindOpt(opt.algo);
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            % loop on number of synergies and repetitions
            
            data = obj.data;
            
            % zero-clip negative data
            if opt.clip
                data = data .* (data>0);
            end
            
            % demean data
            datam = mean(data,2);
            datazero = data - datam*ones(1,size(data,2));
            if opt.zeromeandata
                data = datazero;
            end
            
            
            [nch,nsamp] = size(data);
            N = opt.N;
            nset = length(N);
            for iset=1:nset
                for irep=1:opt.nrep
                    if opt.verbose
                        disp(sprintf('synergy extraction: set %i, rep %i',iset,irep))
                    end
                    
                    switch opt.algo
                        
                        case 'mmf' % Scano et al 2020 Mixed-Matrix Factorization algorithm
                            
                            % initialize synergies and coefficients to uniform random
                            % values in [0 1]
                            Wini = rand(nch,N(iset));
                            Cini = rand(N(iset),nsamp);
                            
                            % run algorithm
                            %               [W,C,R]=find_nmf(obj.data,Wini,Cini,opt);
                            [W,C,R]=find_mmf(data,Wini,Cini,struct(opt));
                            
                            if ~opt.bestrsqrep | irep==1
                                obj.W{iset,irep} = W;
                                obj.C{iset,irep} = C;
                                obj.R(iset,irep) = R(end);
                                obj.Riter{iset,irep} = R;
                            elseif R(end) > obj.R(iset,1) % last extraction has higher Rsq
                                obj.W{iset,1} = W;
                                obj.C{iset,1} = C;
                                obj.R(iset,1) = R(end);
                                obj.Riter{iset,1} = R;
                            end
                            
                            
                            
                            
                        case 'nmf' % Lee & Seung NMF algorithm with Eucledian cost
                            
                            % Change this manually
                            if ~isempty(opt.Wini) && size(opt.Wini{N(iset)},2)==N(iset)
                                disp('Using provided Wini')
                                Wini = opt.Wini{N(iset)};
                                Cini = rand(N(iset),nsamp);
                            else
                                % initialize synergies and coefficients to uniform random
                                % values in [0 1]
                                Wini = rand(nch,N(iset));
                                Cini = rand(N(iset),nsamp);
                            end
                            
                            % run algorithm
                            [W,C,R]=find_nmf(obj.data,Wini,Cini,opt);
                            
                            if ~opt.bestrsqrep | irep==1
                                obj.W{iset,irep} = W;
                                obj.C{iset,irep} = C;
                                obj.R(iset,irep) = R(end);
                                obj.Riter{iset,irep} = R;
                            elseif R(end) > obj.R(iset,1) % last extraction has higher Rsq
                                obj.W{iset,1} = W;
                                obj.C{iset,1} = C;
                                obj.R(iset,1) = R(end);
                                obj.Riter{iset,1} = R;
                            end
                            
                        case 'pca' % Principal Component Analysis
                            
                            % run PCA on demeaned data
                            [pc,lat,expl] = pcacov(cov(data'));
                            
                            if opt.zeromeandata
                                W = pc(:,1:N(iset));
                                C = W'*data;
                            else
                                W = [pc(:,1:N(iset)),datam]; % add mean as last synergy
                                C = [pc(:,1:N(iset))'*datazero;ones(1,size(datazero,2))];
                            end
                            
                            R = sum(expl(1:N(iset)))/100;
                            
                            obj.W{iset,irep} = W;
                            obj.C{iset,irep} = C;
                            obj.R(iset,irep) = R;
                            
                    end
                end
                
                obj.opt = opt;
                
            end
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefFindOpt(obj,algo)
            % default extraction algorithm options
            
            if nargin<2, algo = 'nmf'; end
            
            opt.algo = algo;
            opt.N    = 1;
            opt.nrep = 1;
            opt.clip  = 0;
            opt.zeromeandata = 0;
            
            % specific options
            switch algo
                case 'nmf'
                    opt.niter    = [5 5 1e-4 100];  % number of iterations or termination condition and number of max iterations
                    %                     opt.nmaxiter = 100;
                    opt.bestrsqrep  = 1; % 0-> save all reps; 1-> save only best rep
                    opt.plot = 1;
                    opt.updateW = 1;  % 0-> for data fitting
                    opt.clip  = 1;
                    opt.Wini = []; % Change manually, fix W for data fitting
                    opt.print = 0;  % print message at each iteration
                case 'mmf'
                    opt.niter    = [100 10 10^-6 100];  % number of iterations or termination condition and number of max iterations
                    %                     opt.nmaxiter = 100; %10000
                    opt.bestrsqrep  = 1; % 0-> save all reps; 1-> save only best rep
                    opt.plot = 1;
                    opt.updateW = 1;  % 0-> for data fitting
                    opt.print = 0;  % print message at each iteration
                    opt.muC       = 0.1; % >0 -> use gradient descent (clipped to zero for nnW>0)
                    opt.muW       = 0.1; % >0 -> use gradient descent (clipped to zero for nnW>0)
                    opt.lambdaW   = 50; % weight in cost function for suppression of negative synergy components
                    opt.ind_nonneg = obj.ind_nonneg;
                    
                case 'pca'
                    opt.zeromeandata = 1;
            end
            
            opt.verbose = 1;
            
        end
        
        %------------------------------------------------------------------
        function [datahat,inds] = reconstruct(obj,opt)
            % recostruct EMG data with synergies
            
            % set options or use defaults
            defopt = obj.getDefReconstructOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            if isempty(obj.W)
                datahat = zeros(size(obj.data));
                warning('No synergies found')
                return
            end
            
            W = obj.W{opt.iset,opt.irep};
            C = obj.C{opt.iset,opt.irep};
            
            
            switch obj.type
                
                case {'spatial','temporal','spatiotemporal'}
                    ind0 = 0;
                    ind = [];
                    for ii=1:length(opt.isect)
                        i = opt.isect(ii);
                        inds{ii} = ind0 + [1:length(obj.inds{i})];
                        ind0 = ind0 + length(obj.inds{i});
                        ind = [ind obj.inds{i}];
                    end
                    datahat = W*C(:,ind);
                    
            end
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefReconstructOpt(obj)
            % default reconstruction options
            
            opt.isect = 1;
            opt.irep  = 1;
            opt.isect = 1;
            
        end
        
        %------------------------------------------------------------------
        function plot(obj,opt)
            
            % set options or use defaults
            defopt = obj.getDefPlotOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            % activate figure
            figure(opt.h_fig)
            
            % plot
            switch obj.type
                
                case 'spatial'
                    
                    switch opt.type
                        case 'W'
                            W = obj.W{opt.iset,opt.irep};
                            [M,N] = size(W);
                            
                            % sort synergies
                            if isempty(opt.isort)
                                isort = [1:N];
                            else
                                isort = opt.isort(iset,:);
                                W = W(:,isort);
                                N = length(isort);
                            end
                            
                            % prepare axes
                            for i=1:N
                                opt.h_axes(i,1) = axes('Position',ArrayBox(opt.posaxes,1,N,1,i,[.02 .08]));
                            end
                            
                            
                            switch opt.subtype
                                case 'barh'
                                    for i=1:N
                                        labels{i} = sprintf('W_{%i}',isort(i));
                                    end
                                    b = Bars(W,obj.chlabels,labels);
                                    b.opt.type = opt.subtype;
                                    b.opt.color = opt.syncolor;
                                    b.opt.chcolor = opt.chcolor;
                                    plot(b,opt.h_axes);
                                    
                            end
                    end
                    
                case 'temporal'
                    
                    switch opt.type
                        
                        case 'W'
                            C = obj.W{opt.iset,opt.irep};
                            [ntime,N] = size(C);
                            
                            % sort synergies
                            if isempty(opt.isort)
                                isort = [1:N];
                            else
                                isort = opt.isort;
                                C = C(:,isort);
                                N = length(isort);
                            end
                            
                            % limits
                            maxC = max(max(C));
                            minC = min(min(C));
                            for i=1:N
                                % create axes
                                h_axes(i) = axes('Position',ArrayBox(opt.posaxes,N,1,i,1,[.02 .02]));
                                hb = plot(C(:,i),'k-');
                                if ~isempty(opt.syncolor) && iscell(opt.syncolor) && length(opt.syncolor)==N
                                    set(hb,'color',opt.syncolor{i})
                                end
                                box off
                                set(gca,'XLim',[1 ntime],'YLim',[minC maxC])
                                ylabel(sprintf('C_{%i}',isort(i)))
                            end
                    end
                    
                case 'spatiotemporal'
                    
                    switch opt.type
                        
                        case 'W'
                            WW = obj.W{opt.iset,opt.irep};
                            [nchtime,N] = size(WW);
                            M = obj.nch;
                            ntime = nchtime/M;
                            for i=1:N
                                W(:,i,:) = reshape(WW(:,i),M,ntime);
                            end
                            
                            % sort synergies
                            if isempty(opt.isort)
                                isort = [1:N];
                            else
                                isort = opt.isort;
                                W = W(:,isort,:);
                            end
                            
                            % ranges and number of plots
                            maxW = max(W(:));
                            minW = min(W(:));
                            
                            switch opt.subtype
                                
                                case 'color'
                                    for i=1:N
                                        % create axes
                                        h_axes(i) = axes('Position',ArrayBox(opt.posaxes,1,N,1,i,[.02 .08]));
                                        if M==1
                                            plot(squeeze(W(:,i,:)))
                                        else
                                            imagesc(squeeze(W(:,i,:)))
                                        end
                                        xlabel(sprintf('W_{%i}',isort(i)))
                                        
                                        if i==1
                                            set(gca,'ytick',1:size(W,1),'yticklabel',obj.chlabels)
                                        else
                                            set(gca,'ytick',[])
                                        end
                                    end
                                    
                                otherwise
                                    for i=1:N
                                        % create axes
                                        h_axes(i) = axes('Position',ArrayBox(opt.posaxes,1,N,1,i,[.02 .08]));
                                        if M>1
                                            Wi = squeeze(W(:,i,:))';
                                        else
                                            Wi = squeeze(W(:,i,:));
                                        end
                                        tt.data = Wi';
                                        tt.time = [1:size(Wi,1)];
                                        t(i) = Traces(tt);
                                        t(i).opt.autoscale = 0;
                                        t(i).opt.ylim = [minW maxW];
                                        t(i).opt.fill = 1;
                                        if ~isempty(opt.syncolor) && iscell(opt.syncolor) && length(opt.syncolor)==Nplot
                                            t(i).prop.color = opt.syncolor{i};
                                        end
                                        t(i).timeunits = 'samples';
                                        if i==1, t(i).chlabels = obj.chlabels; end
                                        t(i).label = sprintf('W_{%i}',isort(i));
                                        t(i).opt.plottype = 'stack';
                                    end
                                    plot(t,h_axes)
                                    
                            end
                            
                            
                    end
                    
            end
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefPlotOpt(obj)
            % default plot options
            
            opt.type    = 'W';
            opt.iset    = 1;     % first set
            opt.irep    = 1;     % first repetition
            opt.isort   = [];    % synergy sort order
            opt.h_fig   = [];    % figure handle
            opt.h_axes  = [];    % axes handles
            opt.posaxes  = [.02 .02 .96 .96];  % outer axes position
            opt.syncolor = [];
            opt.chcolor  = [];
            
            switch obj.type
                
                case 'spatial'
                    opt.subtype = 'barh';    % horizontal bars
                    
            end
            
        end
        
        %------------------------------------------------------------------
        function [iset,irep] = indmaxR(obj,N)
            % find set and repetition with max R^2 for N synergies
            
            Ns = obj.num;
            
            if isempty(Ns), iset=[]; irep=[]; return, end
            
            if nargin<2, N=Ns(1); end
            
            for i=1:length(N)
                iN = find(N(i)==Ns);
                if ~isempty(iN) && ~isempty(obj.R)
                    iset(i,1) = iN;
                else % find set with closest N
                    [mm,iset(i,1)] = min(abs(Ns-N(i)));
                end
                [rr,irep(i,1)] = max(obj.R(iset(i),:));
            end
        end
        
        %------------------------------------------------------------------
        function N = num(obj,iset)
            % return number of synergies
            
            [nset,nrep] = size(obj.W);    % number of sets and repetitions
            
            if nargin<2, iset=1:nset; end
            
            for i=1:length(iset)
                switch obj.type
                    case {'spatial','temporal','spatiotemporal'}
                        N(i) = size(obj.W{iset(i),1},2);
                end
            end
            
        end
        
        %------------------------------------------------------------------
        function Nsel = numsel(obj,opt)
            % return optimal number of synergies selected according to criteria
            % specified in opt
            %
                        

            % set options or use defaults
            defopt = obj.getDefNumselOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            N = num(obj);
            nset = length(N);

            switch opt.type
              case 'R2thresh'
                [iset,irep] = indmaxR(obj,N);
                for i=1:nset
                  R(i) = obj.R(iset(i),irep(i));
                end
                isel = find(R>=opt.val);
                if not(isempty(isel))
                  Nsel = N(isel(1));
                else
                  Nsel = [];
                  warning('Selection criteria not met!')
                end
              case 'R2fit'
                i0=1;
                iend = nset-1;
                n = N(i0:iend);
                [iset,irep] = indmaxR(obj,N);
                for i=1:nset
                  R(i) = obj.R(iset(i),irep(i));
                end
                for i=i0:iend
                  x = N(i:end);
                  y = R(i:end);
                  p  = polyfit(x,y,1);
                  yhat = polyval(p,x);
                  sse = (y-yhat)*(y-yhat)';
                  mse(i) = sse/length(x);
                end
                isel = find(mse<=opt.val);
                if not(isempty(isel))
                  Nsel = n(isel(1));
                else
                  Nsel = [];
                  warning('Selection criteria not met!')
                end
            end

        end
         %------------------------------------------------------------------
        function opt = getDefNumselOpt(obj)
            % default options to select the number of synergies
            
            % type         val
            % ----------------------
            % 'R2thresh'  value of R^2 threshold
            % 'R2fit'     value of linear fit mse threshold
            % ----------------------
           
            opt.type    = 'R2thresh';   %R2fit
            opt.val    = 0.8;           % 1e-4
            
        end
        
        %------------------------------------------------------------------
        
    end
    
    methods (Static)
    end
    
end

%------------------------------------------------------------------
% subfunctions
%------------------------------------------------------------------
function [W,C,R] = find_nmf(data,W,C,opt)
% extract non-negative synergies with Lee & Seung algorithm

[M,K]  = size(data); % images
N   = size(W,2); % channels x number of synergies


% data matrix V
V = data;

% SST
Vm = mean(V,2);
resid = V - Vm*ones(1,size(V,2));
SST = trace(resid*resid');


% monitor figures
if opt.plot>1
    h_1 = figure('units','normalized','position',[.05 .12 .5 .8]);
end
if opt.plot>0
    h_2 = figure('units','normalized','position',[.6 .12 .3 .4]);
end
errtot = [];
rsq = [];
inderr = [];

% timing info
tic

%
% general loop
%
niter = opt.niter(1);
if length(opt.niter)==4 % [min iter, monitor window, min R^2 diff]
    niter = opt.niter(1);
    iterwin = opt.niter(2);
    errtol = opt.niter(3);
    nmaxiter = opt.niter(4);
else
    iterwin = 1;
    errtol = Inf;
end

% loop while the rsq (abs) difference in the last iterwin iterations is less then errtol
it=0;
while (it<niter | any(abs(diff(rsq(inderr)))>errtol) ) && it<nmaxiter
    it = it+1;
    if opt.print
        disp(sprintf('Iteration %i - time elspsed: %6.4f',it,toc))
        tic
        tlast = toc;
    end
    
    %
    % update C
    %
    den = (W'*W)*C;
    den = den .* (den>0) + (den<=0);
    num = W'*V;
    C = C .* num ./ den;
    
    %
    % update W
    %
    if opt.updateW
        num = V*C';
        den = W*C*C';
        den = den .* (den>0) + (den<=0);
        W = W .* num ./ den;
    end
    
    %
    %display iteration
    %
    if opt.plot>1
        find_nmf_disppar(W,C,h_1)
    end
    
    
    Vhat = W*C;
    
    resid = V-Vhat;
    ee = trace(resid*resid');
    errtot = [errtot ee];
    
    if opt.plot
        find_nmf_disperr(h_2,errtot,SST,errtol);
    end
    inderr = max(1,length(errtot)-iterwin):length(errtot);
    rsq = 1 - errtot/SST;
    
end

R = rsq;

if opt.plot
    close(h_2);
end
if opt.plot>1
    close(h_1);
end

end

%------------------------------------------------------------------
function find_nmf_disppar(W,C,h,labels)
% diplay parameters
figure(h)
N     = size(W,2);
for irow = 1:N
    subplot(N,2,2*(irow-1)+1)
    %imagesc(W(:,irow))
    barh(W(:,irow)), axis ij, axis tight
    title(sprintf('W range= %5.2f %5.2f',min(min(W(:,irow))),max(max(W(:,irow)))))
    if nargin>3
        set(gca,'ytick',1:size(W,1),'yticklabel',labels)
    end
    subplot(N,2,2*(irow-1)+2)
    bar(C(irow,:)), axis tight
end
drawnow
end

%------------------------------------------------------------------
function find_nmf_disperr(h,errtot,sst,errtol)
% diplay error
figure(h)
subplot(2,1,1)
rsq = 1 - errtot/sst;
plot(rsq)
title(sprintf('R^2 = %f',rsq(end)))
h1 = gca;
subplot(2,1,2)
drsq = diff(rsq);
if ~isempty(drsq)
    plot(2:length(rsq),abs(drsq),'r')
    set(gca,'xlim',xlim,'yscale','log')
    title(sprintf('diff R^2 = %e',drsq(end)))
    line(xlim,errtol*[1 1])
end
drawnow
end


%------------------------------------------------------------------
function [W,C,R] = find_mmf(data,Wini,Cini,opt)
%FIND_MMF extract MMF synergies with gradient descent
%
%   [W,C,R] = find_mmf(data,Wini,Cini,opt)
%

[M,K]  = size(data); % images

% if length(size(Wini)) == 3
%     [M,N,mT]   = size(Wini); % channels x number of synergies
%
% elseif length(size(Wini)) == 2
[M,N]   = size(Wini); % channels x number of synergies
% end

opt.lambdaW = opt.lambdaW/(M*N);
opt.muW     = opt.muW/norm(data,'fro');
opt.muC = opt.muW;
% initialize C
if nargin<3 | isempty(Cini)
    C = find_mmf_initC(N,K,1);
else
    C = Cini;
end

% data matrix V
V = data;

% SST
Vm = mean(V,2);
resid = V - Vm*ones(1,size(V,2));
SST = trace(resid*resid');

% W
if length(size(Wini)) > 2
    mT = size(Wini,3);
    for i=1:N
        Wtemp(:,:) = Wini(:,i,:);
        W(:,i) = reshape(Wtemp',1,mT*M); % in case of fit, synergies must be inizialized
    end
else
    W = Wini;
    
end
% W = Wini;

% monitor figures
if opt.plot>1
    h_1 = figure('units','normalized','position',[.05 .12 .5 .8]);
end
if opt.plot>0
    h_2 = figure('units','normalized','position',[.6 .12 .3 .4]);
end
errtot = [];
rsq = [];
inderr = [];

% timing info
tic

%
% general loop
%
niter = opt.niter(1);
if length(opt.niter)==4 % [min iter, monitor window, min R^2 diff]
    niter = opt.niter(1);
    iterwin = opt.niter(2);
    errtol = opt.niter(3);
    nmaxiter = opt.niter(4);
else
    iterwin = 1;
    errtol = Inf;
end

% loop while the rsq (abs) difference in the last iterwin iterations is less then errtol
it=0;
while (it<niter | any(abs(diff(rsq(inderr)))>errtol) ) & it<nmaxiter
    it = it+1;
    if opt.print
        disp(sprintf('Iteration %i - time elspsed: %6.4f',it,toc))
        tic
        tlast = toc;
    end
    
    
    %
    % update W
    %
    if opt.updateW
        
        if opt.muW>0  % gradient descent
            if opt.lambdaW>0
                %           W = W + 2*opt.muW * (data-W*C)*C' - opt.lambdaW * (W<0) .* W; % cost ||W|| only for negative Ws
                W = W + 2*opt.muW*((data-W*C)*C' - opt.lambdaW*W); %correction term added only in spatial optimization
                % enforce non-negativity
                if not(isempty(opt.ind_nonneg)) && (strcmp(opt.algo, 'mmf') )
                    W(opt.ind_nonneg,:) = (W(opt.ind_nonneg,:)>0).*W(opt.ind_nonneg,:);
                end
            else
                W = W + 2*opt.muW*(data-W*C)*C';
            end
        else % least squares
            W = data/C;
        end
        
        
        
    end
    
    %
    % update C
    %
    %         if opt.nnC  % lsqnonneg
    if opt.muC>0  % gradient descent
        
        %                C = C + 2*mu*(W'*(dataset-W*C));
        %         C = (C>0).*C;
        
        C = C + 2*opt.muC*(W'*(data-W*C));
        if not(isempty(opt.ind_nonneg)) && strcmp(opt.type, 'mmf-temporal')
            C(opt.ind_nonneg) = (C(opt.ind_nonneg)>0).*C(opt.ind_nonneg);
            
        else
            C = C.*(C>0); % clip to zero
        end
    else
        for i=1:K
            C(:,i) = lsqnonneg(W,data(:,i));
        end
    end
    %         else
    %             if opt.muC >0 % gradient descent
    %                 C = C + 2*opt.muC*W'*(data-W*C);
    %                 %               C = C.*(C>0); % clip to zero
    %             else
    %                 % least-squares
    %                 C = W\data;
    %             end
    %         end
    
    
    %
    %display iteration
    %
    if opt.plot>1
        find_mmf_disppar(W,C,h_1)
    end
    
    
    Vhat = W*C;
    errtot = find_mmf_err(errtot,SST,V,Vhat);
    if opt.plot
        find_mmf_disperr(h_2,errtot,SST,errtol);
    end
    inderr = max(1,length(errtot)-iterwin):length(errtot);
    rsq = 1 - errtot/SST;
    
end


% if opt.normW
%       normW = sqrt(ones(1,M)*W.^2);
%       W = W ./ (ones(M,1)*normW);
for in = 1:N
    normW = norm(W(:,in)); %%% TO BE CHECKED
    Wnorm = W(:,in)/normW;
    Cnorm = C(in,:)*normW;
    
    W(:,in) = Wnorm;
    C(in,:) = Cnorm;
    
end
%      for i=1:order
%         W_new_norm(:,i) = W_new(:,i)/norm(W_new(:,i));
%         norms_syn = norm(W_new(:,i));
%         C_new_norm(i,:) = C_new(i,:) * norms_syn;
%      end

% end


R = rsq;

if opt.plot
    close(h_2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = find_mmf_initC(N,K,initype);
% initialize C
switch initype
    case 1 % random in [0,1]
        C = rand(N,K);
end
end

function find_mmf_disppar(W,C,h,labels)
% diplay parameters
figure(h)
N     = size(W,2);
for irow = 1:N
    subplot(N,2,2*(irow-1)+1)
    %imagesc(W(:,irow))
    barh(W(:,irow)), axis ij, axis tight
    title(sprintf('W range= %5.2f %5.2f',min(min(W(:,irow))),max(max(W(:,irow)))))
    if nargin>3
        set(gca,'ytick',1:size(W,1),'yticklabel',labels)
    end
    subplot(N,2,2*(irow-1)+2)
    bar(C(irow,:)), axis tight
end
drawnow
end

function errtot = find_mmf_err(errtot,sst,V,Vhat)
% compute error
resid = V-Vhat;
ee = trace(resid*resid');
errtot = [errtot ee];
end

function find_mmf_disperr(h,errtot,sst,errtol)
% diplay error
figure(h)
subplot(2,1,1)
rsq = 1 - errtot/sst;
plot(rsq)
title(sprintf('R^2 = %f',rsq(end)))
h1 = gca;
subplot(2,1,2)
drsq = diff(rsq);
if ~isempty(drsq)
    plot(2:length(rsq),abs(drsq),'r')
    set(gca,'xlim',xlim,'yscale','log')
    title(sprintf('diff R^2 = %e',drsq(end)))
    line(xlim,errtol*[1 1])
end
drawnow
end
