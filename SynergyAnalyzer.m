% SynergyAnalyzer: class to extract muscle synergies from EMG/kin data
%
% obj = SynergyAnalyzer(data,info,par)
%
% list of methods
%-----------------------------------------------
% dataFilter    filter data
% normalize     normalize EMG and kin data amplitude
% emgPhasic     subtract tonic activity (based on linear ramp from onset to end) to extract phasic activity
% zeromean      subtract mean value of each emg channel in each trial
% average       average data across trials
% groupTrials   group trials according to experimental conditions
% findTrials    find trial number according to conditions
% find          run synergy extraction algorithm
% plot          plot synergies, coefficient, and EMG data reconstruction
% reconstruct   synergy reconstruct of data
% coeftrace     synergy coefficient traces
% nch           number of channels in data matrix
% nonnegch      number of non-negative channels in data matrix
%
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%


classdef SynergyAnalyzer
    
    properties
        data
        syn
        isect
        info
        opt
    end
    
    methods
        %------------------------------------------------------------------
        function obj = SynergyAnalyzer(data,info,par)
            % creates object
            
            defpar = getDefPar;
            % set options or use defaults
            if nargin>2 && isstruct(par)
                fname = fieldnames(par);
                for i=1:length(fname)
                    defpar.(fname{i}) = par.(fname{i});
                end
            end
            par = defpar;
            
            
            type = par.type;
            delay = par.delay;
            order = par.order;
            chlabels = par.chlabels;
            
            
            if nargin <2
                warning('Data Info is missing, please provide info - see getInfo.m function')
                return
            else
                obj.info = info;
            end
            
            
            switch type
                case  'EmgData'
                    obj.data = EmgData(data,[],chlabels);
                    
                case 'KinData'
                    obj.data = KinData(data,[],chlabels);
                case 'EmgKinData'
                    
                    nchemg = size(data(1).emg,1);
                    nchkin = size(data(1).pos,1);
                    if not(isempty(chlabels)) && length(chlabels) == nchemg+nchkin
                        emg = EmgData(data,[],chlabels(1:nchemg));
                        kin = KinData(data,[],chlabels(nchemg+1:nchkin+nchemg));
                        obj.data = EmgKinData(emg,kin,[],order,delay);
                    elseif not(isempty(chlabels)) && length(chlabels) == nchemg
                        emg = EmgData(data,[],chlabels);
                        kin = KinData(data,[]);
                        obj.data = EmgKinData(emg,kin,[],order,delay);
                    elseif not(isempty(chlabels)) && length(chlabels) == nchkin
                        emg = EmgData(data,[]);
                        kin = KinData(data,[],chlabels);
                        obj.data = EmgKinData(emg,kin,[],order,delay);
                    else
                        emg = EmgData(data,[]);
                        kin = KinData(data,[]);
                        obj.data = EmgKinData(emg,kin,[],order,delay);
                    end
            end
            
            obj.syn = Syn;
            
            obj.opt = getDefOpt;
            
        end
        
        %------------------------------------------------------------------
        function obj = dataFilter(obj)
            
            switch class(obj.data)
                case 'EmgData'
                    % filter data
                    if obj.opt.verbose
                        disp(sprintf('filtering EMG data using %s filter',obj.opt.emgFilter.type))
                    end
                    obj.data = filter(obj.data,obj.opt.emgFilter);
                case 'EmgKinData'
                    % filter data
                    if obj.opt.verbose
                        disp(sprintf('filtering EMG data using %s filter',obj.opt.emgFilter.type))
                        disp(sprintf('filtering KIN data using %s filter',obj.opt.kinFilter.type))
                    end
                    obj.data = filter(obj.data,obj.opt);
                    
                case 'KinData'
                    % filter data
                    if obj.opt.verbose
                        disp(sprintf('filtering KIN data using %s filter',obj.opt.kinFilter.type))
                    end
                    obj.data = filter(obj.data,obj.opt.kinFilter);
            end
                        
            if obj.opt.emgFilter.resample
                if obj.opt.verbose
                    disp(sprintf('resampling EMG and KIN data using %5.3f [s] period',obj.opt.emgFilter.resample_period))
                end
                obj.data = resample(obj.data,obj.opt.emgFilter.resample_period);
            end
            
        end
        
        %------------------------------------------------------------------
        function obj = normalize(obj)
            % normalize EMG and kin data amplitude
            if isempty(obj.opt.normalize.isect)
                obj.opt.normalize.isect = [1:length(obj.data)];
            end
            obj.data = normalize(obj.data,obj.opt.normalize);
            
        end
        
        %------------------------------------------------------------------
        function obj = emgPhasic(obj)
            % subtract tonic activity (based on linear ramp from onset to end) to
            % extract phasic activity
            
            if ~isfield(obj.info,'events')
                warning('Movement onset and end events required for phasic EMGs')
                return
            end
            
            opt.type = 'tonic';
            
            nsect = length(obj.info);
            for i=1:nsect
                ind = find(obj.info(i).events.code==13); % onset
                opt.t_onset(i) = obj.info(i).events.time(ind);
                ind = find(obj.info(i).events.code==14); % offset
                opt.t_end(i) = obj.info(i).events.time(ind);
            end
            
            obj.data = subtract(obj.data,opt);
            
        end
        
        %------------------------------------------------------------------
        function obj = zeromean(obj)
            % subtract mean value of each emg channel in each trial
            opt.type = 'mean';
            obj.data = subtract(obj.data,opt);
            
        end
        
        %------------------------------------------------------------------
        function obj = average(obj)
            
            % trial selection
            nsect = length(obj.info);
            if isempty(obj.opt.average.gr)
                if isempty(obj.opt.average.grcond)
                    gr{1}=[1:nsect]; % a single group with all trials
                else
                    gr = obj.groupTrial(opt.grcond{:});
                end
            else
                gr = obj.opt.average.gr;
            end
            
            isect = unique([gr{:}]); % all included sections
            if isempty(isect)
                warning('No groups selected!')
                return
            end
            if obj.opt.verbose
                disp(sprintf('averaging trials using %i groups',length(gr)))
            end
            
            % tref
            tref = zeros(nsect,1);  % assumed data have been aligned
            
            % trange
            trange = obj.opt.average.trange;
            
            % average emg data over trials
            opt.gr = gr;
            opt.tref = tref;
            opt.trange = trange;
            obj.data = average(obj.data,opt);
            
            % update info
            ngr = length(opt.gr);
            for i=1:ngr
                info(i).id = i;
                info(i).type = i;
                info(i).selected = 1;
                if isfield(obj.info(i),'nonnegch')
                    info(i).nonnegch = obj.info(i).nonnegch;
                end
            end
            obj.info = info;
            
        end
        
        %------------------------------------------------------------------
        function gr = groupTrials(obj,varargin)
            % group trial number according to conditions
            %
            %   gr = groupTrials(obj,condtype,condval[,condtype2,condval,...])
            %
            %   gr -> cell array {ngr}
            %
            %   each condition must have exactely ngr componets
            %
            %   contype     condval        notes
            %   ----------------------------------------------------------
            %   'type'     [ngr,ntype]     matches the info.type vector
            %   'type1'    [ngr,1]         matches the info.type(1) value
            %   'type2'    [ngr,1]         matches the info.type(2) value
            %   'type3'    [ngr,1]         matches the info.type(3) value
            %   'typei'    [ngr,2] (i,val) matches the info.type(i) value
            %   'selected' [ngr,1]         matches the info.selected value
            %   'id'       {ngr}           matches info.id value
            
            if nargin<2 | ~isequal(rem(length(varargin),2),0)
                ncond = 0;
            else
                ncond = length(varargin)/2;
                condtype = varargin(1:2:ncond*2);
                condval  = varargin(2:2:ncond*2);
            end
            
            ninfo = length(obj.info);
            
            ngr = size(condval{1},1);
            
            for k=1:ngr
                gr{k} = [];
                for i=1:ninfo
                    ind_j = [];
                    for j=1:ncond
                        nclass = size(condval{j},1);
                        switch condtype{j}
                            case 'type'
                                if isequal(obj.info(i).type,condval{j}(k,:)), ind_j = [ind_j j]; end
                            case 'type1'
                                if isequal(obj.info(i).type(1),condval{j}(k)), ind_j = [ind_j j]; end
                            case 'type2'
                                if isequal(obj.info(i).type(2),condval{j}(k)), ind_j = [ind_j j]; end
                            case 'type3'
                                if isequal(obj.info(i).type(3),condval{j}(k)), ind_j = [ind_j j]; end
                            case 'typei'
                                if isequal(obj.info(i).type(condval{j}(k,1)),condval{j}(k,2)), ind_j = [ind_j j]; end
                            case 'ok'
                                if isequal(obj.info(i).ok,condval{j}(k)), ind_j = [ind_j j]; end
                            case 'selected'
                                if isequal(obj.info(i).selected,condval{j}(k)), ind_j = [ind_j j]; end
                            case 'id'
                                if ismember(obj.info(i).id,condval{j}{k}), ind_j = [ind_j j]; end
                        end % switch
                    end %j
                    if isequal(length(ind_j),ncond), gr{k} = [gr{k} i]; end
                end % i
            end % k
            
        end
        
        %------------------------------------------------------------------
        function ind = findTrials(obj,varargin)
            % find trial number according to conditions
            %
            %   ind = findTrials(obj,condtype,condval[,condtype2,condval,...])
            %
            %   contype     condval        notes
            %   ----------------------------------------------------------
            %   'id'       [nsect,1]       matches the info.id
            %   'type'     [nclass,ntype]  matches the info.type vector
            %   'type1'    [nclass,1]      matches the info.type(1) value
            %   'type2'    [nclass,1]      matches the info.type(2) value
            %   'type3'    [nclass,1]      matches the info.type(3) value
            %   'selected' [nclass,1]      matches the info.selected value
            %
            %    logical OR between classes for each condition
            %    logical AND between conditions
            
            if nargin<2 | ~isequal(rem(length(varargin),2),0)
                ncond = 0;
            else
                ncond = length(varargin)/2;
                condtype = varargin(1:2:ncond*2);
                condval  = varargin(2:2:ncond*2);
            end
            
            ntrial = length(obj.info);
            ind = [];
            for i=1:ntrial
                ind_j = [];
                for j=1:ncond
                    switch condtype{j}
                        case 'id'
                            if ismember(obj.info(i).id,condval{j}), ind_j = [ind_j j]; end
                        case 'type'
                            nclass = size(condval{j},1);
                            for k=1:nclass
                                if isequal(obj.info(i).type,condval{j}(k,:)), ind_j = [ind_j j]; break, end
                            end
                        case 'type1'
                            nclass = size(condval{j},1);
                            for k=1:nclass
                                if isequal(obj.info(i).type(1),condval{j}(k)), ind_j = [ind_j j]; break, end
                            end
                        case 'type2'
                            nclass = size(condval{j},1);
                            for k=1:nclass
                                if isequal(obj.info(i).type(2),condval{j}(k)), ind_j = [ind_j j]; break, end
                            end
                        case 'type3'
                            nclass = size(condval{j},1);
                            for k=1:nclass
                                if isequal(obj.info(i).type(3),condval{j}(k)), ind_j = [ind_j j]; break, end
                            end
                        case 'typei'
                            ii = condval{j}{1};
                            classes = condval{j}{2};
                            nclass = size(classes,1);
                            for k=1:nclass
                                if isequal(obj.info(i).type(ii),classes(k)), ind_j = [ind_j j]; break, end
                            end
                        case 'selected'
                            if isequal(obj.info(i).selected,condval{j}), ind_j = [ind_j j]; end
                            
                    end % switch
                end % j
                if isequal(length(ind_j),ncond), ind = [ind i]; end
            end % i
            
        end
        
        %------------------------------------------------------------------
        function obj = find(obj)
            % run synergy extraction algorithm
            
            % select sections
            if isempty(obj.isect)
                obj.isect = obj.findTrials('selected',1);
            end
            
            % prepare data matrix
            type = obj.opt.find.type;
            [ddata,datapar] = obj.data.getData(type,obj.isect);
            
            
            % create Syn object
            obj.syn = Syn(ddata,type,datapar);
            
            % run algorithm
            if obj.opt.verbose
                disp(sprintf('extracting %s synergies using %s alogrithm',obj.opt.find.type,obj.opt.find.algo))
            end
            
            obj.syn = find(obj.syn,obj.opt.find);
            
            
        end
        
        %------------------------------------------------------------------
        function plot(obj)
            % plot synergies, coefficient, and EMG data reconstruction
            
            if isempty(obj.syn), return, end
            
            opt = obj.opt.plot;
            
            switch opt.type
                
                case 'rsq'
                    % plot R^2 vs. N
                    
                    % prepare figure
                    figname = sprintf('R^2 vs. N for %s synergies (%s algorithm)',obj.syn.opt.type,obj.syn.opt.algo);
                    opt.h_fig = figure(...
                        'NumberTitle','off',...
                        'Name',figname,...
                        'PaperType',opt.papertype );
                    
                    plot(num(obj.syn),obj.syn.R,'k.-')
                    xlabel('Number of synergies')
                    ylabel('R^2')
                    
                case 'W'
                    % plot one set of synergies in one figure
                    
                    % select repetition with N elements with max R^2 (if opt.N is specified)
                    if ~isempty(opt.N)
                        [opt.iset,opt.irep] = obj.syn.indmaxR(opt.N);
                    end
                    N = num(obj.syn,opt.iset);
                    
                    % specific options
                    opt.posaxes = [.08 .05 .9 .9]; % axes box
                    
                    % prepare figure
                    figname = sprintf('%s synergies (%s algorithm)',obj.syn.opt.type,obj.syn.opt.algo);
                    opt.h_fig = figure(...
                        'NumberTitle','off',...
                        'Name',figname,...
                        'PaperType',opt.papertype);
                    
                    % call syn plotting method
                    plot(obj.syn,opt);
                    
                case 'rec'
                    % plot the original data and the reconstruction as synerergy combinations for opt.isect sections
                    
                    % prepare figure
                    figname = sprintf('%s syn. reconstruction (%s algorithm)',obj.syn.opt.type,obj.syn.opt.algo);
                    opt.h_fig = figure(...
                        'NumberTitle','off',...
                        'Name',figname,...
                        'PaperType',opt.papertype,...
                        'Position', [440 377 560 620],...
                        'PaperOrientation','Landscape');
                    
                    % set specific options
                    opt.subtype = '';
                    
                    if isa(obj.data,'EmgKinData')
                        opt.posaxes = [.08 .55 .9 .40; .08 .30 .9 .20]; % 3 axes box (2 for traces and 1 for coefs)
                        coefaxes = [.08 .05 .9 .20];
                    else
                        
                        opt.posaxes = [.08 .40 .9 .55]; % 2 axes box (1 for traces and 1 for coefs)
                        coefaxes = [.08 .05 .9 .30];
                        
                    end
                    
                    opt.autoscale = 0; % autoscale traces ylim
                    
                    % select repetition with N elements with max R^2 (if opt.N is specified)
                    if ~isempty(opt.N)
                        [opt.iset,opt.irep] = indmaxR(obj.syn,opt.N);
                    end
                    % select specific sections
                    if isempty(opt.isect)
                        opt.isect = obj.isect;
                    end
                    % select sections used for extraction
                    opt.isect = intersect(opt.isect,obj.isect);
                    if isempty(opt.isect)
                        warning('Specified sections not available')
                        return
                    end
                    nsect = length(opt.isect);
                    
                    % compute reconstruction of selected sections
                    data = obj.data(opt.isect);
                    datarec = reconstruct(obj,opt);
                    dur = duration(data);
                    
                    % plot data and reconstruction
                    optdata.color = 'k';
                    optdata.fill = 1;
                    optdata.figure = opt.h_fig;
                    for i=1:nsect
                        pos = SizBox(opt.posaxes(1,:),1,dur,1,i,[.02 .02]);
                        optdata.axes(1,i) = axes('Position',pos);
                    end
                    if isa(obj.data,'EmgKinData')
                        for i=1:nsect
                            pos = SizBox(opt.posaxes(2,:),1,dur,1,i,[.02 .02]);
                            optdata.axes(2,i) = axes('Position',pos);
                        end
                    end
                    ha = plot(data,optdata);
                    
                    optdata.color = 'k';
                    optdata.linewidth = 2;
                    optdata.fill = 0;
                    optdata.usetitle = 0;
                    optdata.figure = opt.h_fig;
                    optdata.axes = ha;
                    optdata.rec = 1;
                    plot(datarec,optdata);
                    
                    % compute coefficeint time courses
                    tcoef = coeftrace(obj,opt);
                    
                    % plot coefs on second axes
                    switch obj.syn.type
                        
                        case 'spatial'
                            [tcmin,tcmax] = range(tcoef);
                            tcylim = [min(tcmin) max(tcmax)]; % if autoscale is off
                            for i=1:nsect
                                tcoef(i).opt.fill = 1;
                                tcoef(i).opt.autoscale = opt.autoscale;
                                tcoef(i).opt.ylim = tcylim;
                                tcoef(i).opt.profile = 1;
                                pos = SizBox(coefaxes,1,dur,1,i,[.02 .02]);
                                h_axes_2(i) = axes('Position',pos);
                            end
                            plot(tcoef,h_axes_2);
                            
                        case 'temporal'
                            N = size(obj.syn.W{opt.iset,opt.irep},2);
                            for i=1:nsect
                                pos = SizBox(coefaxes,1,dur,1,i,[.02 .02]);
                                for j=1:N
                                    h_axes_2(j,i) = axes('Position',ArrayBox(pos,1,N,1,j,[.02 .08]));
                                end
                            end
                            plot(tcoef,h_axes_2);
                            
                        case 'spatiotemporal'
                            N = size(obj.syn.W{opt.iset,opt.irep},2);
                            for i=1:nsect
                                h_axes_2(i) = axes('Position',SizBox(coefaxes,1,dur,1,i,[.02 .02]));
                            end
                            plot(tcoef,h_axes_2);
                            
                    end
                    
            end
            
        end
        
        %------------------------------------------------------------------
        function datarec = reconstruct(obj,opt)
            % synergy reconstruct of data
            
            isect = opt.isect;
            nsect = length(isect);
            
            [datahat,inds] = reconstruct(obj.syn,opt);
            
            for ii=1:nsect
                datarec(ii) = obj.data(isect(ii));
                % reshape datahat if necessary
                switch obj.syn.type
                    case 'spatial'
                        datarec(ii).data = datahat(:,inds{ii});
                    case 'temporal'
                        datarec(ii).data = datahat(:,inds{ii})';
                    case 'spatiotemporal'
                        nch = size(datarec(ii).data,1);
                        ntime = size(datahat,1)/nch;
                        datarec(ii).data = reshape(datahat(:,inds{ii}),nch,ntime);
                end
            end
            
        end
        
        %------------------------------------------------------------------
        function t = coeftrace(obj,opt)
            % synergy coefficient traces
            
            isect = opt.isect;
            nsect = length(isect);
            
            C = obj.syn.C{opt.iset,opt.irep};
            N = size(C,1);
            
            % sort synergies
            if isempty(opt.isort) || ~isequal(length(opt.isort),N)
                isort = [1:N];
            else
                isort = opt.isort;
                C = C(isort,:);
            end
            
            for ii=1:nsect
                i = opt.isect(ii);
                tt = [];
                
                switch obj.syn.type
                    
                    case 'spatial'
                        tt.data = C(:,obj.syn.inds{i});
                        tt.time = obj.data(i).time;
                        for j=1:N, chstr{j} = sprintf('C_%i',isort(j)); end
                        if ii==1, tt.chlabels = chstr; end
                        t(ii) = Traces(tt);
                        if ~isempty(opt.syncolor) && length(opt.syncolor)==N
                            t(ii).opt.chcolor = opt.syncolor;
                        end
                        
                    case 'temporal'
                        for j=1:N
                            labels{j} = sprintf('C_{%i}',isort(j));
                        end
                        t(ii) = Bars(C(:,obj.syn.inds{i})',obj.syn.chlabels,labels);
                        if ii>1, t(ii).chlabels = {}; end
                        if ~isempty(opt.syncolor) && length(opt.syncolor)==N
                            t(ii).opt.color = opt.syncolor;
                        end
                        
                    case 'spatiotemporal'
                        for j=1:N
                            chlabels{j} = sprintf('C%i',isort(j));
                        end
                        labels{1} = '';
                        t(ii) = Bars(C(:,obj.syn.inds{i}),chlabels,labels);
                        if ii>1, t(ii).chlabels = {}; end
                        
                        
                end
            end
            
        end
        %------------------------------------------------------------------
        
        function n = nch(obj)
            % number of channels in data matrix
            n = size(obj.data(1).data,1);
        end
        
        %------------------------------------------------------------------
        function n = nonnegch(obj)
            % number of non-negative channels in data matrix
            n = length(obj.data(1).nonnegch);
        end
    end
    
    
    methods (Static)
    end
    
end

%------------------------------------------------------------------
% subfunctions
%------------------------------------------------------------------

function par = getDefPar

par.type = 'EmgData';
par.order = 0;
par.delay = 0;
par.chlabels = [];

end

%------------------------------------------------------------------
function opt = getDefOpt
% get default options

opt.emgFilter.type = 'fir1';
opt.emgFilter.par = [50 .04]; % 20 Hz @ 1KHz sampling
opt.emgFilter.resample = 1;
opt.emgFilter.resample_period = .01; % resampling period [s]

opt.kinFilter.type = 'butter';
opt.kinFilter.par  = [2 10/(100/2)]; %2 order 10Hz @ 100Hz sampling rate

opt.normalize.type = 2;  % max absolute value of each channels
opt.normalize.isect = []; % sections to use for computing max

opt.average.gr = {};  % each cell contains the ids of the trials to average
opt.average.grcond = {}; % group selection criteria type (e.g. {'type',[1 1 1;1 2 1]})
opt.average.trange = [-.5 1];

opt.find.type = 'spatial';
opt.find.algo = 'nmf';
opt.find.N    = [1:8];
% opt.find.nrep = 5;
% opt.find.bestrsqrep  = 1; % 0-> save all reps; 1-> save only best rep
% opt.find.niter = [5 5 1e-4];
% opt.find.plot = 0;

opt.plot.type    = 'W';
opt.plot.subtype = 'barh';    % plot subtype
opt.plot.N       = [];    % number of synergies to plot (choose rep with max R^2)
opt.plot.iset    = 1;     % first set
opt.plot.irep    = 1;     % first repetition
opt.plot.isort   = [];    % synergy sort order
opt.plot.isect   = 1;
opt.plot.syncolor = {};
opt.plot.papertype = 'a4';

opt.verbose = 0;  % print additional messages

end