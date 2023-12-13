% EmgKinData: class to load, preprocess,and plot EMG and kinematics data
%
% obj = EmgKinData(varargin)
%
% list of methods
% -----------------------------------------------------------------
% EmgKinData            class constructor
% subtract              subtract tonic or baseline level to EMG data
% resample              average emgs over time interval dt and adjust kin time accordingly
% filter                filter data
% getDefFilterOpt       define the defaults of filter function
% selectCh              select channels
% merge                 merge kin and emg data in one data matrix
% shift                 introduce a delay between emg and kin data
% split                 divide data into kin and emg objects
% average               average data across trials
% getDefAverageOpt      define the defaults of average function
% normalize             normalize data across trials
% getDefNormalizeOpt    define the defaults of normalize function
% getDefSubtractOpt     define the defaults of subtract function
% getData               get data matrix
% getNsamp              get number of samples
% plot                  plot EMG and KIN data
% getDefPlotOpt         define the defaults of plot function
% duration              compute duration of data for each trials
% datalim               return mim and max value of data
% tsamp                 returns data mean sampling interval (rounded to us)
% timerange             get time range of data of each trial
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%

classdef EmgKinData
    
    properties
        data;
        emg;
        kin;
        time;
        trialId = [];
        chlabels = {};
        nonnegch = [];
        order;
        delay;
        normalized;
    end
    
    methods
        %------------------------------------------------------------------
        function obj = EmgKinData(varargin)
            
            if nargin<1
                return
            end
            
            if isstruct(varargin{1}) && isfield(varargin{1},'data') && isfield(varargin{1},'time')
                %data,trialId,chlabels,order
                emgkin = varargin{1};
                
                ntrial = length(emgkin);
                
                for i=1:ntrial
                    obj(i).data = emgkin(i).data;
                    obj(i).time = emgkin(i).time;
                    if nargin>1 && isequal(length(varargin{2}),ntrial)
                        trialId = varargin{2};
                        obj(i).trialId = trialId(i);
                    else
                        obj(i).trialId = i;
                    end
                    nch = size(obj(i).data,1);
                    if nargin>2 && iscell(varargin{3}) && isequal(length(varargin{3}),nch)
                        chlabels = varargin{3};
                        obj(i).chlabels = chlabels;
                    else
                        for j=1:nch
                            chlabels{j} = sprintf('emg%02i',j);
                        end
                        obj(i).chlabels = chlabels;
                    end
                    
                    if nargin>3
                        obj(i).order = varargin{4};
                    else
                        obj(i).order = 0;
                    end
                    if nargin>4
                        obj(i).nonnegch = varargin{5};
                    else
                        obj(i).nonnegch = [1:18];   %%%%%% TO BE CHECKED
                    end
                    
                end
            elseif isa(varargin{1},'EmgData') && isa(varargin{2},'KinData')
                %emg,kin,trialId,chlabels,order
                emg = varargin{1};
                kin = varargin{2};
                
                if length(emg) ~= length(kin)
                    warning('Emg and Kin data must have same number of trials')
                    return
                end
                ntrial = length(emg);
                
                if nargin>4 && not(isempty(varargin{5}))
                    delay = varargin{5};
                else
                    delay = 0;
                    
                end
                
                if nargin>3 && not(isempty(varargin{4}))
                    order = varargin{4};
                else
                    order = 0;
                    
                end
                
                for i=1:ntrial
                    obj(i).order = order;
                    obj(i).normalized = 0;
                    
                    if nargin>2 && isequal(length(varargin{3}),ntrial)
                        trialId = varargin{3};
                        obj(i).trialId = trialId(i);
                    else
                        obj(i).trialId = i;
                    end
                    
                    obj(i).emg = emg(i);
                    obj(i).kin = kin(i);
                    
                    obj(i) = shift(obj(i),delay);
                    obj(i) = merge(obj(i));
                    
                    nch = size(emg(i).data,1) + size(kin(i).data,1);
                    nemgch = size(emg(i).data,1);
                    nkinch = size(kin(i).data,1);
                    
                    for j=1:nemgch
                        chlabels{j} = emg(i).chlabels{j};
                    end
                    for j=(1+nch-nkinch):nch
                        jj = j - (nch-nkinch);
                        
                        if order>1
                            chlabels{j} = sprintf('d^%i/dt^%i %s',order,order,kin(i).chlabels{jj});
                        elseif order>0
                            chlabels{j} = sprintf('d/dt % s', kin(i).chlabels{jj});
                        else
                            chlabels{j} = kin(i).chlabels{jj};
                        end
                        
                    end
                    obj(i).chlabels = chlabels;
                    
                    obj(i).nonnegch = [1:size(emg(i).data,1)];
                end
                
            end
        end
        %------------------------------------------------------------------
        function obj = subtract(obj,opt)
            % subtract tonic or baseline level to EMG data
            
            % set options or use defaults
            defopt = obj.getDefSubtractOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            ntrial = length(obj);
            e = [obj.emg];
            e = subtract(e,opt);
            
            for i=1:ntrial
                obj(i).emg = e(i);
                
            end
            obj = merge(obj);
        end
        
        %------------------------------------------------------------------
        function obj = resample(obj,dt)
            % average emgs over time interval dt and adjust kin time
            % accordingly
            
            ntrial = length(obj);
            
            for i=1:ntrial
                
                obj(i).emg = resample(obj(i).emg,dt);
                obj(i).kin = resample(obj(i).kin,dt);
                
            end
            obj = merge(obj);
        end
        
        
        %------------------------------------------------------------------
        function obj = filter(obj,opt)
            % set options or use defaults
            defopt = obj.getDefFilterOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            ntrial = length(obj);
            
            for i=1:ntrial
                obj(i).emg = filter(obj(i).emg,opt.emgFilter);
                obj(i).kin = filter(obj(i).kin,opt.kinFilter);
                
            end
            obj = merge(obj);
        end
        
        %------------------------------------------------------------------
        function opt = getDefFilterOpt(obj)
            opt.emgFilter.type = 'fir1';
            opt.emgFilter.par = [50 .04];  % 20 Hz @ 1KHz EMG sampling rate
            opt.emgFilter.resample = 1;
            opt.emgFilter.resample_period = .01; % resampling period [s]
            
            opt.kinFilter.type = 'butter';
            opt.kinFilter.par  = [2 10/(100/2)]; %2 order 10Hz @ 100Hz sampling rate
            
        end
        
        %------------------------------------------------------------------
        function obj = selectCh(obj,chind)
            % select channels
            nkin = length(obj);
            for i=1:nkin
                obj(i).data = obj(i).data(chind,:);
            end
            obj(i).chlabels = obj(i).chlabels(chind);
            
        end
        
        %------------------------------------------------------------------
        function obj = merge(obj,order)
            
            if (nargin<2 || isempty(order)) && isempty(obj(1).order)
                
                order = 0;
            elseif not(isempty(obj(1).order))
                order = obj(1).order;
                
            end
            
            ntrial = length(obj);
            
            for i=1:ntrial
                obj(i).order = order;
                e = obj(i).emg;
                k = obj(i).kin;
                
                
                if ~isequal(tsamp(e),tsamp(k))
                    if obj(i).trialId ==1, disp('size(emg,2) must be equal to length(kin), resampling'), end
                    e = filter(e);
                    dt = 0.01; %10 ms
                    e = resample(e,dt);
                    k = resample(k,dt);
                end
                
                if ~isequal(duration(e),duration(k))
                    if obj(i).trialId ==1, disp('size(emg,2) must be equal to length(kin), slicing'), end
                    ktr = timerange(k);
                    etr = timerange(e);
                    
                    %find smallest common interval
                    if ktr(1)<0 && etr(1)<0
                        tr(1) = max([ktr(1) etr(1)]);
                    elseif ktr(1)<0 || etr(1)<0
                        tr(1) = max([ktr(1) etr(1)]);
                    elseif ktr(1)>0 && etr(1)>0
                        tr(1) = min([ktr(1) etr(1)]);
                    end
                    tr(2) = min([ktr(2) etr(2)]);
                    
                    timek = k.time;
                    indk = find(timek>=tr(1) & timek<=tr(end)); % find samples contained into chosen interval
                    timee = e.time;
                    inde = find(timee>=tr(1) & timee<=tr(end)); % find samples contained into chosen interval
                    
                    if not(isequal(length(indk),length(inde)))
                        
                        [mm,im] = min(abs(timee-timek(indk(1)))); % find sample closest to first data sample
                        inde = im + [1:length(indk)] - 1;
                        
                    end
                    
                    k.data = k.data(:,indk);
                    k.time = k.time(indk);
                    e.data = e.data(:,inde);
                    e.time = e.time(inde);
                end
                
                if not(obj(i).normalized)
                    if order == 0
                        obj(i).data = [e.data; k.data];
                    elseif order == 1
                        obj(i).data = [e.data; velocity(k,k.time).data];
                    elseif order == 2
                        obj(i).data = [e.data; acceleration(k,k.time).data];
                    end
                else
                    obj(i).data = [e.data; k.data];
                end
                obj(i).time = k.time;
            end
        end
        
        
        %------------------------------------------------------------------
        function obj = shift(obj,delay)
            
            if (nargin<2 || isempty(delay)) && isempty(obj(1).delay)
                
                delay = 0;
            elseif not(isempty(obj(1).delay))
                delay = obj(1).delay;
                
            end
            
            
            ntrial = length(obj);
            
            for i=1:ntrial
                obj(i).delay = delay;
                e = obj(i).emg;
                k = obj(i).kin;
                
                k = shift(k,delay);
                
                obj(i).kin = k;
                
            end
        end
        
        %------------------------------------------------------------------
        function obj = split(obj)
            
            ntrial = length(obj);
            
            for i=1:ntrial
                nch = size(obj(i).data,1); % number of channels
                emgch = obj(i).nonnegch;
                nemgch = length(emgch);
                kinch = [(nemgch+1):nch];
                
                obj(i).emg.data = obj(i).data(emgch,:);
                obj(i).kin.data = obj(i).data(kinch,:) ;
                obj(i).emg.time = obj(i).time;
                obj(i).kin.time = obj(i).time ;
                obj(i).emg.chlabels = obj(i).chlabels(emgch);
                obj(i).kin.chlabels = obj(i).chlabels(kinch);
                
            end
        end
        
        
        
        %------------------------------------------------------------------
        function objav = average(obj,opt)
            % average EMG/KIN data across trials
            
            % set options or use defaults
            defopt = obj.getDefAverageOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            e = [obj.emg];
            k = [obj.kin];
            
            eav = average(e,opt);
            kav = average(k,opt);
            
            
            objav = EmgKinData(eav,kav,[],obj(1).order,0); %set delay to zero as data are already delayed
            
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefAverageOpt(obj)
            
            nkin = length(obj);
            opt.gr = {[1:nkin]};
            
            opt.tref = zeros(1,nkin);
            
            tr = timerange(obj);
            trc = [max(tr(:,1)) min(tr(:,2))];
            if diff(trc)>0
                opt.trange = trc;
            else
                opt.trange = [];
            end
            
        end
        
        %------------------------------------------------------------------
        function [obj,objnorm] = normalize(obj,opt)
            %normalize data in emg/kin amplitude
            %
            %   type  action
            %   -------------------------------------------------------------
            %   0         use normdata [nch,1] for normalization of each channel
            %   1         normalize to max of any channel
            %   2         normalize each channel to max in that channel
            %   32        max of each abs(channel) and same range for emg [0 2] and trq [-1 1]
            
            
            % set options or use defaults
            defopt = obj.getDefNormalizeOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            ntrial = length(obj);
            
            % check for order
            for i=1:ntrial
                obj(i).normalized = 1;
                if obj(i).order == 1
                    obj(i).kin.data = velocity(obj(i).kin,obj(i).kin.time).data;
                elseif obj(i).order == 2
                    obj(i).kin.data = acceleration(obj(i).kin,obj(i).kin.time).data;
                end
            end
            e = [obj.emg];
            k = [obj.kin];
            
            e = normalize(e,opt);
            k = normalize(k,opt);
            
            for i=1:ntrial
                obj(i).emg = e(i);
                obj(i).kin = k(i);
            end
            obj = merge(obj);
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefNormalizeOpt(obj)
            
            nkin = length(obj);
            
            opt.type = 32;  % max absolute value of each channels
            opt.isect = [1:nkin]; % sections to use for computing max
            opt.normdata = [];
            opt.nonnegch = [];
        end
        
        %------------------------------------------------------------------
        function opt = getDefSubtractOpt(obj)
            
            opt.type = 'tonic';    % subtract tonic activity to get phasic EMG data
            opt.t_pre   = [-.4 -.2];    % interval before onset for initial level
            opt.t_post  = [.2 .4];      % interval after end for final level
            opt.t_onset = [];
            opt.t_end = [];
            opt.clip    = 1;         % clip to zero after subtraction
            
        end
        %------------------------------------------------------------------
        function [data,datapar] = getData(obj,type,isect)
            % get data matrix
            
            nkin = length(obj);
            if nargin<2, type = 'spatial'; end
            if nargin<3, isect = [1:nkin]; end
            isect = intersect([1:nkin],isect);
            nsect = length(isect);
            
            nch = size(obj(isect(1)).data,1);
            
            switch type
                
                case 'spatial'
                    % rows are channels, columns are time samples x trials
                    data = zeros(nch,obj.getNsamp(isect));
                    inds = cell(1,nsect);
                    isamp = 0;
                    for ii=1:nsect
                        i = isect(ii);
                        nsamp = length(obj(i).time);
                        data(:,isamp+[1:nsamp]) = obj(i).data;
                        inds{ii} = isamp+[1:nsamp];
                        isamp = isamp + nsamp;
                    end
                    
                    %
                    % other options to be implemented
                    %
                    %                 case 'temporal'
                    %                     % rows are time samples, columns are channels x trials
                    %                     ntime = length(k(isect(1)).time);
                    %                     data = zeros(ntime,nch*nsect);
                    %                     inds = cell(1,nsect);
                    %                     if ~isequalinterval(k(isect))
                    %                         warning('data must have the same interval for temporal synergies')
                    %                         return
                    %                     end
                    %
                    %                     nsamp = nch;
                    %                     isamp = 0;
                    %                     for ii=1:nsect
                    %                         i = isect(ii);
                    %                         data(:,isamp+[1:nsamp]) = k(i).data';
                    %                         inds{ii} = isamp+[1:nsamp];
                    %                         isamp = isamp + nsamp;
                    %                     end
                    %
                    %                 case 'spatiotemporal'
                    %                     % rows are channels x time samples, columns trials
                    %                     ntime = length(k(isect(1)).time);
                    %                     data = zeros(ntime*nch,nsect);
                    %                     inds = cell(1,nsect);
                    %                     if ~isequalinterval(k(isect))
                    %                         warning('data must have the same interval for spatiotemporal synergies')
                    %                         return
                    %                     end
                    %
                    %                     nsamp = 1;
                    %                     isamp = 0;
                    %                     for ii=1:nsect
                    %                         i = isect(ii);
                    %                         data(:,isamp+[1:nsamp]) = k(i).data(:);
                    %                         inds{ii} = isamp+[1:nsamp];
                    %                         isamp = isamp + nsamp;
                    %                     end
                    
            end
            
            datapar.nch = nch;
            datapar.inds = inds;
            datapar.nonnegch = obj(1).nonnegch;
            datapar.chlabels = obj(1).chlabels;
        end
        
        %------------------------------------------------------------------
        function nsamp = getNsamp(obj,isect)
            % get number of samples
            
            nkin = length(obj);
            if nargin<2, isect = [1:nkin]; end
            isect = intersect([1:nkin],isect);
            
            nsamp = 0;
            for i=isect
                nsamp = nsamp + length(obj(i).time);
            end
            
        end
        
        %------------------------------------------------------------------
        function hha = plot(obj,opt)
            % plot EMG/KIN data
            
            % set options or use defaults
            defopt = obj.getDefPlotOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            opte = opt;
            optk = opt;
            
            obj = split(obj);
            
            e = [obj.emg];
            k = [obj.kin];
            
            %
            % figure and axes
            %
            
            if isempty(opt.axes)
                opte.pos = opt.pos(1,:);
                optk.pos = opt.pos(2,:);
            else
                opte.axes = opt.axes(1,:);
                optk.axes = opt.axes(2,:);
            end
            
            if isempty(opt.figure)
                hf = figure;
                opte.figure = hf;
                optk.figure = hf;
                
            end
            
            % color
            if not(opt.rec)
                opte.color = zeros(1,3);
                optk.color = .8*ones(1,3);
                
            end
            hae = plot(e,opte);
            optk.usetitle = 0;
            hak = plot(k,optk);
            
            if nargout>0, hha = [hae;hak]; end
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefPlotOpt(obj)
            
            nkin = length(obj);
            nnonnegch = length(obj(1).nonnegch);
            nch = size(obj(1).data,1); % number of channels
            
            % chsel
            opt.emgsel = obj.nonnegch;
            opt.kinsel = [(nnonnegch+1):nch]-nnonnegch;
            % isect
            opt.isect = [1:nkin];
            
            % tref
            opt.tref = [];
            
            % events
            opt.event_code = {};
            opt.event_time = {};
            opt.events_color = {};
            opt.events_style = {};
            
            % title
            opt.usetitle = 1;
            opt.chtitle = {};
            
            % plotting options
            opt.figure = [];
            opt.axes   = [];
            opt.pos = [.08 .40 .9 .55; .08 .05 .9 .30];
            
            opt.overlap = 0;
            opt.spacing = .01;
            
            opt.xlim = [];
            opt.ylim = []; % for each individual trace
            opt.chscale = 0;
            opt.chscalelabel = '';
            opt.fill = 0;
            opt.color = 'k';
            
            opt.rec = 0;
            
        end
        
        %------------------------------------------------------------------
        function du = duration(obj)
            % compute duration of data for each trials
            nkin = length(obj);
            for i=1:nkin
                du(i) = range(obj(i).time);
            end
            
        end
        
        
        %------------------------------------------------------------------
        function val = datalim(obj)
            %return mim and max value of data
            
            nkin = length(obj);
            for i=1:nkin
                valmin(i) = min(obj(i).data(:));
                valmax(i) = max(obj(i).data(:));
            end
            val = [min(valmin) max(valmax)];
        end
        
        %------------------------------------------------------------------
        function t = tsamp(obj,prec)
            %returns data mean sampling interval (rounded to us)
            if nargin<2, prec=10^-6; end
            
            nkin = length(obj);
            for i=1:nkin
                t(i) = mean(round(diff(obj(i).time)/prec))*prec;
            end
            
        end
        
        %------------------------------------------------------------------
        function tr = timerange(obj)
            %get time range of data of each trial
            
            nkin = length(obj);
            for i=1:nkin
                tr(i,:) = obj(i).time([1 end]);
            end
            
        end
        
    end
    
    
    %------------------------------------------------------------------
    methods (Static)
        
    end
    
end

%------------------------------------------------------------------
% subfunctions
%------------------------------------------------------------------
