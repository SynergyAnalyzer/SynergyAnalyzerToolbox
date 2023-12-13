% EmgData: class to load, preprocess,and plot EMG data
%
% obj = EmgData(data,trialId,chlabels)
%
% list of methods
% -----------------------------------------------------------------
% EmgData               class constructor
% subtract              subtract tonic or baseline level to EMG data
% resample              average emgs over time interval dt
% filter                filter EMG data
% getDefFilterOpt       define the defaults of filter function
% selectCh              select channels
% average               average data across trials
% getDefAverageOpt      define the defaults of average function
% normalize             normalize data across trials
% getDefNormalizeOpt    define the defaults of normalize function
% getDefSubtractOpt     define the defaults of subtract function
% getData               get data matrix
% getNsamp              get number of samples
% plot                  plot EMG data
% getDefPlotOpt         define the defaults of plot function
% duration              compute duration of data for each trials
% emglim                return mim and max value of EMG data
% tsamp                 returns data mean sampling interval (rounded to us)
% timerange             get time range of data of each trial
% max                   compute max of each EMG channel
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%


classdef EmgData
    
    properties
        data;
        time;
        trialId = [];
        chlabels = {};
    end
    
    methods
        %------------------------------------------------------------------
        function obj = EmgData(data,trialId,chlabels)
            if nargin<1
                return
            end
            if ~isstruct(data) || ~isfield(data,'emg') || ~isfield(data,'emgtime')
                warning('Data input must be a structure with .emg and .emgtime fields')
                return
            end
            ntrial = length(data);
            for i=1:ntrial
                if ~isequal(size(data(i).emg,2),length(data(i).emgtime))
                    warning('size(data(%i).emg,2) must be equal to length(data(%i).emgtime)',i,i)
                    data(i).emgtime = [1:size(data(i).emg,2)];
                end
                obj(i).data = data(i).emg;
                obj(i).time = data(i).emgtime;
                if nargin>1 && isequal(length(trialId),ntrial)
                    obj(i).trialId = trialId(i);
                else
                    obj(i).trialId = i;
                end
                nch = size(data(i).emg,1);
                if nargin>2 && iscell(chlabels) && isequal(length(chlabels),nch)
                    obj(i).chlabels = chlabels;
                else
                    for j=1:nch
                        emglabels{j} = sprintf('emg%02i',j);
                    end
                    obj(i).chlabels = emglabels;
                end
            end
            
        end
        
        %------------------------------------------------------------------
        function obj = filter(obj,opt)
            % filter EMG data
            %
            %   type       par       notes
            %   =============================================================
            %   'fir1'     [N Wn]    low pass finite impulse response filter
            %   'fir1'     [N W1 W2] band pass finite impulse response filter
            %   'butter'   [N Wn]    low pass Nth order Butterworth filter
            %   'butter'   [N W1 W2] band pass 2Nth order Butterworth filter
            %   'rectify'  []        rectification
            %   'rectify'  [n]       rectification and resampling
            %   'submean'  []        rectification after mean subtraction
            %   'resample' [N WN n]  resample
            %   'rms'      [N n]     root mean square and resample
            %   'average'  [N n]     moving average and resample
            %   'high'     [N W2]    high pass FIR1 on NON-rectified EMGs
            %   'notch'    [W0 Q]    IIR notch filter with notch freqency (W0*Fs/2)
            %                        quality factor Q (Q = W0/bw)
            %   'cliptozero'         clip negative values to zero
            
            
            % set options or use defaults
            defopt = obj.getDefFilterOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            nemg = length(obj);
            
            for i=1:nemg
                
                [nch,nsamp] = size(obj(i).data);
                
                switch opt.type
                    
                    case {'fir1','butter'}
                        if length(opt.par)<2
                            warning('fir1 type requires two parms (N,Wn)')
                            return
                        end
                        N = opt.par(1);
                        if length(opt.par)==2
                            Wn = opt.par(2);
                        else
                            Wn = opt.par(2:3);
                        end
                        if N>(3*nsamp)
                            warning('filter order too large for given data!')
                            return
                        end
                        switch opt.type
                            case 'fir1'
                                B = fir1(N,Wn);
                                A = 1;
                            case 'butter'
                                % note N is order!!
                                [B,A] = butter(N,Wn);
                        end
                        obj(i).data = filtfilt(B,A,abs(obj(i).data'))';
                        
                    case 'high'
                        if length(opt.par)<2
                            warning('high type requires two parms (N,Wn)')
                            return
                        end
                        N = opt.par(1);
                        if length(opt.par)==2
                            Wn = opt.par(2);
                        else
                            Wn = opt.par(2:3);
                        end
                        if N>(3*nsamp)
                            warning('filter order too large for given data!')
                            return
                        end
                        B = fir1(N,Wn,'high');
                        obj(i).data = filtfilt(B,1,obj(i).data')';
                        
                    case 'stop'
                        if length(opt.par)<2
                            warning('high type requires two parms (N,Wn)')
                            return
                        end
                        N = opt.par(1);
                        Wn = opt.par(2:3);
                        if N>(3*nsamp)
                            warning('filter order too large for given data!')
                            return
                        end
                        B = fir1(N,Wn,'stop');
                        obj(i).data = filtfilt(B,1,obj(i).data')';
                        
                    case 'notch'
                        if length(opt.par)<2
                            warning('notch type requires two parms (W0,Q)')
                            return
                        end
                        if ~exist('iirnotch','file')
                            warning('notch type requires Filter Design Toolbox')
                            return
                        end
                        nfreq = size(opt.par,1);
                        for j=1:nfreq
                            W0 = opt.par(j,1);
                            bw = W0/opt.par(j,2);
                            [B,A] = iirnotch(W0,bw);
                            obj(i).data = filtfilt(B,A,obj(i).data')';
                        end
                        
                    case 'resample'
                        if length(opt.par)<3
                            warning('resample type requires 3 parms (N,Wn,n)')
                            return
                        end
                        N  = opt.par(1);
                        Wn = opt.par(2);
                        n  = opt.par(3);
                        if N>(3*nsamp)
                            warning('filter order too large for given data!')
                            return
                        end
                        if rem(N,2)==0, N=N+1; end  % resample requires odd filter
                        B = fir1(N,Wn);
                        obj(i).data = resample(abs(obj(i).data'),1,n,B)';
                        ndata = size(obj(i).data,2);
                        time = zeros(1,ndata);
                        ntime = length(obj(i).time);
                        for j=1:ndata
                            time(j) = mean(obj(i).time(1+(j-1)*n:min(2*j,ntime)));
                        end
                        obj(i).time = time;
                        
                    case 'rectify'
                        if ~isempty(opt.par)
                            n = opt.par(1);
                        else
                            n = 1;
                        end
                        datatemp = abs(obj(i).data);
                        ind = [1:n:nsamp];
                        obj(i).data = datatemp(:,ind);
                        obj(i).time = obj(i).time(ind);
                        
                    case 'submean'
                        if ~isempty(opt.par)
                            n = opt.par(1);
                        else
                            n = 1;
                        end
                        datatemp = obj(i).data;
                        datatemp = data - mean(datatemp,2)*ones(1,nsamp);
                        ind = [1:n:nsamp];
                        obj(i).data = datatemp(:,ind);
                        obj(i).time = obj(i).time(ind);
                        
                    case {'rms','average'}
                        N  = opt.par(1);
                        n  = opt.par(2);
                        nresamp = floor(nsamp/n);
                        
                        datatemp = zeros(nch,nresamp);
                        for j=1:nresamp
                            if ~rem(N,2) % even
                                N = N+1;
                            end
                            ind = [(1-N)/2+1:(N+1)/2];
                            if strcmp(opt.type,'rms')
                                datatemp(:,j) = sqrt(mean(obj(i).data(:,min(nsamp,max(1,n*(j-1)+ind))).^2,2));
                            elseif strcmp(opt.type,'average')
                                datatemp(:,j) = mean(abs(obj(i).data(:,min(nsamp,max(1,n*(j-1)+ind)))),2);
                            end
                        end
                        obj(i).data = datatemp;
                        obj(i).time = obj(i).time(1:n:nresamp*n);
                        
                    case 'cliptozero'
                        obj(i).data = obj(i).data .* (obj(i).data>0);
                        
                        
                    otherwise
                        warning('Unknown filter type! EMGs were not filtered.')
                        
                end
                
            end
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefFilterOpt(obj)
            
            opt.type = 'rectify';
            opt.par  = [];
            
        end
        
        
        %------------------------------------------------------------------
        function obj = selectCh(obj,chind)
            % select channels
            nemg = length(obj);
            for i=1:nemg
                obj(i).data = obj(i).data(chind,:);
            end
            obj(i).chlabels = obj(i).chlabels(chind);
            
        end
        
        %------------------------------------------------------------------
        function obj = resample(obj,dt)
            % average emgs over time interval dt
            %
            %   obj = resample(obj,dt)  computes mean of adjacent intervals of duration dt [s]
            %
            %      OR
            %
            %   m = resample(obj,t_range) compute mean values between t_range(1) and
            %   t_range(2); m -> [nch,nemg]
            
            if nargin<2
                dt = timerange(obj); % 0-> dt, 1->t_range
            end
            
            nemg = length(obj);
            
            if length(dt)==1
                type = 0;
            elseif size(dt,2)==2
                type = 1;
                if isequal(size(dt,1),nemg)
                    t_range = dt;
                else
                    t_range = ones(nemg,1)*dt(1,:);
                end
            else
                warning('invalid input!')
                return
            end
            
            switch type
                
                case 0 % dt
                    for i=1:nemg
                        [nch,nsamptot] = size(obj(i).data);
                        t_sample = tsamp(obj(i)); % sampling interval
                        nsamp = round(dt/t_sample); % number of original samples to sum for each integrated sample
                        nintervals = floor(nsamptot/nsamp); % number of integrated intervals
                        datatemp = zeros(nch,nintervals);       % allocate space
                        for k=1:nintervals
                            ind = nsamp*(k-1)+[1:nsamp];
                            datatemp(:,k) = mean(obj(i).data(:,ind),2);
                        end
                        obj(i).data = datatemp;
                        obj(i).time = obj(i).time([1:nsamp:nsamp*nintervals])+t_sample*(nsamp-1)/2;  % if nsamp==1 -> obj.time does not change
                        
                    end
                    
                case 1 % mean over t_range
                    nch = size(obj(1).data,1);
                    m = zeros(nch,nemg);
                    for i=1:nemg
                        ind = find(obj(i).time>=t_range(i,1) & obj(i).time<=t_range(i,2));
                        m(:,i) = mean(obj(i).data(:,ind),2);
                    end
                    obj = m;
                    
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
            
            % check options
            nemg = length(obj);
            switch opt.type
                case 'tonic'
                    if ~isequal(length(opt.t_onset),nemg)
                        warning('onset times must be provided for tonic subtraction')
                        return
                    end
                    
            end
            
            % loop on trials
            for i = 1:nemg
                
                switch opt.type
                    
                    case 'tonic'
                        t_on = opt.t_onset(i);
                        t_off = opt.t_end(i);
                        
                        tr = timerange(obj(i));
                        on_range = t_on+opt.t_pre;
                        off_range = t_off+opt.t_post;
                        
                        if (on_range(2) <tr(1) || off_range(1) > tr(2))
                            warning('wrong tails selection for computing tonic subtraction, proceeding without subtracting tonic component')
                            
                        else
                            if (on_range(1) <tr(1) )
                                on_range(1) = tr(1);
                                warning('computing tonic component on pre-portion of data available, not entire 200ms selected')
                            end
                            
                            if ( off_range(2) > tr(2))
                                off_range(2) = tr(2);
                                warning('computing tonic component on post-portion of data available, not entire 200ms selected')
                            end
                            
                            val_on = resample(obj(i),on_range);
                            val_off = resample(obj(i),off_range);
                            
                            datatonic = zeros(size(obj(i).data));
                            time = obj(i).time;
                            ntime = length(time);
                            
                            [mm,ind_on]  = min(abs(time-t_on));
                            [mm,ind_off] = min(abs(time-t_off));
                            
                            datatonic(:,1:ind_on-1) = val_on * ones(1,ind_on-1);
                            datatonic(:,ind_off+1:ntime) = val_off * ones(1,ntime-ind_off);
                            indramp  = [ind_on:ind_off];
                            datatonic(:,indramp) = interp1([ind_on ind_off],[val_on,val_off]',indramp)';
                            
                            obj(i).data = obj(i).data - datatonic;  % phasic
                        end
                        if opt.clip
                            obj(i).data = obj(i).data.*(obj(i).data>0); % clip negative values to zero
                        end
                        
                    case 'mean'
                        obj(i).data = obj(i).data - mean(obj(i).data,2)*ones(1,size(obj(i).data,2));
                        
                end
                
            end
            
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
        function obj = average(obj,opt)
            % average EMG data across trials
            
            % set options or use defaults
            defopt = obj.getDefAverageOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            % loop on groups
            gr = opt.gr;
            tref = opt.tref;
            trange = opt.trange;
            
            ngr = length(gr);
            
            ts    = tsamp(obj(gr{1}(1)));
            nch   = size(obj(1).data,1);      % number of channels
            
            tav = [trange(1):ts:trange(2)];  % times of averaged emgs
            nsampav = length(tav);  % number of samples in average
            
            for i=1:ngr
                
                datatemp = zeros(nch,nsampav);
                ndata    = zeros(1,nsampav);
                
                chlabels = obj(gr{i}(1)).chlabels;
                eavtemp.emgtime = tav;
                
                for j=1:length(gr{i})
                    jj = gr{i}(j);
                    time = obj(jj).time - tref(jj);
                    ind = find(time>=tav(1) & time<=tav(end)); % find samples contained into chosen interval
                    if ~isempty(ind)
                        if isequal(nsampav,length(ind))
                            indav = [1:nsampav];
                        else
                            [mm,im] = min(abs(tav-time(ind(1)))); % find sample closest to first data sample
                            indav = im + [1:length(ind)] - 1;
                        end
                        data_j = obj(jj).data(:,ind);
                        datatemp(:,indav) = datatemp(:,indav) + data_j;
                        ndata(indav) = ndata(indav) + 1;
                    end
                end
                
                % mean
                datam = datatemp ./ (ones(nch,1)*(ndata + (ndata==0)));
                eavtemp.emg = datam;
                
                eav(i) = EmgData(eavtemp,i,chlabels);
                
            end
            
            obj = eav;
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefAverageOpt(obj)
            
            nemg = length(obj);
            opt.gr = {[1:nemg]};
            
            opt.tref = zeros(1,nemg);
            
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
            %normalize data in emg amplitude
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
            
            nemg = length(obj);
            
            opt.isect = intersect([1:nemg],opt.isect);
            if isempty(opt.isect)
                warning('empty isect, using all sections')
                opt.isect = [1:nemg];
            end
            nch = size(obj(opt.isect(1)).data,1);
            switch opt.type
                case 1
                    normdata = max(max(obj(opt.isect)));
                case 2
                    normdata = max(obj(opt.isect),2);
                case 32
                    
                    normdata = .5*max(obj(opt.isect),2);
                    
            end
            
            for i=1:nemg
                
                [nch,nsamp] = size(obj(i).data);
                
                switch opt.type
                    
                    case 0
                        if  ~isequal(size(opt.normdata,1),nch)
                            warning('normdata missing or not valid!')
                            return
                        end
                        obj(i).data = obj(i).data ./ (opt.normdata*ones(1,nsamp));
                        objnorm(i).data = (opt.normdata*ones(1,nsamp));
                        
                    case 1
                        obj(i).data = obj(i).data / normdata;
                        objnorm(i).data = normdata .* ones(size(obj(i).data));
                        
                    case {2,32}
                        obj(i).data = obj(i).data ./ (normdata*ones(1,nsamp));
                        objnorm(i).data = normdata*ones(1,nsamp);
                        
                        
                        
                end
                
            end
        end
        
        %------------------------------------------------------------------
        function opt = getDefNormalizeOpt(obj)
            
            nemg = length(obj);
            
            opt.type = 2;  % max absolute value of each channels
            opt.isect = [1:nemg]; % sections to use for computing max
            opt.normdata = [];
            
        end
        
        %------------------------------------------------------------------
        function [data,datapar] = getData(obj,type,isect)
            % get data matrix
            
            nemg = length(obj);
            if nargin<2, type = 'spatial'; end
            if nargin<3, isect = [1:nemg]; end
            isect = intersect([1:nemg],isect);
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
                    
                case 'temporal'
                    % rows are time samples, columns are channels x trials
                    ntime = length(obj(isect(1)).time);
                    data = zeros(ntime,nch*nsect);
                    inds = cell(1,nsect);
                    if ~isequalinterval(obj(isect))
                        warning('data must have the same interval for temporal synergies')
                        return
                    end
                    
                    nsamp = nch;
                    isamp = 0;
                    for ii=1:nsect
                        i = isect(ii);
                        data(:,isamp+[1:nsamp]) = obj(i).data';
                        inds{ii} = isamp+[1:nsamp];
                        isamp = isamp + nsamp;
                    end
                    
                case 'spatiotemporal'
                    % rows are channels x time samples, columns trials
                    ntime = length(obj(isect(1)).time);
                    data = zeros(ntime*nch,nsect);
                    inds = cell(1,nsect);
                    if ~isequalinterval(obj(isect))
                        warning('data must have the same interval for spatiotemporal synergies')
                        return
                    end
                    
                    nsamp = 1;
                    isamp = 0;
                    for ii=1:nsect
                        i = isect(ii);
                        data(:,isamp+[1:nsamp]) = obj(i).data(:);
                        inds{ii} = isamp+[1:nsamp];
                        isamp = isamp + nsamp;
                    end
                    
            end
            datapar.nch = nch;
            datapar.inds = inds;
            datapar.chlabels = obj(1).chlabels;
            
        end
        
        %------------------------------------------------------------------
        function nsamp = getNsamp(obj,isect)
            % get number of samples
            
            nemg = length(obj);
            if nargin<2, isect = [1:nemg]; end
            isect = intersect([1:nemg],isect);
            
            nsamp = 0;
            for i=isect
                nsamp = nsamp + length(obj(i).time);
            end
            
        end
        
        %------------------------------------------------------------------
        function hha = plot(obj,opt)
            % plot EMG data
            
            % set options or use defaults
            defopt = obj.getDefPlotOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            %
            % figure and axes
            %
            nsect = length(opt.isect);
            if ~isfield(opt,'axes') | any(~ishandle(opt.axes)) | any(~strcmp(get(opt.axes,'type'),'axes'))
                if isempty(opt.figure), hf = figure; else hf = figure(opt.figure); end
                if opt.overlap
                    ha = axes('Position',opt.pos);
                else
                    dur = duration(obj(opt.isect));
                    for i=1:nsect
                        pos_i = SizBox(opt.pos,1,dur,1,i,[opt.spacing 0]);
                        ha(i) = axes('Position',pos_i);
                    end
                end
            else
                hf = get(opt.axes(1),'parent');
                ha = opt.axes;
            end
            
            if length(opt.ylim)<2
                yl = emglim(obj(opt.isect));
            else
                yl = opt.ylim;
            end
            
            %
            % loop on trials
            %
            for ii=1:nsect
                i = opt.isect(ii);
                
                tt.data = obj(i).data(opt.emgsel,:);
                tt.time = obj(i).time;
                if ~isempty(opt.tref) && length(opt.tref)==nemg, tt.time = tt.time - opt.tref(i); end
                t(ii) = Traces(tt);
                if ii==1, t(ii).chlabels = obj(i).chlabels(opt.emgsel); end
                if opt.usetitle
                    if isempty(opt.emgtitle) || isempty(opt.emgtitle{i})
                        t(ii).label = sprintf('%i',obj(i).trialId);
                    else
                        t(ii).label = sprintf('%s (%i)',opt.emgtitle{i},obj(i).trialId);
                    end
                end
                if ~isempty(opt.xlim)
                    if isequal(size(opt.xlim,1),nemg)
                        xl = opt.xlim(i,:);
                    else
                        xl = opt.xlim(1,:);
                    end
                    t(ii).opt.xlim = xl;
                    t(ii).opt.autotrange = 0;
                end
                t(ii).opt.autoscale = 0;
                t(ii).opt.ylim = yl;
                t(ii).opt.fill = opt.fill;
                if (opt.overlap & ii==1) | (~opt.overlap & ii==nsect)
                    t(ii).opt.yscale = opt.emgscale;
                    t(ii).opt.yscalelabel = opt.emgscalelabel;
                end
                if ~isempty(opt.color), t(ii).prop.color = opt.color; end
                if ~isempty(opt.linewidth), t(ii).prop.linewidth = opt.linewidth; end
                
            end
            
            plot(t,ha);
            
            %
            % events
            %
            if ~isempty(opt.event_code)
                evtcolor = opt.events_color;
                evtstyle = opt.events_style;
                for ii=1:nsect
                    i = opt.isect(ii);
                    evtcode = opt.event_code{i};
                    evttime = opt.event_time{i};
                    nevt = length(evtcode);
                    if opt.overlap, axes(ha(1)), else, axes(ha(ii)); end
                    yl = ylim;
                    for j=1:nevt
                        line([1;1]*evttime(j),yl'*[1 1],'color',evtcolor{evtcode(j)},'linestyle',evtstyle{evtcode(j)})
                    end
                end
            end
            
            if nargout>0, hha = ha; end
            
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefPlotOpt(obj)
            
            nemg = length(obj);
            nch = size(obj(1).data,1); % number of channels
            
            % emgsel
            opt.emgsel = [1:nch];
            
            % isect
            opt.isect = [1:nemg];
            
            % tref
            opt.tref = [];
            
            % events
            opt.event_code = {};
            opt.event_time = {};
            opt.events_color = {};
            opt.events_style = {};
            
            % title
            opt.usetitle = 1;
            opt.emgtitle = {};
            
            % plotting options
            opt.figure = [];
            opt.axes   = [];
            opt.pos = [.09 .08 .88 .86];
            opt.overlap = 0;
            opt.spacing = .01;
            
            opt.xlim = [];
            opt.ylim = []; % for each individual trace
            opt.emgscale = 0;
            opt.emgscalelabel = '';
            opt.fill = 0;
            opt.color = 'k';
            opt.linewidth = [];
            
        end
        
        %------------------------------------------------------------------
        function du = duration(obj)
            % compute duration of EMG data for each trials
            nemg = length(obj);
            for i=1:nemg
                du(i) = range(obj(i).time);
            end
            
        end
        
        %------------------------------------------------------------------
        function val = isequalinterval(obj)
            % check if trials have the same time samples
            nemg = length(obj);
            val = true;
            timeref = obj(1).time;
            for i=2:nemg
                if ~isequal(obj(i).time,timeref)
                    val = false;
                    return
                end
            end
            
        end
        
        %------------------------------------------------------------------
        function val = emglim(obj)
            % return mim and max value of emg data
            
            nemg = length(obj);
            for i=1:nemg
                valmin(i) = min(obj(i).data(:));
                valmax(i) = max(obj(i).data(:));
            end
            val = [min(valmin) max(valmax)];
        end
        
        %------------------------------------------------------------------
        function t = tsamp(obj,prec)
            %returns data mean sampling interval (rounded to us)
            if nargin<2, prec=10^-6; end
            
            nemg = length(obj);
            for i=1:nemg
                t(i) = mean(round(diff(obj(i).time)/prec))*prec;
            end
            
        end
        
        %------------------------------------------------------------------
        function tr = timerange(obj)
            %get time range of EMG data of each trial
            
            nemg = length(obj);
            for i=1:nemg
                tr(i,:) = obj(i).time([1 end]);
            end
            
        end
        
        
        %------------------------------------------------------------------
        function val = max(obj,type)
            % compute max of each channel
            %
            % type = 0 => max of all trials
            % type = 1 => max of each trial
            % type = 2 => max of all trials abs
            
            if nargin<2, type = 0; end
            
            nemg = length(obj);
            
            for i=1:nemg
                [nch,nsamptot] = size(obj(i).data);
                if i==1, val = NaN * ones(nch,1); end
                switch type
                    case 0
                        val = max(val,max(obj(i).data,[],2));
                    case 1
                        val(:,i) = max(obj(i).data,[],2);
                    case 2
                        val = max(val,max(abs(obj(i).data),[],2));
                end
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
