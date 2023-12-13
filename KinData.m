% KinData: class to load, preprocess,and plot kinematics data
%
% obj = KinData(data,trialId,chlabels)
%
% list of methods
% -----------------------------------------------------------------
% KinData               class constructor
% subtract              subtract tonic or baseline level to kin data
% resample              resample kin object at timenew times
% filter                filter kin data
% getDefFilterOpt       define the defaults of filter function
% selectCh              select channels
% average               average data across trials
% getDefAverageOpt      define the defaults of average function
% normalize             normalize data across trials
% getDefNormalizeOpt    define the defaults of normalize function
% getDefSubtractOpt     define the defaults of subtract function
% getData               get data matrix
% getNsamp              get number of samples
% plot                  plot kin data
% getDefPlotOpt         define the defaults of plot function
% duration              compute duration of data for each trials
% isequalinterval       check if trials have the same time samples
% tsamp                 returns data mean sampling interval (rounded to us)
% kinlim                return mim and max value of kin data
% timerange             get time range of data of each trial
% max                   compute max of each kin channel
% velocity              return velocity at time t
% acceleration          return velocity at time t
% shift                 shift time base by tshift (t <- t + tshift)
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%


classdef KinData
    
    properties
        data;
        time;
        trialId = [];
        chlabels = {};
    end
    
    methods
        %------------------------------------------------------------------
        function obj = KinData(data,trialId,chlabels)
            if nargin<1
                return
            end
            if ~isstruct(data) || ~isfield(data,'pos') || ~isfield(data,'postime')
                warning('Data input must be a structure with .pos and .postime fields')
                return
            end
            ntrial = length(data);
            for i=1:ntrial
                if ~isequal(size(data(i).pos,2),length(data(i).postime))
                    warning('size(data(%i).pos,2) must be equal to length(data(%i).postime)',i,i)
                    data(i).postime = [1:size(data(i).pos,2)];
                end
                obj(i).data = data(i).pos;
                obj(i).time = data(i).postime;
                if nargin>1 && isequal(length(trialId),ntrial)
                    obj(i).trialId = trialId(i);
                else
                    obj(i).trialId = i;
                end
                nch = size(data(i).pos,1);
                if nargin>2 && iscell(chlabels) && isequal(length(chlabels),nch)
                    obj(i).chlabels = chlabels;
                else
                    for j=1:nch
                        chlabels{j} = sprintf('dof%02i',j);
                    end
                    obj(i).chlabels = chlabels;
                end
                
            end
            
        end
        
        %------------------------------------------------------------------
        function obj = filter(obj,opt)
            % filter kin data
            %
            %   type       par       notes
            %   =============================================================
            %   'fir1'     [N Wn]    low pass finite impulse response filter
            %   'fir1'     [N W1 W2] band pass finite impulse response filter
            %   'butter'   [N Wn]    low pass Nth order Butterworth filter
            %   'butter'   [N W1 W2] band pass 2Nth order Butterworth filter
            %   'submean'  []        rectification after mean subtraction
            %   'resample' [N WN n]  resample
            %   'average'  [N n]     moving average and resample
            %   'high'     [N W2]    high pass FIR1 on NON-rectified kins
            %   'notch'    [W0 Q]    IIR notch filer with notch freqency (W0*Fs/2)
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
            
            nkin = length(obj);
            
            for i=1:nkin
                
                [knch,nsamp] = size(obj(i).data);
                
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
                        obj(i).data = filtfilt(B,A,obj(i).data')';
                        
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
                        
                    case 'average'
                        N  = opt.par(1);
                        n  = opt.par(2);
                        nresamp = floor(nsamp/n);
                        datatemp = zeros(nch,nresamp);
                        for j=1:nresamp
                            if ~rem(N,2) % even
                                N = N+1;
                            end
                            ind = [(1-N)/2+1:(N+1)/2];
                            
                            datatemp(:,j) = mean(abs(obj(i).data(:,min(nsamp,max(1,n*(j-1)+ind)))),2);
                            
                        end
                        obj(i).data = datatemp;
                        obj(i).time = obj(i).time(1:n:nresamp*n);
                        
                    case 'cliptozero'
                        obj(i).data = obj(i).data .* (obj(i).data>0);
                        
                        
                    otherwise
                        warning('Unknown filter type! kins were not filtered.')
                        
                end
                
            end
            
        end
        
        %------------------------------------------------------------------
        function opt = getDefFilterOpt(obj)
            
            opt.type = 'butter';
            opt.par  = [2 10/(100/2)]; %2 order 10Hz @ 100Hz sampling rate
            
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
        function obj = average(obj,opt)
            % average kin data across trials
            
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
            
            tav = [trange(1):ts:trange(2)];  % times of averaged kins
            nsampav = length(tav);  % number of samples in average
            
            for i=1:ngr
                
                datatemp = zeros(nch,nsampav);
                ndata    = zeros(1,nsampav);
                
                chlabels = obj(gr{i}(1)).chlabels;
                kavtemp.postime = tav;
                
                for j=1:length(gr{i})
                    jj = gr{i}(j);
                    ktime = obj(jj).time - tref(jj);
                    ind = find(ktime>=tav(1) & ktime<=tav(end)); % find samples contained into chosen interval
                    if ~isempty(ind)
                        if isequal(nsampav,length(ind))
                            indav = [1:nsampav];
                        else
                            [mm,im] = min(abs(tav-ktime(ind(1)))); % find sample closest to first data sample
                            indav = im + [1:length(ind)] - 1;
                        end
                        
                        kdata = obj(jj).data(:,ind);
                        if all(~isnan(kdata(:)))
                            data_j = interp1(ktime(ind),obj(jj).data(:,ind)',tav(indav),'spline')';
                            
                            
                            datatemp(:,indav) = datatemp(:,indav) + data_j;
                            
                            %                         to avoid discontinuity at boundaries
                            indav_pre = find(tav<tav(indav(1))); % indexes of samples before averaging interval
                            indav_post = find(tav>tav(indav(end))); % indexes of samples before averaging interval
                            % set data ouside averaging interval to boundary value
                            if ~isempty(indav_pre), datatemp(:,indav_pre) = datatemp(:,indav_pre) + data_j(:,1)*ones(1,length(indav_pre)); end
                            if ~isempty(indav_post), datatemp(:,indav_post) = datatemp(:,indav_post) + data_j(:,end)*ones(1,length(indav_post)); end
                            ndata = ndata + 1;
                            
                        end
                    end
                end
                
                % mean
                datam = datatemp ./ (ones(nch,1)*(ndata + (ndata==0)));
                kavtemp.pos = datam;
                
                kav(i) = KinData(kavtemp,i,chlabels);
                
            end
            
            obj = kav;
            
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
            %normalize data in kin amplitude
            %
            %   type  action
            %   -------------------------------------------------------------
            %   0         use normdata [nch,1] for normalization of each channel
            %   1         normalize to max of any channel
            %   2         normalize each channel to max in that channel
            %   32        max of each abs(channel)
            
            % set options or use defaults
            defopt = obj.getDefNormalizeOpt;
            if nargin>1 && isstruct(opt)
                fname = fieldnames(opt);
                for i=1:length(fname)
                    defopt.(fname{i}) = opt.(fname{i});
                end
            end
            opt = defopt;
            
            nkin = length(obj);
            
            opt.isect = intersect([1:nkin],opt.isect);
            if isempty(opt.isect)
                warning('empty isect, using all sections')
                opt.isect = [1:nkin];
            end
            nch = size(obj(opt.isect(1)).data,1);
            switch opt.type
                case 1
                    normdata = max(max(obj(opt.isect)));
                case 2
                    normdata = max(obj(opt.isect));
                case 32
                    
                    normdata = max(obj(opt.isect),2);
            end
            
            for i=1:nkin
                
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
            
            nkin = length(obj);
            
            opt.type = 2;  % max absolute value of each channels
            opt.isect = [1:nkin]; % sections to use for computing max
            opt.normdata = [];
            
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
            % plot KIN data
            
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
                yl = kinlim(obj(opt.isect));
            else
                yl = opt.ylim;
            end
            
            
            %
            % loop on trials
            %
            for ii=1:nsect
                i = opt.isect(ii);
                
                
                
                tt.data = obj(i).data(opt.kinsel,:);
                tt.time = obj(i).time;
                if ~isempty(opt.tref) && length(opt.tref)==nkin, tt.time = tt.time - opt.tref(i); end
                t(ii) = Traces(tt);
                if ii==1, t(ii).chlabels = obj(i).chlabels(opt.kinsel); end
                if opt.usetitle
                    if isempty(opt.kintitle) || isempty(opt.kintitle{i})
                        t(ii).label = sprintf('%i',obj(i).trialId);
                    else
                        t(ii).label = sprintf('%s (%i)',opt.kintitle{i},obj(i).trialId);
                    end
                end
                if ~isempty(opt.xlim)
                    if isequal(size(opt.xlim,1),nkin)
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
                    t(ii).opt.yscale = opt.kinscale;
                    t(ii).opt.yscalelabel = opt.kinscalelabel;
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
            
            nkin = length(obj);
            nch = size(obj(1).data,1); % number of channels
            
            %order
            opt.order = 0;
            
            % kinsel
            opt.kinsel = [1:nch];
            
            % isect
            opt.isect = [1:8];
            
            % tref
            opt.tref = [];
            
            % events
            opt.event_code = {};
            opt.event_time = {};
            opt.events_color = {};
            opt.events_style = {};
            
            % title
            opt.usetitle = 1;
            opt.kintitle = {};
            
            % plotting options
            opt.figure = [];
            opt.axes   = [];
            opt.pos = [.09 .08 .88 .86];
            opt.overlap = 0;
            opt.spacing = .01;
            
            opt.xlim = [];
            opt.ylim = []; % for each individual trace
            opt.kinscale = 0;
            opt.kinscalelabel = '';
            opt.fill = 0;
            opt.color = 'k';
            opt.linewidth = [];
            
        end
        
        %------------------------------------------------------------------
        function du = duration(obj)
            % compute duration of KIN data for each trials
            nkin = length(obj);
            for i=1:nkin
                du(i) = range(obj(i).time);
            end
            
        end
        
        %------------------------------------------------------------------
        function val = isequalinterval(obj)
            % check if trials have the same time samples
            nkin = length(obj);
            val = true;
            timeref = obj(1).time;
            for i=2:nkin
                if ~isequal(obj(i).time,timeref)
                    val = false;
                    return
                end
            end
            
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
            %get time range of KIN data of each trial
            
            nkin = length(obj);
            for i=1:nkin
                tr(i,:) = obj(i).time([1 end]);
            end
            
        end
        %------------------------------------------------------------------
        function val = kinlim(obj)
            nkin = length(obj);
            for i=1:nkin
                valmin(i) = min(obj(i).data(:));
                valmax(i) = max(obj(i).data(:));
            end
            val = [min(valmin) max(valmax)];
        end
        
        
        %------------------------------------------------------------------
        function val = max(obj,type)
            % compute max of each channel
            %
            % type = 0 => max of all trials
            % type = 1 => max of each trial
            % type = 2 => max of all trials abs
            
            if nargin<2, type = 0; end
            
            nkin = length(obj);
            
            for i=1:nkin
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
        
        %------------------------------------------------------------------
        
        function v = velocity(obj,t)
            % return velocity at time t
            
            [nkin,nsect] = size(obj);  % nkin -> number of kin objects per section; nsect -> number of sections
            
            tr = timerange(obj);
            dt = tsamp(obj);
            
            if nargin<2 | ~isequal(size(t,1),nsect)
                t = tr;
                
            end
            
            for i=1:nsect
                for ikin=1:nkin
                    if isempty(obj(i))
                        v(ikin,i).data = [NaN NaN NaN]';
                        v(ikin,i).time = NaN;
                    else
                        % data
                        ndof = size(obj(i).data,1);
                        
                        kdata = obj(i).data;
                        
                        ktime = obj(i).time;
                        
                        if all(tr(i,1)>t(i,:)) | all(tr(i,2)<t(i,:))
                            warning(sprintf('time samples out of time range for sect %i kin %i',i,ikin))
                            v(ikin,i).data = NaN;
                            v(ikin,i).time = NaN;
                        else
                            ind = find(t(i,:)>=tr(i,1) & t(i,:)<=tr(i,2));
                            timei = t(i,ind);
                            
                            dtime =  diff(ktime,1,2);
                            kdata = diff(kdata,1,2)./(ones(ndof,1)*dtime);  % along rows
                            ktime = ktime(1:end-1) + dtime/2;
                            for kk =1:ndof
                                pp = csaps(ktime, kdata(kk,:));
                                v(ikin,i).data(kk,:) = fnval(pp,timei);
                                v(ikin,i).time = timei;
                            end
                            
                        end
                    end
                end
            end
            
            
        end
        
        %------------------------------------------------------------------
        
        function a = acceleration(obj,t)
            % return velocity at time t
            
            [nkin,nsect] = size(obj);  % nkin -> number of kin objects per section; nsect -> number of sections
            
            tr = timerange(obj);
            dt = tsamp(obj);
            
            if nargin<2 | ~isequal(size(t,1),nsect)
                t = tr;
            end
            
            for i=1:nsect
                for ikin=1:nkin
                    
                    % data
                    ndof = size(obj(i).data,1);
                    
                    kdata = obj(i).data;
                    
                    ktime = obj(i).time;
                    
                    if all(tr(i,1)>t(i,:)) | all(tr(i,2)<t(i,:))
                        warning(sprintf('time samples out of time range for sect %i kin %i',i,ikin))
                        
                    else
                        % acceleration
                        for j=1:2
                            dtime =  diff(ktime,1,2);
                            kdata = diff(kdata,1,2)./(ones(ndof,1)*dtime);  % along rows
                            ktime = ktime(1:end-1) + dtime/2;
                        end
                        
                        % find valid time range and interpolate
                        ind = find(t(i,:)>=tr(i,1) & t(i,:)<=tr(i,2));
                        timei = t(i,ind);
                        
                        for kk =1:ndof
                            pp = csaps(ktime, kdata(kk,:));
                            a(ikin,i).data(kk,:) = fnval(pp,timei);
                            a(ikin,i).time = timei;
                        end
                        
                    end
                    
                end
            end
            
        end
        
        %------------------------------------------------------------------
        
        function obj = resample(obj,timenew)
            %resample kin object at timenew times
            
            
            if nargin<1
                return
            end
            
            [nkin,nsect] = size(obj);  % nkin -> number of kin objects per section; nsect -> number of sections
            
            [ntime,nsamp] = size(timenew);
            if ntime==1 & nsamp==1  % single value -> timenew is the resampling period, use full range
                tr = timerange(obj);
                freq = 1/timenew;
                timenew = [min(tr(:,1)):1/freq:max(tr(:,2))];
                timenew = ones(nsect,1)*timenew;
            else
                if ~isequal(ntime,nsect)  % use first row of timenew for all sections
                    timenew = ones(nsect,1)*timenew(1,:);
                end
            end
            
            for i=1:nsect
                for ikin=1:nkin
                    
                    tr = timerange(obj(ikin,i));
                    if all(tr(1)>timenew(i,:)) | all(tr(2)<timenew(i,:))
                        warning(sprintf('new times out of time range for sect %i kin %i',i,ikin))
                    else
                        
                        % find valid range
                        ind = find(timenew(i,:)>=tr(1) & timenew(i,:)<=tr(2));
                        tresamp = timenew(i,ind);
                        datafit = spline(obj(ikin,i).time,obj(ikin,i).data);
                        obj(ikin,i).time = tresamp;
                        
                        obj(ikin,i).data = fnval(datafit,tresamp);
                        
                    end
                    
                end
            end
            
        end
        
        %------------------------------------------------------------------
        
        function obj = shift(obj,tshift)
            % shift time base by tshift (t <- t + tshift)
            
            
            if nargin<2 | isempty(tshift), return, end
            
            nkin = length(obj);
            
            if ~isequal(length(tshift),nkin)  % only one value specified
                tshift = ones(nkin,1)*tshift;
            end
            
            for i=1:nkin
                
                if ~isempty(obj(i))
                    
                    obj(i).time = obj(i).time + tshift(i);
                    
                    
                    
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
