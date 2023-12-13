% SynergyAnalyzer demo: extraction of mixed kinematic-muscular spatial synergies using MMF
%
% EMG and kinematic data from fast reaching movement in the frontal plane,
% 1 start conditions (center), 8 directions; 10
% repetitions for each condition; time is referred to movement onset
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%

%% check Matlab version (for compatibility with new Matlab class-definition syntax)
if verLessThan('matlab','7.6')
    error('Matlab 7.6 or later releases are required for running this demo')
end

%% load data from file
disp('Loading raw EMG data & kinematic data (reaching to 8 target in frontal and sagittal planes)')
load('data.mat')

%% create object EMG/KIN data

info = getInfo(data);
par.type = 'EmgKinData';
par.chlabels = [emgchannels(:)', kinchannels(:)'];
par.order = 2; % determines the order of the kinematic derivatives (0->position, 1->velocity, 2->acceleration) 
par.delay = -0.05; % 50ms -> KIN data were anticipated of 50 ms to account for an estimated electromechanical delay between EMG and KIN (Scano et al. 2022)

sa = SynergyAnalyzer(data,info,par);
% sa = SynergyAnalyzer(data,'EmgKinData',[emgchannels(:)', kinchannels(:)'],order,delay);
sa.opt.verbose = 1;


%% preprocess EMG & KIN data

% filter and resample
sa.opt.emgFilter.type = 'fir1';
sa.opt.emgFilter.par = [50 20/(1000/2)];  % 20 Hz @ 1KHz EMG sampling rate
sa.opt.emgFilter.resample = 1;
sa.opt.emgFilter.resample_period = .01; % resampling period [s]

sa.opt.kinFilter.type = 'butter';
sa.opt.kinFilter.par  = [2 10/(100/2)]; % 2nd order 10Hz @ 100Hz sampling rate
saf = sa.dataFilter;

%% plot sample filtered EMG data
ntrial = 8;
for i=1:ntrial
    itrial = sa.findTrials('type',[1 1 i]); % [plane  center-out direction]
    ind(i) = itrial(1);
end

disp('Press any key to plot a sample of aligned and filtered EMG and KIN data')
pause
opt.fill = 1;
opt.figure = figure('NumberTitle','off','Name','Filtered EMG and KIN data');
plot(saf.data(ind),opt) % frontal plane, 8 directions

%% get phasic activity by subtracting tonic activity
saf = saf.emgPhasic;

%% average data in groups of trials with same plane, start condition, direction
saf.opt.average.gr = saf.groupTrials('type3',[1:8]' );  % type3 is the target number  
saf.opt.average.trange = [-.3 .7]; % time interval in s before and after onset for averaging
sav = saf.average;

disp('Press any key to plot a sample of averaged EMG and KIN data')
pause
opt.figure = figure('NumberTitle','off','Name','Averaged EMG and KIN data');
plot(sav.data(1:8),opt) % frontal plane, 8 directions

%% normalize amplitude
disp('Press any key to normalize data for extracting synergies')
disp('NB: once normalized data will be modified in the main structure')
pause
sav.opt.normalize.type = 32;
sav.opt.normalize.nonnegch = sav.data(1).nonnegch;
sav = sav.normalize;

%% plot sample filtered, averaged, and normalized EMG data
opt.figure = figure('NumberTitle','off','Name','Normalized EMG and KIN data');
plot(sav.data(1:8),opt) % frontal plane, 8 directions


%% extract spatial synergies
disp('Press any key to extract')
pause
sav.opt.find.algo = 'mmf';
sav.opt.find.N = [1:14];
sav.opt.find.nrep = 10;
sav.opt.find.niter = [100 10 10^-6 10000];
sav.opt.find.plot = 0;
s1 = sav.find;

%% plot R^2
disp('Press any key to plot extracted spatial synergies')
pause
s1.opt.plot.type = 'rsq';
plot(s1)

%% plot synergies
s1.opt.plot.N = input('Enter the number of synergies\n');
cols = jet(s1.opt.plot.N);
for i=1:s1.opt.plot.N
    colors{i} = cols(i,:);
end
nemg = nonnegch(sav);
nemgkin = nch(sav);

s1.opt.plot.syncolor = colors(1:s1.opt.plot.N);
for ic = 1:nemg, s1.opt.plot.chcolor{ic} = 1; end
for ic = (nemg+1):nemgkin, s1.opt.plot.chcolor{ic} = .75; end
s1.opt.plot.type = 'W';
plot(s1)

%% plot data recontruction
disp('Press any key to plot EMG data reconstruction by spatial synergies')
pause
s1.opt.plot.type = 'rec';
s1.opt.plot.isect = [1:8];
plot(s1)


