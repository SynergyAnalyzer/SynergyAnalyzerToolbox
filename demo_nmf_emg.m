% SynergyAnalyzer demo: extraction of spatial synergies using NMF
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
disp('Loading raw EMG data (reaching to 8 target in frontal and sagittal planes)')
load('data.mat')

%% create object and preprocess EMG data

info = getInfo(data);
par.type = 'EmgData';
par.chlabels = emgchannels;

sa = SynergyAnalyzer(data,info,par);
sa.opt.verbose = 1;

% filter and resample
sa.opt.emgFilter.type = 'fir1';
sa.opt.emgFilter.par = [50 .04];  % 20 Hz @ 1KHz EMG sampling rate
sa.opt.emgFilter.resample = 1;
sa.opt.emgFilter.resample_period = .01; % resampling period [s]
saf = sa.dataFilter;

%% plot sample raw and filtered EMG data
disp('Press any key to plot a sample of raw EMG data')
pause

for i=1:8
  itrial = sa.findTrials('type',[1 1 i]); % [plane center-out direction]
  ind(i) = itrial(1);
end
optemg.figure = figure('NumberTitle','off','Name','Raw EMG data');

plot(sa.data(ind),optemg) % frontal plane, 8 directions

disp('Press any key to plot a sample of filtered EMG data')
pause
optemg.fill = 1;
optemg.figure = figure('NumberTitle','off','Name','Filtered EMG data');

plot(saf.data(ind),optemg) % frontal plane, 8 directions

%% get phasic activity by subtracting tonic activity
saf = saf.emgPhasic;

%% average data in groups of trials with same plane, start condition, direction
saf.opt.average.gr = saf.groupTrials('type3',[1:8]');
saf.opt.average.trange = [-.3 .7];
sav = saf.average;

% normalize amplitude
sav.opt.normalize.type = 2;
sav = sav.normalize;

%% plot sample filtered, averaged, and normalized EMG data
disp('Press any key to plot a sample of filtered, averaged, and normalized EMG data')
pause
optemg.figure = figure('NumberTitle','off','Name','Averaged EMG data');
plot(sav.data(1:8),optemg) % frontal plane, 8 directions


%% extract spatial synergies
sav.opt.find.N = [1:8];
sav.opt.find.nrep = 5; % number of repetitions
sav.opt.find.niter = [5 5 1e-4 100]; % number of iterations or termination condition
sav.opt.find.plot = 0;
s1 = sav.find;

%% plot spatial synergies
disp('Press any key to plot extracted spatial synergies')
pause

s1.opt.plot.type = 'rsq';
plot(s1)

%%
s1.opt.plot.N = input('Enter the number of synergies\n');
cols = jet(s1.opt.plot.N);
for i=1:s1.opt.plot.N
    colors{i} = cols(i,:);
end
s1.opt.plot.syncolor = colors(1:s1.opt.plot.N);
s1.opt.plot.type = 'W';
plot(s1)

%%
disp('Press any key to plot EMG data reconstruction by spatial synergies')
pause
s1.opt.plot.type = 'rec';
s1.opt.plot.isect = [1:8];
plot(s1)

