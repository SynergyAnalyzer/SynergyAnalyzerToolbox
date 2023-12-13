function info = getInfo(data)
% get trial type information contained into data
%
% for reaching data assumes data.info has plane, start, and target codes
%  plane: 1 -> frontal, 2 -> sagittal
%  start, target: 0 -> center, 1 -> medial/back, 3 -> down, 5 ->
%                 lateral/forward, 7-> up
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git
% License: GNU GPL v3
%

ntrial = length(data);


for i=1:ntrial
    % id
    info(i).id = i;
    
    % type
    info(i).type(1) = data(i).info.plane;   % plane
    
    if data(i).info.start==0   % center-out/out-center, direction
        info(i).type(2) = 1;
        info(i).type(3) = data(i).info.target;
    else
        info(i).type(2) = 0;
        info(i).type(3) = mod(data(i).info.start+3,8)+1;  % assumes 8 targets
    end
    
    % selected
    info(i).selected = 1;
    
    % events
    info(i).events.code = [13 14];
    info(i).events.time = [data(i).info.t_onset data(i).info.t_end];
    
    
end

end