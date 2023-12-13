function box=PositionBox(outbox,inbox)
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/MartaRussoPhD/SynergyAnalyzerToolbox.git
%
box(1) = outbox(1)+outbox(3)*inbox(1);
box(2) = outbox(2)+outbox(4)*inbox(2);
box(3) = outbox(3)*inbox(3);
box(4) = outbox(4)*inbox(4);
end
