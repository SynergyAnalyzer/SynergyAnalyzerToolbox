function box=ArrayBox(outbox,nrow,ncol,row,col,BSp)
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/MartaRussoPhD/SynergyAnalyzerToolbox.git
%

if nargin<6, BSp  = .01; end
if length(BSp)>1
  SpH = BSp(1);
  SpV = BSp(2);
else
  SpH = BSp(1);
  SpV = BSp(1);
end  

BWi  = (1-(ncol+1)*SpH)/ncol;
BHe  = (1-(nrow+1)*SpV)/nrow;
box  = PositionBox(outbox,[SpH+(col(1)-1)*(BWi+SpH),...
    1-row(end)*(BHe+SpV),...
    BWi*(col(end)-col(1)+1),...
    BHe*(row(end)-row(1)+1)]);
end
