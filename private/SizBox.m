function sbox = SizBox(outbox,sizHe,sizWi,irow,icol,Sp) 
% find position box depending on size
% 
%   sbox = SizBox(outbox,sizHe,sizWi,irow,icol,Sp)
%
% Synergy Analyzer Toolbox for MATLAB: https://github.com/MartaRussoPhD/SynergyAnalyzerToolbox.git
%

if nargin<6, Sp = .05; end
nrow=size(sizHe,1);
ncol=size(sizWi,2);

if length(Sp)>1
  SpH = Sp(1);
  SpV = Sp(2);
else
  SpH = Sp(1);
  SpV = Sp(1);
end  

maxHe = max(sizHe,[],2);
maxWi = max(sizWi,[],1);
maxHe = maxHe/sum(maxHe);
maxWi = maxWi/sum(maxWi);
BWi  = (1-(ncol+1)*SpH)*maxWi(icol);
BHe  = (1-(nrow+1)*SpV)*maxHe(irow);
if BWi<0 
  warning('Too many columns!')
  BWi  = 1/ncol;
end  
if BHe<0
  warning('Too many rows!')
  BHe  = 1/nrow;
end

cumWi = [0,(1-(ncol+1)*SpH)*cumsum(maxWi)];
cumHe = (1-(nrow+1)*SpV)*cumsum(maxHe);

pos(1) = icol*SpH+cumWi(icol);
pos(2) = 1-irow*SpV-cumHe(irow);
pos(3) = BWi;
pos(4) = BHe;
sbox  = PositionBox(outbox,pos);
end

