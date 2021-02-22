function [whichleave,tsamp] = sampleHopTime_tab(P,tlist,cumdist,nt)
% use tabulated values to sample hopping time from a given node

% uniform variable
w = rand(1);

% decide which node you hop to next, according to P_ik*
deg = length(P);
whichleave = datasample(1:deg,1,'Weights',P);

% use tabulated cumulative distribution to sample
tsamp = interp1(cumdist(1:nt(whichleave),whichleave),tlist(1:nt(whichleave),whichleave),w);
end
