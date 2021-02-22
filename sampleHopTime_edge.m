function [whichleave,tsamp,success,edgehopinfo] = sampleHopTime_edge(x0,L,xi)
% for a particle starting on an edge, decide which node to hop to
% then sample time for it to first hit that node.
% x0 = starting position along the edge (from 0 to L)
% L = total length of edge

% use the general sampleHopTime routine for sampling, while passing the
% correct splitting probabilities, roots, residues for the degree-2 case

% xi = -log(epsilon) where epsilon is the tolerance, used for setting
% transition between short time and long time limits

% assumes D = 1

% returns:
% whichleave = which node was hit first (1 for node at 0, 2 for node at L)
% tsamp = sampled leaving time
% success = successfully solved for leaving time?
% edgehopinfo = structure with fields uroots, rpu2, lens (for later use in
% sampling position along edge after a certain time interval)

if (~exist('xi')) % set default
    xi = - log(1e-14);
end

Lx0 = L-x0;
minLx0 = min(x0,Lx0);

if (minLx0 < 1e-8) % too close to edge, will not sample
    tsamp = 0;
    P = [1-x0/L,x0/L];
    whichleave = datasample([1,2],1,'Weights',P);
    success = true;
    edgehopinfo = struct('uroots',[],'rpu2',[],'lens',[],'tstar',[]);
    return
end
    

% time below which to switch to short-t limit
tstar = 25*minLx0^2/(4*xi);
%%

% get the appropriate roots and residues for starting at this position
% on this edge

% assume we will not be sampling t below epsilon precision
mmax = min(ceil(2*xi*L/(5*pi*min(x0,Lx0))), xi^2*L/sqrt(eps)/pi);
mlist = 1:mmax;

uroots = (mlist*pi/L)';
uroots2 = uroots.^2;

sinx0 = sin(uroots'*x0);
rpu2  = 2/pi*sinx0./mlist;
rpu2 = [rpu2; -rpu2.*(-1).^mlist];

P = [1-x0/L,x0/L];

lens = [x0,Lx0];
[whichleave,tsamp, success] = sampleHopTime(P,uroots2,rpu2,lens,tstar);

edgehopinfo = struct('uroots',uroots,'rpu2',rpu2,'lens',lens,'tstar',tstar);

end
