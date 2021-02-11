function [foundroots,invrp,nroot,tstar,poles] = nodePropagatorRoots(lens,options)
% for a particular node neighborhood, calculate the roots of the propagator
% function denominator
% lens = row vector of edge lengths in the neighborhood of the node
% options = structure for optional parameters
% assumes D=1
% output:
% foundroots = u values for roots (positive real numbers)
% invrp = inverse residues (1/r_p) evaluated at the roots
% tstar = time cutoff to be used in switching from short-time to long-time
% limit
% poles = poles used as boundaries to search for u values

opt = struct();
% fractional shift of boundaries between which we look for roots
opt.boundshift = 1e-3;
opt.minsep = 1e-10; % minimal separation between boundaries to look for roots
opt.epsilon = 1e-14; % convergence tolerance for determining max number of roots to use

if (exist('options','var'))
    opt = copyStruct(options,opt);
end

% max prefactor for calculating roots
xi = -log(opt.epsilon);
mmax = 2*xi*max(lens)/(5*pi*min(lens));

xi = -log(opt.epsilon);
tstar = 25*min(lens)^2/(4*xi);

%%
% find the poles of the *denominator* for the Pik. Roots of the denominator will be in between
% these


poles = [];
lmin = min(lens);
for ec = 1:length(lens)
    mmax = 2*xi/(5*pi*lmin)*lens(ec);
    mvals = 1:mmax;
    poles = [poles mvals*pi/lens(ec)];
end
poles = sort(poles);

% remove poles that are too close together
diffs = diff(poles);
while (min(diffs)<opt.minsep)
    poles(diffs<opt.minsep) = [];
    diffs = diff(poles);
end
poles = [0 poles];

%% look for a root between every pair of poles
% function for finding roots
rootfunc = @(u) sum(cot(lens'*u),1);

for uc = 1:length(poles)-1
    del = opt.boundshift*(poles(uc+1)-poles(uc));
    foundroots(uc) = fzero(rootfunc,[poles(uc)+del poles(uc+1)-del]);
end

nroot = length(foundroots);
%% calculate the inverse residues (denominator derivatives) at the denominator roots
UL = lens'*foundroots;


invrp = 0.5*sin(UL)./foundroots.*(lens*csc(UL).^2);
