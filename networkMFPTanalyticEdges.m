function [MFPTs,Vars,P,Q,dPds,dQds,PE,QE,dPEds,dQEds] = networkMFPTanalyticEdges(NT,absrates,targets)
% Gives the mfpt and variance in first passage time for diffusion on a 
% network using an analytic approach

% --------------
% Inputs:
% --------------
% NT is a network object (created with NetworkObj)
% targets = list of target nodes (gives MFPT to any of the targets)
% assumes a diffusivity of 1

% absrates = vector of size [NT.nedge, 1], gives absorbance rates on each
% edge (can be 0 if non-absorbing)

% targets = list of perfectly absorbing target nodes (can be empty)

% ------------
% Outputs: 
% MFPTs = vector of size [NT.nedge, 1] giving the mean first passage time
% from each starting edge
% Vars = vector of size [NT.nnode, 1] giving the variance in first passage time
% from each starting edge
% ------------
% for debugging only:
% P = matrix of splitting probabilities into each adjacent node
% Q = matrix of mean waiting times at each node
% PE = matrix of splitting probabilities from edge (row) to bounding nodes (column)
% QE = matrix of mean waiting times at each edge
% dPds = derivative of P(s) wrt to s at s=0
% dQds = derivative of Q(s) wrt to s at s=0
% dPEds = derivative of PE(s) wrt to s at s=0
% dQEds = derivative of QE(s) wrt to s at s=0

% TODO: current system is slow, get rid of the extra calculations and
%       (zeroind) resets. Both here and networkMFPTanalytic.m
% TODO: finding values for each edge is done multiple times, since loop is over
%        nodes. Maybe no way to avoid this.

nottargets=1:NT.nnode;
nottargets(targets)=[];

% square roots of gammas
sqrtgam = sqrt(absrates);

P = zeros(NT.nnode,NT.nnode);
Q = zeros(NT.nnode,1);
dPds = zeros(NT.nnode,NT.nnode);
dQds = zeros(NT.nnode,1);
PE = zeros(NT.nedge,NT.nnode);
QE = zeros(NT.nedge,1);
dPEds = zeros(NT.nedge,NT.nnode);
dQEds = zeros(NT.nedge,1);

% Loop through all nodes to construct:
% P_{nc,j}, Q_nc, and their derivatives for node nc and neighbors=j
% PE_{m,j}, QE_m, and their derivatives for edge m and bounding nodes=j
for nc = 1:NT.nnode
    deg = NT.degrees(nc);                   % degree of node
    edges = NT.nodeedges(nc,1:deg);         % edges adjacent to this node
    lens = NT.edgelens(edges);              % edge lens from this node
    
    alphaLs = sqrtgam(edges).*lens;
    % deal with alpha=0 cases separately
    zeroind = (alphaLs==0);
    
    den = sqrtgam(edges).*coth(alphaLs);
    den(zeroind) = 1./lens(zeroind);
    
    % Find P values, then correct the gamma=0 cases
    sumden = sum(den);
    numP = sqrtgam(edges)./(sinh(alphaLs));
    numP(zeroind) = 1./lens(zeroind);
    
    % Find Q values prior to sum, correcting gamma=0 cases
    numQ = (1./sqrtgam(edges)).*(coth(alphaLs)-csch(alphaLs)); 
    numQ(zeroind) = lens(zeroind)/2;
    
    % Find dPds, similar process
    dPds1 = (1./sqrtgam(edges)).*csch(alphaLs).*(1 - alphaLs.*coth(alphaLs));
    dPds1(zeroind) = -lens(zeroind)/3;
    dPds2 = (1./sqrtgam(edges)).*coth(alphaLs) - lens.*(csch(alphaLs)).^2;
    dPds2(zeroind) = 2*lens(zeroind)/3;
    
    % Find dQds
    dQds1 = -(1./sqrtgam(edges).^3).*(coth(alphaLs) - csch(alphaLs)) + (lens./sqrtgam(edges).^2).*(coth(alphaLs).*csch(alphaLs) - (csch(alphaLs)).^2);
    dQds1(zeroind) = -(lens(zeroind).^3)/12;
    
    % Find PE and QE
    PE1 = (1./lens)./sqrtgam(edges).*tanh(lens/2.*sqrtgam(edges));
    PE1(zeroind) = 1/2;
    QE1 = (1./sqrtgam(edges)).*(1-(2./alphaLs).*tanh(alphaLs/2));
    QE1(zeroind) = (lens(zeroind).^2)/12;
    
    % Find dPEds and dQEds
    dPEds1 = (1./sqrtgam(edges).^2)/2.*( (1-(1./alphaLs).*sinh(alphaLs))./(1+cosh(alphaLs)) );
    dPEds1(zeroind) = -lens(zeroind).^2/24;
    dQEds1 = (1./sqrtgam(edges).^4).*(3./alphaLs.*tanh(alphaLs/2) - 1/2.*sech(alphaLs/2).^2 - 1);
    dQEds1(zeroind) = -lens(zeroind).^4/120;
    
    % Save to P, Q, dPds, dQds, dPEds, dQEds
    P(nc,NT.nodenodes(nc,1:deg)) = numP/sumden;
    Q(nc) = sum(numQ)/sumden;
    dPds(nc,NT.nodenodes(nc,1:deg)) = (1/2)*(dPds1/sumden - (numP/sumden^2)*sum(dPds2));
    dQds(nc) = (1/2)*(sum(dQds1)/sumden - (sum(numQ)/sumden^2)*sum(dPds2));
    PE(edges,nc) = PE1;
    QE(edges) = QE1;
    dPEds(edges,nc) = dPEds1;
    dQEds(edges) = dQEds1;
end

% Mean first passage time is: QE + PE.(I-P)^(-1).Q
% Mean square first passage time is:
% -2*[dQEds + dPEds.(I-P)^(-1).Q + PE.(I-P)^(-1).dPds.(I-P)^(-1).Q 
%     + PE.(I-P)^(-1).dQds]

% Removing target nodes from matrix inversion--they cause a singularity
if targets
    P = P(nottargets,nottargets);
    Q = Q(nottargets);
    dPds = dPds(nottargets,nottargets);
    dQds = dQds(nottargets);
    PE = PE(:,nottargets);
    dPEds = dPEds(:,nottargets);
    
    % inverse without targets
    IPinv = inv(eye(length(nottargets))-P);
else
    % inverse with targets
    IPinv = inv(eye(NT.nnode)-P);
end

MFPTs = QE + PE*(IPinv*Q);
SquMFPTs = -2*(dQEds + dPEds*IPinv*Q + PE*IPinv*dPds*IPinv*Q + PE*IPinv*dQds);

Vars = SquMFPTs - MFPTs.^2;

end