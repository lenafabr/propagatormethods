function [MFPTs,Vars,P,Q,dPds,dQds] = networkMFPTanalytic(NT,absrates,targets)
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
% MFPTs = vector of size [NT.nnode, 1] giving the mean first passage time
% from each starting node
% Vars = vector of size [NT.nnode, 1] giving the variance in first passage time
% from each starting node
% ------------
% for debugging only:
% P = matrix of splitting probabilities into each adjacent node
% Q = matrix of mean waiting times at each node
% dPds = derivative of P(s) wrt to s at s =0
% dQds = derivative of Q(s) wrt to s at s =0

nottargets=1:NT.nnode;
nottargets(targets)=[];

% square roots of gammas
sqrtgam = sqrt(absrates);

P = zeros(NT.nnode,NT.nnode);
Q = zeros(NT.nnode,1);
dPds = zeros(NT.nnode,NT.nnode);
dQds = zeros(NT.nnode,1);

% Loop through all nodes to construct:
% P_{nc,j}, Q_nc, and their derivatives for node nc and neighbors=j
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
    
    % Save to P, Q, dPds, and dQds
    P(nc,NT.nodenodes(nc,1:deg)) = numP/sumden;
    Q(nc) = sum(numQ)/sumden;
    dPds(nc,NT.nodenodes(nc,1:deg)) = (1/2)*(dPds1/sumden - (numP/sumden^2)*sum(dPds2));
    dQds(nc) = (1/2)*(sum(dQds1)/sumden - (sum(numQ)/sumden^2)*sum(dPds2));
end

% Mean first passage time is given by (I-P)^(-1).Q
% Mean of square is -2(I-P)^(-1).dPds.(I-P)^(-1).Q - 2(I-P)^(-1).dQds

% initialize variables
MFPTs = zeros(NT.nnode,1);
SquMFPTs = zeros(NT.nnode,1);

% Removing target nodes from matrix inversion--they cause a singularity
if targets
    P = P(nottargets,nottargets);
    Q = Q(nottargets);
    dPds = dPds(nottargets,nottargets);
    dQds = dQds(nottargets);
    
    % inverse without targets
    IPinv = inv(eye(length(nottargets))-P);
else
    % inverse with targets
    IPinv = inv(eye(NT.nnode)-P);
end

if targets
    MFPTs(nottargets)=IPinv*Q;
    SquMFPTs(nottargets)=-2*(IPinv*dPds*IPinv*Q + IPinv*dQds);
else
    MFPTs=IPinv*Q;
    SquMFPTs=-2*(IPinv*dPds*IPinv*Q + IPinv*dQds);
end

Vars = SquMFPTs - MFPTs.^2;

end