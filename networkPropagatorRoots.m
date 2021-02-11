function nethopinfo = networkPropagatorRoots(NT,options)
% tabulate all the propagator roots for all the network nodes
% also get splitting probability Pvals, and denominator derivatives (invrp)
% to use in sampling leaving times
% returns a structure with the following fields:
% nethopinfo.Pvals = splitting probabilities
% nethopinfo.uroots(:,i) = propagator roots (u_p) for node i
% nethopinfo.nroots(i) = number of roots for node i
% nethopinfo.rpu2(:,:,i) = (r_p)/u_p^2 for each adjacent edge and each
% nethopinfo.tstar = time to switch from short to long time expansions, for
% each node
% root, at node i

opt.epsilon = 1e-14;

if (exist('options','var'))
    opt = copyStruct(options,opt,1);
end


uroots = zeros(500,NT.nnode);
hprime = zeros(max(NT.degrees),500,NT.nnode);
Pvals = zeros(max(NT.degrees),NT.nnode);

for nc = 1:NT.nnode
    deg= NT.degrees(nc);
    lens = NT.edgelens(NT.nodeedges(nc,1:deg))';
        
    invlens = 1./lens;
    Pvals(1:deg,nc) = invlens/sum(invlens);
    
    [tmproots,tmpinvrp,tmpnroots,tmptstar] = nodePropagatorRoots(lens,opt);
    nroots(nc) = tmpnroots;
    tstar(nc) = tmptstar;
    
    uroots(1:nroots(nc),nc) = tmproots;
    invrp(1:deg,1:nroots(nc),nc) = tmpinvrp;
    invrpu2(1:deg,1:nroots(nc),nc) = tmpinvrp.*tmproots.^2;
end

uroots = uroots(1:max(nroots),:);
invrp = invrp(:,1:max(nroots),:);
invrpu2 = invrpu2(:,1:max(nroots),:);

nethopinfo = struct();
nethopinfo.uroots = uroots;
nethopinfo.nroots = nroots;
nethopinfo.rpu2 = 1./invrpu2;
nethopinfo.Pvals = Pvals;
nethopinfo.tstar = tstar;
nethopinfo.epsilon = opt.epsilon;

%if tabulatetimes
%end
end