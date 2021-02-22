function [whichedge,xsamp,success] = sampleHopPosition(lens,uroots,rpu2,dt,shortlim)
% given particle starts on a node, sample its position along one of the
% adjacent edges, given it has not passed any other nodes
% use erfc approximation for dt<t* and series over roots for dt>t*
% uroots= precalculated propagator roots for the large t series
% rpu2 = precalculated r_p/u_p^2 for the large t series
% lens = row vector of edge lengths attached to this node
% shortlim = use short-time limit?
% assumes D = 1;

if (isempty(lens))
    error('bad lengths in sampleHopPosition')
elseif abs(dt) <= 0
    error('bad dt in sampleHopePosition')
end

%%
if (shortlim) % short-time limit
    
    deg = length(lens);
    sdt = sqrt(dt);
    erfcvals = erfc(lens/sdt);
    
    sumj= 2/deg^2*sum(erfcvals);
    
    % weight of ending up on each edge
    Y0 = 1/deg*(1 + 2*erfcvals - 2*erfc(lens/(2*sdt)) - 2*erfc(3*lens/(2*sdt))) ...
        - sumj + 4/deg^2*sum(erfc((2*lens' + lens)/(2*sdt)),1) - 2/deg^2*sum(erfc((lens'+lens)/sdt),1);
    % decide what edge you end up on
    whichedge = datasample(1:length(lens),1,'Weights',Y0);
    lik = lens(whichedge);
    
    % sample position along that edge
    w = rand(1);
    
    constterms = 1/deg*(1 + 2*erfcvals(whichedge)) - sumj - 2/deg^2*sum(erfc((lens+lik)/sdt));
    func = @(x) (constterms - 1/deg*(erfc((lik-x)/(2*sdt)) + erfc((lik+x)/(2*sdt)) ...
        + erfc((3*lik-x)/(2*sdt)) + erfc((3*lik+x)/(2*sdt)))...
        + 2/deg^2*sum(erfc((2*lens' + lik - x)/(2*sdt))+ erfc((2*lens' + lik + x)/(2*sdt)),1))/Y0(whichedge) - w;
    
    
    [xsamp,fval,exitflag] = fzero(@(x) func(x),[0,lens(whichedge)]);
    success = (exitflag==1);
    
else % use long-time asymptotic series
    UL = lens'*uroots';
    cosUL = cos(UL);
    expvals = exp(-uroots.^2*dt);
    
    % Un-normalized probability of ending step on each edge:
    Y0 = rpu2.*(1-cosUL)*expvals;
    
    % decide what edge you end up on
    whichedge = datasample(1:length(lens),1,'Weights',Y0);
    
    % sample position along that edge
    w = rand(1);
    
    coeff = rpu2(whichedge,:).*expvals';
    func = @(x) coeff*(cos(uroots*x) - cosUL(whichedge,:)')/Y0(whichedge) - w;
    
    [xsamp,fval,exitflag] = fzero(@(x) func(x),[0,lens(whichedge)]);
    success = (exitflag==1);
end

end