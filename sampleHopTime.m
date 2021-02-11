function [whichleave,tsamp, success] = sampleHopTime(P,uroots2,rpu2,lens,tstar)
% decide which node to hop to, weighted by P_ik*
% then, sample time to leave node and hit that neighbor
% use erfc approximation for t<t* and series over roots for t>t*
% uroots2= precalculated propagator roots for the large t series, squared
% rpu2 = precalculated r_p/u_p^2 for the large t series
% lens = row vector of edge lengths attached to this node
% tstar = cutoff time to switch between long-time vs short-time limits


opt = struct();
% how many expansions to try for bracketing
opt.maxtry = 10;
% how much to multiply when bracketing
opt.scl = 10;

% uniform variable
w = rand(1);

% make sure it's a row vector
if (size(lens,1)>size(lens,2))
    lens = lens';
end

% decide which node you hop to next, according to P_ik*
deg = length(lens);
whichleave = datasample(1:deg,1,'Weights',P);


lk= lens(whichleave);
Pk = P(whichleave);

% large t series for survival probability, for times above t*
largetfunc = @(t) rpu2(whichleave,:)*exp(-uroots2*t)/Pk - w;

% small t approximation for survival probability, for times below t*
smalltfunc = @(t) (Pk - ((2/deg)*(erfc(lk./(2*sqrt(t))) + erfc(3*lk./(2*sqrt(t))))...
    - 4*sum(erfc((lk+2*lens')*1./(2*sqrt(t))),1)./deg^2))/Pk - w;


% bracket the solution interval
minT = min(lens)^2/10;
maxT = max(lens)^2*2;
if (minT< tstar)
    fvals(1) = smalltfunc(minT);
else
    fvals(1) = largetfunc(minT);
end
if (maxT< tstar)
    fvals(2) = smalltfunc(maxT);
else
    fvals(2) = largetfunc(maxT);
end


for trial = 1:opt.maxtry
    if (fvals(1)>0 & fvals(2)<0); break; end
    if (fvals(1)<0); minT = minT/opt.scl; end
    if (fvals(2)>0); maxT = maxT*opt.scl; end
    if (minT< tstar)
        fvals(1) = smalltfunc(minT);
    else
        fvals(1) = largetfunc(minT);        
    end
    if (maxT< tstar)
        fvals(2) = smalltfunc(maxT);      
    else
        fvals(2) = largetfunc(maxT);
    end 
end

% decide which limit to use to pinpoint root
if (maxT < tstar)
    uselim = 1; % small t limit
elseif (minT>tstar)
    uselim = 2; % large t limit
else
    % spans across the transition
    fvalstar = largetfunc(tstar);
    if (fvalstar>0)
        minT = tstar; uselim = 2;
    else
        maxT = tstar; uselim = 1;
    end
end
    

if trial>=opt.maxtry
    error('failed to bracket function in sampleHopTime')
end

if (uselim==1)
    [tsamp,fval,exitflag] = fzero(@(t) smalltfunc(t),[minT,maxT]);
else
    [tsamp,fval,exitflag] = fzero(@(t) largetfunc(t),[minT,maxT]);
end
success = (exitflag==1);

end