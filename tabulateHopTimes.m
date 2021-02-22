function nethopinfo = tabulateHopTimes(nethopinfo,NT)

% tabulate cumulative distribution for conditional hopping time 
% to each neighbor from each node
%%
opt = struct()
% size of central region of times around l^2/2 
opt.rangescl = 10; % scaling factor for how far out to go in leaving times
% number of points to tabulate over, logarithmically spaced
opt.nt = 1000; 

% outer tabulation region, with fewer points
opt.rangesclouter = 100; 
opt.ntouter = 100; 
 
if (exist('options','var'))
    opt = copyStruct(options,opt);
end

nt = opt.nt + 2*opt.ntouter-2;
hopTimeCumDist = zeros(nt,max(NT.degrees),NT.nnode);
hoptimelist = hopTimeCumDist;
nhoptime = zeros(max(NT.degrees),NT.nnode);

for nc = 1:NT.nnode
    deg = NT.degrees(nc);
    
    lens = NT.edgelens(NT.nodeedges(nc,1:deg))';
    
    nr = nethopinfo.nroots(nc);
    P = nethopinfo.Pvals(1:deg,nc);
    uroots2 = nethopinfo.uroots(1:nr,nc).^2;                        
    rpu2 = nethopinfo.rpu2(1:deg,1:nr,nc);
          
    % min and max values for central range with lots of sample points
    tmin1 = min(lens)^2/2/opt.rangescl; tmax1 = max(lens)^2/2*opt.rangescl;
    % min and max values for outer range with fewer sample points
    tmin2 = min(lens)^2/2/opt.rangesclouter; tmax2 = max(lens)^2/2*opt.rangesclouter;
    tlist1 = logspace(log10(tmin1),log10(tmax1),opt.nt);
    if (opt.ntouter==0)
        tlist = tlist1;
    else
        tlist2 = logspace(log10(tmin2),log10(tmin1),opt.ntouter);
        tlist3 = logspace(log10(tmax1),log10(tmax2),opt.ntouter);
        tlist = [tlist2(1:end-1) tlist1 tlist3(2:end)];
    end       
       
    smallind = tlist<nethopinfo.tstar(nc);
    largeind = tlist>=nethopinfo.tstar(nc);
    smallt = tlist(smallind);
    larget = tlist(largeind);
    
    for whichleave = 1:deg
        
        lk= lens(whichleave);
        Pk = P(whichleave);
        
        % large t series for survival probability, for times above t*
        largetfunc = @(t) rpu2(whichleave,:)*exp(-uroots2*t)/Pk;

        % small t approximation for survival probability, for times below t*
        smalltfunc = @(t) 1 - ((2/deg)*(erfc(lk./(2*sqrt(t))) + erfc(3*lk./(2*sqrt(t))))...
         - 4*sum(erfc((lk+2*lens')*1./(2*sqrt(t))),1)./deg^2)/Pk;      
       
       smallvals = smalltfunc(smallt)';
       largevals = largetfunc(larget)';
       allvals = [smallvals; largevals];
       allt = [smallt larget];
       
       % get rid of redundancies and set edge values
       % first deal with the low t side
       % last instance of non-decreasing values or value above 1
       ind = find(diff(allvals)>=0 | allvals(1:end-1)>1,1,'last');       
       %ind = max(ind,find(allvals>1,1,'last'));
       if (~isempty(ind))            
           allt = allt(ind+1:end);
           allvals(ind+1) = 1;
           allvals = allvals(ind+1:end);           
       end
     
       % now deal with the high t side
       ind = find(allvals<=0,1);
       if (isempty(ind) |ind==length(largevals) )
           allvals(end) = 0; % force tmax2 to be smallest t value ever sampled
       else           
           allt = larget(1:ind);
           allvals = allvals(1:ind);
           allvals(ind) = 0;           
       end
                 
       nvals = length(allvals);
       
       hopTimeCumDist(1:nvals,whichleave,nc) = allvals;              
       hoptimelist(1:nvals,whichleave,nc) = allt';
       nhoptime(whichleave,nc) = nvals;
    end
end

% adjust outer edges to always sample between min and max time
hopTimeCumDist(1,:) = 1;
hopTimeCumDist(end,:) = 0; 

% save results into nethopinfo structure
nethopinfo.hoptimecumdist = hopTimeCumDist;
nethopinfo.hoptimelist = hoptimelist;
nethopinfo.nhoptime = nhoptime;

end