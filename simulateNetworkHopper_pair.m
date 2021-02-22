function [reacttimes] = simulateNetworkHopper_pair(NT,nethopinfo,npart,options)
% simulate *multiple interacting particles* hopping on a network
% only runs one iteration at a time
% NT = network object
% nethopinfo: precalculated data on roots and splitting probabilities for
% the network
% npart = number of particles to simulate
% savepos = npart x dim x ntime array of saved positions for each particle
% savepos(:,1,:) = edge it is on
% savepos(:,2,:) = position along that edge
%%

opt = struct();
%opt.minT = nethopinfo.tmin;
%opt.maxT = 1e3;
opt.maxsteps = 1e5;
opt.maxtime = 1e5;
opt.hitfirst = true; % look at time to hit the first of many targets
opt.printevery = 10;
% find first hitting times to these target nodes
% by default, check all target nodes
opt.targetnodes = NaN;

% by default, start on random nodes
opt.startnodes = NaN;
% start uniformly on edges
opt.startedgeuniform = false;

% save particle positions at these times
opt.savetimes = [];

% use tabulated distribution of node hopping times
opt.usetabulation = true;
% number of hop times to pre-sample for each node (speeds up code)
opt.npresample = 10000; 

% stop simulation when any of the particles have hit
opt.stopfirsthit = false;

% stop evolving a particle once it has hit all nodes
opt.stopwhenallhit = false; 

% display more notifications
opt.verbose = false;

if (exist('options','var'))
    opt = copyStruct(options,opt);
end

if (isnan(opt.targetnodes))
    opt.hitfirst = false;
    opt.targetnodes = 1:NT.nnode;
end
if (isnan(opt.startnodes))
    opt.startnodes =  randsample(NT.nnode,npart,true);
end

disp(['target nodes: ' sprintf(' %f', opt.targetnodes)])

if (length(opt.startnodes)==1)
    nodepos =opt.startnodes*ones(1,npart);
else
    nodepos = opt.startnodes;
end

% current time for each particle
curtime = zeros(1,npart);
done = false(1,npart);

Pvals = nethopinfo.Pvals;
nroots = nethopinfo.nroots;
uroots = nethopinfo.uroots;
rpu2 = nethopinfo.rpu2;
tstar = nethopinfo.tstar;
xi = -log(nethopinfo.epsilon);

on_edge = false(npart,1);

if (opt.startedgeuniform)
    % start particles uniformly on all edges
    on_edge = true(npart,1);
    
    % sample starting edge id (proportional to edge length)
    lens = NT.edgelens;
    edgeid = datasample(1:NT.nedge,npart*2,'Weights',lens);
    edgeid = [edgeid(1:npart); edgeid(npart+1:end)];
    % sample starting edge position uniformly
    edgepos = rand(2,npart).*lens(edgeid);
    nodepos = zeros(2,npart);
else
    error('multi particle sims only set up with startedgeuniform for now')
end

nextsave = ones(npart,1);
if (isempty(opt.savetimes))
    nextsavetime = inf*ones(npart,1);
else
    nextsavetime = opt.savetimes(1)*ones(npart,1);
end

savepos = zeros(npart,2,2,length(opt.savetimes));

%% presample hopping times from all nodes to all neighbors to
% take better advantage of matlab's parallelization
if (opt.npresample > 0 && opt.usetabulation)
    
    presampleInd = zeros(max(NT.degrees),NT.nnode);
    presampleTimes = zeros(opt.npresample,max(NT.degrees),NT.nnode);
    presampleWhichLeaveInd = zeros(1,NT.nnode);
    presampleWhichLeave = zeros(opt.npresample,NT.nnode);
    
    % presample many hopping times from each node to each neighbor
    for nc = 1:NT.nnode
        deg = NT.degrees(nc);
        % presample where to hop from each node
        P = Pvals(1:deg,nc);
        presampleWhichLeave(:,nc) = datasample(1:deg,opt.npresample,'Weights',P);
        presampleWhichLeaveInd(nc) = 1;
        
        for whichleave = 1:NT.degrees(nc)
            nt = nethopinfo.nhoptime(whichleave,nc);
            tlist = nethopinfo.hoptimelist(1:nt,whichleave,nc);
            cumdist = nethopinfo.hoptimecumdist(1:nt,whichleave,nc);
            
            w = rand(opt.npresample,1);
            presampleTimes(:,whichleave,nc) = interp1(cumdist,tlist,w);
            % current indices in the presampled arrays
            presampleInd(whichleave,nc) = 1;
        end
    end
end


%% start by establishing next event times for each domain
eventset = false(npart,1); % next event for this particle already set

% sorted list of next event times  and the particles involved
nexteventtime = inf*ones(2,npart);
whichbound = zeros(2,npart);
nextnodepos = zeros(2,npart);
nextedgepos = zeros(2,npart);
nextedgeid = zeros(2,npart);

% establish if particles share a boundary
for pc = 1:npart % cycle over all trials
    sharebound(:,:,pc) = updateShareBound(NT,nodepos(:,pc),edgepos(:,pc),edgeid(:,pc));
end
%%
clear rectinfo edgeinfoboth
for pc = 1:npart
    psample = [1,2]; % sample both particles to start with
    [tsamp, whichbound(:,pc), nextnodepos(:,pc),nextedgeid(:,pc),nextedgepos(:,pc),...
        edgeinfotmp,rectinfotmp,sameedge(pc)] = ...
        sampleNextEvent(NT,psample,nodepos(:,pc),edgeid(:,pc),edgepos(:,pc),nethopinfo);
    
    rectinfo(pc) = rectinfotmp;
    edgeinfoboth(:,pc) = edgeinfotmp;
    nexteventtime(:,pc) = tsamp;
    
    if whichbound(1,pc) == -1
        done(pc) = true;
        reacttimes(pc) = tsamp;
    end
end

%%
for step = 1:opt.maxsteps
    % propagate forward to next event (directly affected particles)
    if (all(done))
        display('all pairs have reacted')
        break
    end
    
    if (mod(step,opt.printevery)==0)
        pleft = find(~done);
        [a,b] = min(curtime(~done));
        disp([step a pleft(b) nnz(done)])
    end
    
    for pc = 1:npart
        if (done(pc)); continue; end
        
        [mintime,nextpart] = min(nexteventtime(:,pc));
        otherpart = 3-nextpart;
        dt = mintime-curtime(pc);
        curtime(pc) = mintime;
        
        if (sameedge(pc)==1)
            edgepos(:,pc) = nextedgepos(:,pc);
            edgeid(:,pc) = nextedgeid(:,pc);
            nodepos(:,pc) = nextnodepos(:,pc);
            %if (any(nodepos(:,pc)>0)); sameedge(pc)=2; end
        else
            % update position of next propagated particle
            edgepos(nextpart,pc) = nextedgepos(nextpart,pc);
            edgeid(nextpart,pc) = nextedgeid(nextpart,pc);
            nodepos(nextpart,pc) = nextnodepos(nextpart,pc);
        end
        
        if (whichbound(1,pc)==-1)
            if (opt.verbose)
                disp(sprintf('Iter %d Particles hit',pc))
            end
            % reaction event
            done(pc) = true;
            reacttimes(pc) = mintime;
            continue
        end
        
        if (sameedge(pc)~=1 && sharebound(nextpart,whichbound(nextpart,pc),pc))
            % particles were not propagated as a pair
            % other particle shares the boundary reached
            if (opt.verbose)
                disp(sprintf('Iter %d Shared boundary: %d %d ',pc,nextpart,otherpart))
            end
            
            % propagate other particle in space           
            % numerical issues if one particle takes infinitessimally tiny step
            % the other should then not be propagated
            if (dt>0) 
                if (nodepos(otherpart,pc)==0) % particle currently on edge
                    edgeinfo = edgeinfoboth(otherpart,pc);
                    if (isempty(edgeinfo.lens))
                        error('missing edge info')
                    end
                    shortlim = (dt<edgeinfo.tstar);
                    [whichedge,xsamp,success] = sampleHopPosition(edgeinfo.lens,edgeinfo.uroots,...
                        edgeinfo.rpu2,dt,shortlim);
                    if (whichedge==1)% where is particle on actual edge
                        edgepos(otherpart,pc) = edgepos(otherpart,pc) - edgeinfo.lens(1)+xsamp;
                    else
                        edgepos(otherpart,pc) = edgepos(otherpart,pc) + edgeinfo.lens(2)-xsamp;
                    end
                else % other particle currently on node
                    nc = nodepos(otherpart,pc);
                    deg = NT.degrees(nc);
                    if (sameedge(pc)==2) % hopping particle is on an edge adjacent to this node
                        edgeinfo = edgeinfoboth(otherpart,pc);
                        shortlim = (dt<edgeinfo.tstar);
                        [whichedge,xsamp,success] = sampleHopPosition(edgeinfo.lens,edgeinfo.uroots,...
                            edgeinfo.rpu2,dt,shortlim);
                        % new distance from node
                        distnode = edgeinfo.lens(whichedge)-xsamp;
                        newedge = randi(deg,1); % pick new edge randomly
                        ec = NT.nodeedges(nc,newedge);
                        if (NT.edgenodes(ec,1)==nc) % outgoing edge
                            edgepos(otherpart,pc) = distnode;
                        elseif (NT.edgenodes(ec,2)==nc) % incoming edge
                            edgepos(otherpart,pc) = NT.edgelens(ec)-distnode;
                        else
                            error('bad connectivity in partial node hop')
                        end
                    else % full node neighborhood available
                        shortlim = (dt<nethopinfo.tstar(nc));
                        nr = nethopinfo.nroots(nc);
                        deg = NT.degrees(nc);
                        lens = NT.edgelens(NT.nodeedges(nc,1:deg))';
                        [whichedge,xsamp,success] = sampleHopPosition(lens,nethopinfo.uroots(1:nr,nc),...
                            nethopinfo.rpu2(1:deg,1:nr,nc),dt,shortlim);
                        if (~success); error('failed to sample hop position from node'); end
                        
                        nodepos(otherpart,pc) = 0;
                        ec = NT.nodeedges(nc,whichedge);
                        edgeid(otherpart,pc) = ec;
                        % track how far along the edge particle is located
                        if (NT.edgenodes(ec,1)==nc) % outgoing edge
                            edgepos(otherpart,pc) = lens(whichedge)-xsamp;
                        elseif (NT.edgenodes(ec,2)==nc)
                            edgepos(otherpart,pc) = xsamp;
                        else
                            error('network connectivity problem. Hopped to nonadjacent edge')
                        end
                    end
                end
            end
            
            % sample new event times for both particles
            psample = [1,2];
            [tsamp, whichbound(:,pc), nextnodepos(:,pc),nextedgeid(:,pc),nextedgepos(:,pc),...
                edgeinfotmp,rectinfotmp,sameedge(pc)] = ...
                sampleNextEvent(NT,psample,nodepos(:,pc),edgeid(:,pc),edgepos(:,pc),nethopinfo);
            
            nexteventtime(:,pc) = curtime(pc)+tsamp;
            
            rectinfo(pc) = rectinfotmp;
            edgeinfoboth(:,pc) = edgeinfotmp;
            
            % update shared boundaries
            sharebound(:,:,pc) = updateShareBound(NT,nodepos(:,pc),edgepos(:,pc),edgeid(:,pc));
            
        elseif (sameedge(pc)==1)
            % particles on same edge
            if (opt.verbose)
                disp(sprintf('Iter %d Same edge',pc))
            end
            % sample new event times for both particles
            psample = [1,2];
            [tsamp, whichbound(:,pc), nextnodepos(:,pc),nextedgeid(:,pc),nextedgepos(:,pc),...
                edgeinfotmp,rectinfotmp,sameedge(pc)] = ...
                sampleNextEvent(NT,psample,nodepos(:,pc),edgeid(:,pc),edgepos(:,pc),nethopinfo);
            nexteventtime(:,pc) = curtime(pc)+tsamp;
            
            rectinfo(pc) = rectinfotmp;
            edgeinfoboth(:,pc) = edgeinfotmp;
            
        else
            if (opt.verbose)
                disp(sprintf('Iter %d One particle jump: %d %d',pc,nextpart,otherpart))
            end
            
            % particle jumped to a boundary that is not shared
            % leave the other particle alone, and only sample a new
            % time for the jumping particle
            psample = nextpart;
            [tsamp, whichbound(nextpart,pc), nextnodepos(nextpart,pc),...
                nextedgeid(nextpart,pc),nextedgepos(nextpart,pc),...
                edgeinfotmp,rectinfotmp,sameedge(pc)] = ...
                sampleNextEvent(NT,psample,nodepos(:,pc),edgeid(:,pc),edgepos(:,pc),nethopinfo);
            edgeinfoboth(nextpart,pc) = edgeinfotmp(1);
            
            nexteventtime(nextpart,pc) = curtime(pc)+tsamp;
        end
        
        % update boundary sharing information
        sharebound(:,:,pc) = updateShareBound(NT,nodepos(:,pc),edgepos(:,pc),edgeid(:,pc));
    end
end


end