function [targethittime,hittimes,savepos,opt,curtime,nodepos,step] = simulateNetworkHopper(NT,nethopinfo,npart,options)
% simulate particles hopping on a network
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
opt.D = 1;
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

curtime = zeros(1,npart);
done = false(1,npart);

% note first hitting times to all nodes
hittimes = inf*ones(npart,NT.nnode);
targethittime = inf*ones(npart,1);

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
    edgeid = datasample(1:NT.nedge,npart,'Weights',lens);
    % sample starting edge position uniformly
    edgepos = rand(npart,1).*lens(edgeid);    
end

nextsave = ones(npart,1);
if (isempty(opt.savetimes))
    nextsavetime = inf*ones(npart,1);
else
    nextsavetime = opt.savetimes(1)*ones(npart,1);
end

savepos = zeros(npart,2,length(opt.savetimes));

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
%%
for step = 1:opt.maxsteps
    if (opt.stopfirsthit & any(done))
        % if any of the particles have hit, set new maxtime and run all particles until then
        opt.maxtime = min(targethittime);
    end
    %[step curtime nextsavetime]
    
    %if (mod(step,opt.printevery)==0); disp([step max(curtime) nnz(~isinf(targethittime))]); end
    if (mod(step,opt.printevery)==0); disp([step min(curtime) nnz(done)]); end
    
    for pc = 1:npart
        if (done(pc)); continue; end
        
        if (on_edge(pc)) % particle is currently on an edge, and needs to hop to a node)
            
            % sample time to hop to node
            ec = edgeid(pc); edgelen = NT.edgelens(ec);
            [whichleave,tsamp, success,edgehopinfo] = sampleHopTime_edge(edgepos(pc),edgelen,xi);
            if (~success)
                error('failed to sample time to leave edge')
            end
            % convert to real time
            tsamp = tsamp/opt.D;
            
            dosave = curtime(pc)+tsamp > nextsavetime(pc);
            if (dosave)                              
                % propagate to the next save point only                
                dt = nextsavetime(pc)-curtime(pc);
                
                % figure out where along the edge the particle is at the
                % save-time
                shortlim = (dt<edgehopinfo.tstar);
                [whichedge,xsamp,success] = sampleHopPosition(edgehopinfo.lens,edgehopinfo.uroots,...
                    edgehopinfo.rpu2,dt*opt.D,shortlim);
                if (~success); error('failed to sample position along edge'); end
                
                if (whichedge==1) % flip coord system if ended up above x0
                    edgepos(pc) = xsamp;
                else
                    edgepos(pc) = edgelen-xsamp;
                end       
            else
                % particle got to node before next save time
                on_edge(pc) = false;
                curtime(pc) = curtime(pc)+tsamp;
                nodepos(pc) = NT.edgenodes(ec,whichleave);
            end
            
        else % particle is currently on a node and needs to hop to another node
            
            % keep track of first passage to all nodes
            if (isinf(hittimes(pc,nodepos(pc))))
                hittimes(pc,nodepos(pc)) = curtime(pc);
                
                if (opt.stopwhenallhit && all(~isinf(hittimes(pc,:))))
                    % stop evolving this particle once it has hit all nodes 
                    done(pc) = true;
                end
            end
            
            % keep track of first passage to one of a set of target nodes
            % stop when all targets hit
            
            targettimes = hittimes(pc,opt.targetnodes);
            
            if (opt.hitfirst & any(~isinf(targettimes)))
                % stop particle after it hits the first of the target nodes
                done(pc) = true;
                
                targethittime(pc) = min(targettimes);
            end
            
            if (~opt.hitfirst & all(~isinf(targettimes)))
                % stop particle after it hits all the targets
                done(pc) = true;
            end
            
            
            nc = nodepos(pc); deg = NT.degrees(nc);
            lens = NT.edgelens(NT.nodeedges(nc,1:deg))';
                       
            % which node to move to next            
            P = Pvals(1:deg,nc);
            
            if (opt.usetabulation)                             
                if (opt.npresample> 0)
                    % use presampled hopping directions and times
                   
                    % select where to hop
                    ind = presampleWhichLeaveInd(nc);      
                    whichleave = presampleWhichLeave(ind,nc);
                    presampleWhichLeaveInd(nc) = ind+1;
                    if (ind==opt.npresample)
                        % ran out of presampled values, sample again
                        presampleWhichLeave(:,nc) = datasample(1:deg,opt.npresample,'Weights',P);
                        presampleWhichLeaveInd(nc) = 1;
                    end
                    
                    % select hopping time
                    ind = presampleInd(whichleave,nc);                    
                    tsamp = presampleTimes(ind,whichleave,nc);
                    presampleInd(whichleave,nc) = ind+1;
                    if (ind==opt.npresample)
                        if (opt.verbose)
                            disp(sprintf('Resampling: %d %d', whichleave, nc))
                        end
                        % ran out of presampled times for this node and neighbor
                        % sample some more
                        nt = nethopinfo.nhoptime(whichleave,nc);
                        tlist = nethopinfo.hoptimelist(1:nt,whichleave,nc);
                        cumdist = nethopinfo.hoptimecumdist(1:nt,whichleave,nc);
                        
                        w = rand(opt.npresample,1);
                        presampleTimes(:,whichleave,nc) = interp1(cumdist,tlist,w);                       
                        presampleInd(whichleave,nc) = 1;
                    end                
                else
                    % use tabulated distribution of hopping times, sampling
                    % directly for this one particle
                    nt = nethopinfo.nhoptime(:,nc);
                    tlist = nethopinfo.hoptimelist(:,:,nc);
                    cumdist = nethopinfo.hoptimecumdist(:,:,nc);
                    [whichleave,tsamp] = sampleHopTime_tab(P,tlist,cumdist,nt);
                end
            else % direct function evaluation
            nr = nroots(nc);
            uroots2 = uroots(1:nr,nc).^2;                        
            [whichleave,tsamp, success] = sampleHopTime(P,uroots2,rpu2(1:deg,1:nr,nc),lens,tstar(nc));                        
            if (~success)
                error('failed to sample time')
            end
            end
            
            %tsamp
            
            nextnode = NT.nodenodes(nc,whichleave);
            % scale to true diffusivity
            tsamp = tsamp/opt.D;
            
            dosave = curtime(pc)+tsamp > nextsavetime(pc);
            
            if (dosave)
                % propagate to the next save point only
                
                dt = nextsavetime(pc)-curtime(pc);                
                
                % figure out where the particle is at the save-point
                shortlim = (dt<tstar(nc));
                nr = nethopinfo.nroots(nc);
                [whichedge,xsamp,success] = sampleHopPosition(lens,uroots(1:nr,nc),rpu2(1:deg,1:nr,nc),dt,shortlim);
                if (~success); error('failed to sample hop position'); end
                               
                on_edge(pc) = true;
                % track which edge it is on, and how far along edge
                ec= NT.nodeedges(nc,whichedge);
                edgeid(pc) = ec;
                
                % track how far along the edge particle is located
                if (NT.edgenodes(ec,1)==nc) % outgoing edge
                    edgepos(pc) = lens(whichedge)-xsamp;
                elseif (NT.edgenodes(ec,2)==nc)
                    edgepos(pc) = xsamp;
                else
                    error('network connectivity problem. Hopped to nonadjacent edge')
                end
            else
                curtime(pc) = curtime(pc)+tsamp;
                nodepos(pc) = nextnode;
            end
            
        end
        
%         if (edgeid(pc)<=0)
%             error('bad edge id')
%         end
        
        
        if (dosave) % particle reached a savepoint
            curtime(pc) = nextsavetime(pc);
            % save positions in terms of edge id and position along edge
            savepos(pc,:,nextsave(pc)) = [edgeid(pc) edgepos(pc)];
            
%             if (edgeid(pc)<=0)
%                 error('bad edge id')
%             end
            
            nextsave(pc) = nextsave(pc)+1;
            if (nextsave(pc)>length(opt.savetimes))
                nextsavetime(pc) = inf;
            else
                nextsavetime(pc) = opt.savetimes(nextsave(pc));
            end
                        
        end
        
        if (curtime(pc) > opt.maxtime)
            done(pc) = true;
        end
    end
    
    if (all(done)); break; end
end
