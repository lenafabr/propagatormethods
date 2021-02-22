function propagateSpatialPair(dt,psample,nodepos,edgeid,edgepos,edgehopinfo,edgehopinfo2)
% propagate particle pair spatially within its domain (no-passage propagator)

if (~exist('xi','var'))
    xi = -log(nethopinfo.epsilon);
end

    % sample next event time for particle pc
    %%
    newedgeid = edgeid;
    newnodepos = nodepos;
       
    if (length(psample)==2 && all(nodepos==0) && edgeid(1)==edgeid(2))
        %% both particles on same edge
        flip = edgepos(2)<edgepos(1);
        
        shortlim = (dt<edgehopinfo.tstar);
        [whichedge,xsamp,success] = sampleHopPosition(edgehopinfo.lens,edgehopinfo.uroots,...
                    edgehopinfo.rpu2,dt,shortlim);
                
        shortlim = (dt<edgehopinfo2.tstar);      
        [whichedge2,xsamp2,success2] = sampleHopPosition(edgehopinfo2.lens,edgehopinfo2.uroots,...
                    edgehopinfo2.rpu2,dt,shortlim);
        
        if (flip)
            
                
                
        x10 = min(edgepos);
        x20 = max(edgepos);
        L= NT.edgelens(edgeid(1));
        
        [tsamp,newx1,newx2,whichbound,tiltrect] = sampleHopPair(x10,x20,L,xi);                
        
        if (flip)
            newedgepos = [newx2 newx1];
            pc1 = 2; pc2 = 1;
        else
            newedgepos = [newx1 newx2]
            pc1 = 1; pc2 = 2;
        end        
        
        if (whichbound==1) % reached a node
            newedgeid(pc1) = 0;
            newnodepos(pc1) = NT.edgenodes(edgeid(1),1)
        elseif (whichbound==2)
            newedgeid(pc2) = 0;
            newnodepos(pc2) = NT.edgenodes(edgeid(1),2)            
        end
        
    else
        error('particles not paired up')
    end     
    
end