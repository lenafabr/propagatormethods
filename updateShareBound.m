function sharebound = updateShareBound(NT,nodepos,edgepos,edgeid)
% for 2 particles, decide whether they share boundaries
% sharebound(i,j) = true if particle i shares its j-th boundary with
% another particle

sharebound = false(2,max(NT.maxdeg,2));

if nodepos(1) == 0
    % particle 1 on an edge
    ec = edgeid(1);
    
    if (nodepos(2)==0)
        %particle 2 also on edge
        if (edgeid(1)==edgeid(2))
            if (edgepos(1)<edgepos(2))
                sharebound(1,2) = true;
                sharebound(2,1) = true;
            else
                sharebound(1,1) = true;
                sharebound(2,2) = true;
            end
        else
            ec2 = edgeid(2);
            for bc1 = 1:2
                for bc2 = 1:2
                    if NT.edgenodes(ec,bc1) == NT.edgenodes(ec2,bc2)
                        sharebound(1,bc1) = true;
                        sharebound(2,bc2) = true;
                    end
                end
            end
        end
    else  
        % particle 2 on node
        nc = nodepos(2);
        %ec = edgeid(1);
                        
        for ecc = 1:NT.degrees(nc)
            if NT.nodeedges(nc,ecc) == ec %particle is on an edge adjacent to this node
                sharebound(2,ecc) = true;
                if (NT.edgenodes(ec,1)==nc)
                    sharebound(1,1)=true;
                elseif (NT.edgenodes(ec,2)==nc)
                    sharebound(1,2) = true;
                end
            else % particle may be on an edge that shares a node boundary with this node
                
                if NT.nodenodes(nc,ecc) ==NT.edgenodes(ec,1)
                    sharebound(2,ecc) = true;
                    sharebound(1,1) = true;
                elseif NT.nodenodes(nc,ecc) ==NT.edgenodes(ec,2)
                    sharebound(2,ecc) = true;
                    sharebound(1,2) = true;
                end
            end
        end
    end              
else
    % particle 1 is on node
    nc = nodepos(1);
    deg = NT.degrees(nc);
    
    if (nodepos(2) == 0) % 2 is on edge    
        ec= edgeid(2);
        
        for ecc = 1:deg
            if NT.nodeedges(nc,ecc) == ec %particle is on an edge adjacent to this node
                sharebound(1,ecc) = true;
                if (NT.edgenodes(ec,1)==nc)
                    sharebound(2,1)=true;
                elseif (NT.edgenodes(ec,2)==nc)
                    sharebound(2,2) = true;
                end
            else % particle may be on an edge that shares a node boundary with this node                                
                if NT.nodenodes(nc,ecc) ==NT.edgenodes(ec,1)
                    sharebound(1,ecc) = true;
                    sharebound(2,1) = true;
                elseif NT.nodenodes(nc,ecc) ==NT.edgenodes(ec,2)
                    sharebound(1,ecc) = true;
                    sharebound(2,2) = true;
                end
            end
        end
    else % 2 is on node
        nc2 = nodepos(2);  
        deg2 = NT.degrees(nc2);
        for cc = 1:NT.degrees(nc)
            if (NT.nodenodes(nc,cc)==nc2)
                % nodes are directly adjacent
                sharebound(1,cc) = true;
            elseif ismember(NT.nodenodes(nc,cc),NT.nodenodes(nc2,1:deg2))
                % nodes share an adjacent node
                sharebound(1,cc) = true;
            end
        end
        for cc = 1:NT.degrees(nc2)
            if (NT.nodenodes(nc2,cc)==nc)
                sharebound(2,cc)=true; 
            elseif ismember(NT.nodenodes(nc2,cc),NT.nodenodes(nc,1:deg))
                sharebound(2,cc) = true; 
            end
        end
    end
end

end