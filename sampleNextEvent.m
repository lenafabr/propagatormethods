function [tsamp, whichbound, newnodepos,newedgeid,newedgepos,saveedgehopinfo,rectinfo,sameedge]...
    = sampleNextEvent(NT,psample,nodepos,edgeid,edgepos,nethopinfo)

%if (~exist('xi','var'))
    xi = -log(nethopinfo.epsilon);
%end
    mindist = 1e-14;

if (any(edgepos<0))
    error('negative edgepos')
end

    % sample next event time for particle pc
    %%
    newedgepos = edgepos(psample);
    newedgeid = edgeid(psample);
    newnodepos = nodepos(psample);
      
     fields = {'uroots','rpu2','lens','tstar'};
        c = cell(length(fields),1);
        info1 = cell2struct(c,fields);
       saveedgehopinfo = [info1; info1];
    rectinfo = struct('H',[],'W',[],'tiltrect',[],'edgeinfo1',info1,'edgeinfo2',info1,'L',[]);
 
   % check if both particles are so close we should consider them reacted
   if (length(psample)==2)
       if (all(nodepos==0) && edgeid(1)==edgeid(2))
           % particles are on same edge
           diff = abs(edgepos(2)-edgepos(1));
           if (diff<mindist) % particles reacted
               tsamp=[0;0];
               whichbound=[-1;-1];
               newedgeid = edgeid;
               newnodepos = [0;0];
               newedgepos = mean(edgepos)*[1;1];               
               sameedge = 1;
               return
           end
       end
       for pc = 1:2
           pother = 3-pc;
           
           if (nodepos(pc)>0 & nodepos(pother)==0) % one particle on node, other on connected edge
               % other is on a connected edge
                if (NT.edgenodes(edgeid(pother),1)==nodepos(pc))
                    diff = edgepos(pother);
                elseif (NT.edgenodes(edgeid(pother),2)==nodepos(pc))
                    diff = NT.edgelens(edgeid(pother)) - edgepos(pother);
                else
                    diff = inf;
                end
                if (diff<mindist)
                        tsamp=[0;0];
                        whichbound=[-1;-1];
                        newedgeid = [0;0];
                        newnodepos = [nodepos(pother);nodepos(pother)];
                        newedgepos = [NaN;NaN];                 
                        sameedge = 2;
                        return               
                end
                
           end
       end
   end
       
    if (length(psample)==2 && all(nodepos==0) && edgeid(1)==edgeid(2))    
        sameedge = 1;               
      
        %% both particles on same edge
        flip = edgepos(2)<edgepos(1);
        
        x10 = min(edgepos);
        x20 = max(edgepos);
        L= NT.edgelens(edgeid(1));
        
        [tsamp,newx1,newx2,whichbound,tiltrect,rectinfo] = sampleHopPair(x10,x20,L,xi);                
        
        if (flip)
            newedgepos = [newx2 newx1];
            pc1 = 2; pc2 = 1;
        else
            newedgepos = [newx1 newx2];
            pc1 = 1; pc2 = 2;
        end        
        
        if (whichbound==1) % reached a node
            newedgeid(pc1) = 0;
            newnodepos(pc1) = NT.edgenodes(edgeid(1),1);
        elseif (whichbound==2)
            newedgeid(pc2) = 0;
            newnodepos(pc2) = NT.edgenodes(edgeid(1),2);            
        end
        
    else
        sameedge = 0;
        fields = {'H','W','tiltrect','edgeinfo1','edgeinfo2','L'};
        c = cell(length(fields),1);
        rectinfo = cell2struct(c,fields);
        
        for cc = 1:length(psample)
            pc = psample(cc); % sample for this particle
            pother = 3-pc; % the other particle
            
            if (nodepos(pc)==0) % particle is on an edge
                ec = edgeid(pc);
                if (nodepos(pother)>0 && NT.edgenodes(ec,1)==nodepos(pother))
                    % other particle is on adjacent node 1
                    sameedge = 2;
                    L = NT.edgelens(ec) - edgepos(pc)/2;
                    x0 = edgepos(pc)/2;                   
                    [whichbound(cc),tsamp(cc), success,edgehopinfo] = sampleHopTime_edge(x0,L,xi);
                    saveedgehopinfo(cc) = edgehopinfo;            
                    
                    if (whichbound(cc)==2) % hopped to node
                        newnodepos(cc) = NT.edgenodes(ec,2);
                        newedgepos(cc) = NaN;
                        newedgeid(cc) = 0;
                    else % hopped to edge
                        newnodepos(cc) = 0;
                        newedgepos(cc) = edgepos(pc)/2;
                    end
                elseif (nodepos(pother)>0 && NT.edgenodes(ec,2)==nodepos(pother))
                    sameedge = 2; % particles are on a node and adjacent edge
                    % other particle is on adjacent node 2
                    L = (NT.edgelens(ec) + edgepos(pc))/2;
                    x0 = edgepos(pc);                   
                    [whichbound(cc),tsamp(cc), success,edgehopinfo] = sampleHopTime_edge(x0,L,xi);
                    saveedgehopinfo(cc) = edgehopinfo;    
                    
                    if (whichbound(cc)==1) % hopped to node
                        newnodepos(cc) = NT.edgenodes(ec,1);
                        newedgepos(cc) = NaN;
                        newedgeid(cc) = 0;
                    else % hopped to edge
                        newnodepos(cc) = 0;
                        newedgepos(cc) = L;                        
                    end
                    
                else % hop on full edge
                    L = NT.edgelens(ec);
                    x0 = edgepos(pc);
                    [whichbound(cc),tsamp(cc), success,edgehopinfo] = sampleHopTime_edge(x0,L,xi);
                    saveedgehopinfo(cc) = edgehopinfo;    
                    
                    newedgepos(cc) = NaN;
                    newedgeid(cc) = 0;
                    newnodepos(cc) = NT.edgenodes(ec,whichbound(cc));
                end
            else % particle is on node
                nc = nodepos(pc);
                deg = NT.degrees(nc);
                
                boundtype = zeros(deg,1);
                lens = zeros(deg,1);
                for bc = 1:NT.degrees(nc)
                    ec = NT.nodeedges(nc,bc);
                    if (nodepos(pother)==0 && edgeid(pother) == ec) % particle on adjacent edge
                        if (NT.edgenodes(ec,1) == nc)
                            % keep track of minimal length from this node
                            lens(bc) = edgepos(pother)/2;
                        else
                            lens(bc) = (NT.edgelens(ec)-edgepos(pother))/2;
                        end
                        boundtype(bc) = 1; % edge boundary type
                    elseif (nodepos(pother)==NT.edgenodes(ec,1) || nodepos(pother)==NT.edgenodes(ec,2))
                        % particle on adjacent node
                        lens(bc) = NT.edgelens(ec)/2;
                        boundtype(bc) = 1;
                    else % no particle on this edge or adjacent node
                        lens(bc) = NT.edgelens(ec);    
                        boundtype(bc) = 0; % node boundary type
                    end
                end
                
                if (all(boundtype==0))
                    sameedge = 0;
                    % full sampling of node with different length edges
                    P = nethopinfo.Pvals(1:deg,nc);
                    nr = nethopinfo.nroots(nc);
                    uroots2 = nethopinfo.uroots(1:nr,nc).^2;                        
                    tstar = nethopinfo.tstar;
                    [whichbound(cc),tsamp(cc), success] = sampleHopTime(P,uroots2,nethopinfo.rpu2(1:deg,1:nr,nc),lens,tstar(nc));  
                    newnodepos(cc) = NT.nodenodes(nc,whichbound(cc));
                    newedgeid(cc) = 0;
                    newedgepos(cc) = NaN;
                else % cut off equal-length edges around this node
                    sameedge = 2; % particles are on node and adjacent edge
                    [minlen,bmin] = min(lens);
                   
                    L = 2*minlen;
                    x0 = minlen;
                    [whichbound(cc),tsamp(cc), success,edgehopinfo] = sampleHopTime_edge(x0,L,xi);
                    saveedgehopinfo(cc) = edgehopinfo;
                    
                    % pick what boundary you go to uniformly at random
                    whichbound(cc) = randi(NT.degrees(nc));
                    newec = NT.nodeedges(nc,whichbound(cc));
                    
                    if (boundtype(bmin)==0 && abs(NT.edgelens(newec) - minlen)<nethopinfo.epsilon)
                        % hopped to a node (minimum length actually goes to
                        % node)
                        newnodepos(cc) = NT.nodenodes(nc,whichbound(cc));
                        newedgeid(cc) = 0;
                        newedgepos(cc) = NaN;
                    else
                        % hopped to edge
                        newnodepos(cc)= 0;
                        newedgeid(cc) = newec;
                        if (NT.edgenodes(newedgeid(cc),1)==nc) % outgoing edge
                            newedgepos(cc) = minlen;
                        else %incoming edge
                            newedgepos(cc) = NT.edgelens(newedgeid(cc)) - minlen;
                        end
                    end
                                                                  
                end
            end
            
        end
    end       
end