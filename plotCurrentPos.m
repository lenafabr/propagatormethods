function plotCurrentPos(NT,nodepos,edgeid,edgepos)
% plot current particle positions
npart = length(edgeid);
savePos = zeros(npart,2,1);
if (size(edgeid,1)<size(edgeid,2))    
    savePos(:,:,1) = [edgeid' edgepos];
else
    savePos(:,:,1) = [edgeid edgepos];
end

% pretend node positions are on edges
for pc = 1:npart
    if nodepos(pc)>0
        nc = nodepos(pc);
        ec = NT.nodeedges(nc,1);
        edgeid(pc) = ec;
        if (NT.edgenodes(ec,1)==nc)
            edgepos(pc) = 0;
        else
            edgepos(pc) = NT.edgelens(ec);
        end
    end
end

savePos = zeros(npart,2,1);
if (size(edgeid,1)<size(edgeid,2))    
    savePos(:,:,1) = [edgeid' edgepos];
else
    savePos(:,:,1) = [edgeid edgepos];
end
        
[tracklist] = savepos2tracklist(NT,savePos)

%%
NT.plotNetwork(struct('labels',1))
hold all
cmat = jet(2);
for tc = 1:length(tracklist)
    track = tracklist{tc};
    plot(track(:,1),track(:,2),'o','MarkerSize',8,'LineWidth',2,'Color',cmat(tc,:))
    text(track(:,1),track(:,2)+0.1*mean(NT.edgelens),sprintf('%d',tc),'Color',cmat(tc,:))
end
hold off

end