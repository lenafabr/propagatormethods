function newedgepos = propagateSpatialPair(dt,rectinfo)
% propagate particle pair spatially within its domain (no-passage propagator)
% newedgepos is in order of initial particle location

%% both particles on same edge

edgeinfo1 = rectinfo.edgeinfo1;
edgeinfo2 = rectinfo.edgeinfo2;

shortlim = (dt<edgeinfo1.tstar);
[whichedge,xsamp,success] = sampleHopPosition(edgeinfo1.lens,edgeinfo1.uroots,...
    edgeinfo1.rpu2,dt,shortlim);

if (whichedge==2) % flip coord system if ended up above x0
    xsamp = rectinfo.W-xsamp;
end

shortlim = (dt<edgeinfo2.tstar);
[whichedge2,xsamp2,success2] = sampleHopPosition(edgeinfo2.lens,edgeinfo2.uroots,...
    edgeinfo2.rpu2,dt,shortlim);

if (whichedge2==2) % flip coord system if ended up above x0
    xsamp2 = rectinfo.H-xsamp2;
end

if (rectinfo.tiltrect) % tilted rectangle
    u = xsamp2;
    v = xsamp+rectinfo.H;
    
    % new positions of the two particles
    newedgepos(1) = (v-u)/sqrt(2);
    newedgepos(2) = (v+u)/sqrt(2);
else
    newedgepos(1) = xsamp;
    newedgepos(2) = rectinfo.L - (rectinfo.H - xsamp2);
end


end