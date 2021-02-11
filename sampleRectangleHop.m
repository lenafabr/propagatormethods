function [tsamp,newpos,whichhit,edgehopinfo1,edgehopinfo2] = sampleRectangleHop(start,W,H,xi)
% propagate a particle in a rectangle
% picking the first time to reach one edge of the rectangle
% and then doing a no-passage propagation along the other dimension
% assumes isotropic diffusivity D=1
% start = (x,y) starting position in rectangle (between 0 and W, between 0
% and H)
% W, H = width, height of rectangle
% xi = -log(epsilon) used to flip from short t to long t limit

% newpos = new particle position
% whichhit = which boundary was hit; in order: low x, high x, low y, high y

% find absorbing time for each dimension


if (~exist('xi')) % set default
    xi = - log(1e-14);
end

[whichleave1,tsamp1,success1,edgehopinfo1] = sampleHopTime_edge(start(1),W,xi);
[whichleave2,tsamp2,success2,edgehopinfo2] = sampleHopTime_edge(start(2),H,xi);

xbound = [0,W];
ybound = [0,H];

if (tsamp1<tsamp2)
    % leave through the sides
    tsamp = tsamp1;
    newpos(1) = xbound(whichleave1);
    whichhit = whichleave1;
    
    shortlim = (tsamp<edgehopinfo2.tstar);
    
    %     if (tsamp == 0) % instantaneous transition, do not move vertically
    %         whichedge = 1;
    %         xsamp =
    % sample vertical position
    if (tsamp<1e-14) % essentially 0 time hop
        %whichedge = datasample(1:2,1);
        %xsamp = edgehopinfo2.lens(whichedge);
        newpos(2) = start(2);
    else
        [whichedge,xsamp,success] = sampleHopPosition(edgehopinfo2.lens,edgehopinfo2.uroots,...
            edgehopinfo2.rpu2,tsamp,shortlim);
        
        if (whichedge==1) % flip coord system if ended up above x0
            newpos(2) = xsamp;
        else
            newpos(2) = H-xsamp;
        end
    end
else
    % leave through the top or bottom
    tsamp = tsamp2;
    newpos(2) = ybound(whichleave2);
    whichhit = whichleave2+2;
    
    shortlim = (tsamp<edgehopinfo1.tstar);
    
    % sample horizontal position
    if (tsamp<1e-14) % essentially 0 time hop
        %whichedge = datasample(1:2,1);
        %xsamp = edgehopinfo1.lens(whichedge);
        newpos(1) = start(1);
    else
        [whichedge,xsamp,success] = sampleHopPosition(edgehopinfo1.lens,edgehopinfo1.uroots,...
            edgehopinfo1.rpu2,tsamp,shortlim);
        
        if (whichedge==1) % flip coord system if ended up above x0
            newpos(1) = xsamp;
        else
            newpos(1) = W-xsamp;
        end
    end
    
end

end