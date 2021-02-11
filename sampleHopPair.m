function [tsamp,newx1,newx2,whichbound,tiltrect,rectinfo] = sampleHopPair(x10,x20,L,xi)
%% sample a hop time and position for a pair of particles
% the hop should propagate the particles to the edge of 
% a tilted rectangle that fits within their triangle of propagation
% assumes *both* particles have equal diffusivity (D=1)
% INPUT:
% x10, x20 = initial particle positions along the interval (btwn 0 and L)
% L = interval length
% xi = cutoff for using short-time limit
% OUTPUT:
% tsamp = sampled time to leave tilted rectangular region
% newx1, newx2 = new particle positions after this time
% whichbound: how the pair exits their domain
% -1 = particles hit each other
% 0 = no boundary hit
% 1 = first particle hit lower boundary
% 2 = last particle hit upper boundary


if (~exist('xi')) % set default
    xi = - log(1e-14);
end

if (x10<0 || x10 > L || x20 <0 || x20>L)
    error(sprintf('bad x values in sampleHopPair: %f %f %f', x10, x20, L))
end

%x1 = min([x10,x20]);
%x2 = max([x10,x20]);
s2 = sqrt(2);

% tilt coordinate system
u0 = 1/s2*abs(x20-x10);
v0 = 1/s2*(x20+x10);

% set up bounding rectangle if tilted
if (v0 < L/s2)
    H = (u0+v0)/2; 
    W = s2*L-2*H;
    
    % distances to rectangular boundaries
    disttilt = [v0-H,W-v0,u0,H-u0];
else
    H = (s2*L-v0+u0)/2;
    W = s2*L-2*H;    
    disttilt = [(s2*L-v0-H),v0-H,u0,H-u0];
end


if (x10<x20);
    flip = 0;
    x1= x10; x2 = x20;
else
    flip = 1;
    x1 = x20; x2 = x10;
end
   
% set up bounding rectangle if not tilted
b = (x2-x1)/2;
Wstraight = x1+b;
Hstraight = L-Wstraight;

diststraight = [x1,Wstraight-x1,L-x2,Hstraight-(L-x2)];
    
% set your bounding rectangle to have the largest min distance 
% from current position to a boundary
if (min(diststraight)>min(disttilt))
    % sample time to hop out of non-tilted rectangle (particles can reach
    % domain edges)
    tiltrect = 0;    
        
    [tsamp,newpos,whichhit,edgeinfo1,edgeinfo2] = ...
        sampleRectangleHop([x1,Hstraight-(L-x2)],Wstraight,Hstraight,xi);
     
    if (flip)
        newx2 = newpos(1); newx1 = L - (Hstraight-newpos(2));      
    else
        newx1 = newpos(1); newx2 = L - (Hstraight-newpos(2));          
    end           
    
   switch whichhit
       case 1
           whichbound = 1; % lower boundary hit by first particle
       case 4
           whichbound = 2; % upper boundary hit by last particle
       otherwise
           whichbound = 0; % no boundary hit
   end
else
    tiltrect = 1;
    % sample time to hop out of tilted rectangle (particles can react
    [tsamp,newpos,whichhit,edgeinfo1,edgeinfo2] = sampleRectangleHop([v0-H,u0],W,H,xi);
    
    u = newpos(2);
    v = newpos(1)+H;
    
    % new positions of the two particles
    newx1 = (v-u)/s2;
    newx2 = (v+u)/s2;
    
    if (whichhit==3)
        whichbound = -1; % particles react
    else
        % particles do not hit a triangle boundary
        whichbound = 0;
    end
        
end

% rectangle information
rectinfo = struct('H',H,'W',W,'tiltrect',tiltrect,'edgeinfo1',edgeinfo1,'edgeinfo2',edgeinfo2,'L',L);
if (~tiltrect)
    rectinfo.H = Hstraight; rectinfo.W = Wstraight;
end
