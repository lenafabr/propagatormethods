% Example code for obtaining mean first passage times on a network using 
% both analytic methods and simulations

addpath('networktools')

%% Load example network and plot
NT=NetworkObj('example.net');
NT.plotNetwork

%% Analytic example 1
% Choose two target nodes and find mean exit time (and variance) for 
% particles starting on all nodes of the network

targets = [4 17];                   % set targets
absrates = zeros(NT.nedge,1);       % no absorbing edges

[MFPTs,Vars,~] = networkMFPTanalytic(NT,absrates,targets);
mean(MFPTs)
mean(sqrt(Vars))

figure
colormap spring
NT.plotNetwork
hold on
scatter(NT.nodepos(:,1),NT.nodepos(:,2),60,MFPTs,'o','Filled')
scatter(NT.nodepos(targets,1),NT.nodepos(targets,2),75,[0 .75 0],'o','Filled')
colorbar
hold off

%% Analytic example 2
% Choose 3 absorbing edges with finite reaction rate, and find MFPTs/Vars
% for particles starting uniformly on all edges of the network

targets=[];                         % no target nodes
absedges = [4 10 22];               % choose which edges are absorbing
absrates = zeros(NT.nedge,1);
absrates(absedges) = [30 10 20];    % each edge can have a different rate

[MFPTs,Vars,~] = networkMFPTanalyticEdges(NT,absrates,targets);
mean(MFPTs)
mean(sqrt(Vars))

cmap=colormap(spring);
figure
NT.plotNetwork
hold on
for i=1:NT.nedge
    plot(NT.nodepos(NT.edgenodes(i,:),1), NT.nodepos(NT.edgenodes(i,:),2),'Color',cmap(round(256*MFPTs(i)/max(MFPTs)),:),'LineWidth',2.5)
end
for i=absedges
    plot(NT.nodepos(NT.edgenodes(i,:),1), NT.nodepos(NT.edgenodes(i,:),2),'Color',[0 .75 0] ,'LineWidth',2.5)
end
hold off

%% Simulation example 1 (no save times, no edge propagation)

% calculate propagator roots
nethopinfo = networkPropagatorRoots(NT,struct('epsilon',1e-14));
% tabulate hopping times for each node and neighbor
nethopinfo = tabulateHopTimes(nethopinfo,NT);

options = struct();
options.targetnodes = [4 17];       % set target nodes
options.startedgeuniform = 0;       % option for starting on edges
npart = 500;                        % number of particles to simulate
[targethittime,~] = simulateNetworkHopper(NT,nethopinfo,npart,options);

% The MFPT and standard deviation should match analytic results
MFPTsSim = mean(targethittime(~isinf(targethittime)))
STDSim   = std(targethittime(~isinf(targethittime)))

%% Particle pair simulation
% Reuse nethopinfo from above

options = struct();
options.startedgeuniform = 1;       % must start on edges for pair sims
npart = 100;                        % number of particle pairs to simulate

[reacttimes] = simulateNetworkHopper_pair(NT,nethopinfo,npart,options);
mean(reacttimes)


