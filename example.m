% Example code for obtaining mean first passage times on a network using 
% both analytic methods and simulations

addpath('networktools')

%% Load example network and plot
NT=NetworkObj('example.net');
NT.setCumEdgeLen();

% plot network, with node indices labeled
NT.plotNetwork(struct('labels',1))
title('Example network with node indices labeled')

%% Analytic example 1
% Choose two target nodes and find mean exit time (and variance) for 
% particles starting on all nodes of the network

targets = [4 17];                   % set target node indices
absrates = zeros(NT.nedge,1);       % no absorbing edges

[MFPTs,Vars,~] = networkMFPTanalytic(NT,absrates,targets);
MFPTanaly = mean(MFPTs)
STDanaly = mean(sqrt(Vars))

figure
colormap spring
opt=struct();                       % struct for plotting options
opt.nodesize=50;
opt.nodecolor=MFPTs;
NT.plotNetwork(opt)
hold on
scatter(NT.nodepos(targets,1),NT.nodepos(targets,2),75,[0 .75 0],'o','Filled')
hold off
colorbar
title('MFPT to targets (green)')

%% Analytic example 2
% Choose 3 absorbing edges with finite reaction rate, and find MFPTs/Vars
% for particles starting uniformly on all edges of the network

targets=[];                         % no target nodes
absedges = [4 10 22];               % choose which edges are absorbing
absrates = zeros(NT.nedge,1);
absrates(absedges) = [30 10 20];    % each edge can have a different absorption rate

[MFPTs,Vars,~] = networkMFPTanalyticEdges(NT,absrates,targets);
mean(MFPTs)
mean(sqrt(Vars))

% plot network and color by MFPT along each edge
figure
cmap=spring(100);colormap(cmap);
opt=struct();
opt.plotnodes=[];
opt.edgeplotopt={'LineWidth',2.5};
opt.edgecolor=cmap(round(100*MFPTs/max(MFPTs)),:);
opt.edgecolor(absedges',:)=ones(length(absedges),3).*[0 .75 0];
NT.plotNetwork(opt)
caxis([min(MFPTs) max(MFPTs)])
colorbar
title('Color = MFPT starting from each edge')

%% Simulation example 1 (no save times)

% calculate propagator roots
nethopinfo = networkPropagatorRoots(NT,struct('epsilon',1e-14));
% tabulate hopping times for each node and neighbor
nethopinfo = tabulateHopTimes(nethopinfo,NT);

options = struct();
options.targetnodes = [4 17];       % set target nodes
options.startedgeuniform = 0;       % option for starting on edges
npart = 500;                        % number of particles to simulate

[targethittime,~] = simulateNetworkHopper(NT,nethopinfo,npart,options);

% the MFPT and standard deviation should match analytic example 1
MFPTsim = mean(targethittime(~isinf(targethittime)));
STDsim   = std(targethittime(~isinf(targethittime)));

disp(sprintf('Simulated and analytic MFPTs: %f %f', MFPTsim, MFPTanaly))
disp(sprintf('Simulated and analytic st.dev.: %f %f',STDsim,STDanaly))
%% Simulation example 2 (with save times)

options = struct();
options.savetimes = 0.01:0.01:1;    % specify save times
options.targetnodes = [];           % no target nodes
options.maxtime = 1;                % max time to run simulation
npart = 5;

[~,~,savepos,~] = simulateNetworkHopper(NT,nethopinfo,npart,options);

% convert to trajectories
tracklist = savepos2tracklist(NT,savepos);

% plot single trajectory over network
figure
NT.plotNetwork()
hold all
colormap spring
track = tracklist{2}; % use trajectory of particle 2
scatter(track(:,1),track(:,2),25,options.savetimes,'filled')
hold off
colorbar
title('Simulated trajectory of a single particle, color = time')

%% Particle pair simulation
% Reuse nethopinfo from above and find mean time to react for two particles
% diffusing on the network.

options = struct();
options.startedgeuniform = 1;       % must start on edges for pair sims, chooses random location along random edge for each particle
npart = 100;                        % number of particle pairs to simulate

[reacttimes] = simulateNetworkHopper_pair(NT,nethopinfo,npart,options);
encountertime = mean(reacttimes)


