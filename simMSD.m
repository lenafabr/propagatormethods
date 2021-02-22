%% load in Zuben's Voronoi network
load('../results/SimsVsAnaly20210119.mat','NT')
NT.interpolateEdgePaths(3);
NT.setEdgeLens();
NT.setCumEdgeLen(1:NT.nedge,1);
NT0 = copy(NT);

%%
NT = copy(NT0);
% remove edges
nrm = 90;
NT=decimateEdges(NT,nrm);
NT.setCumEdgeLen(1:NT.nedge,1);
NT.plotNetwork()

NTdec = copy(NT);
%
npart = 1000; % 
nethopinfodec = networkPropagatorRoots(NT,struct('epsilon',1e-14));
nethopinfodec = tabulateHopTimes(nethopinfo,NT);

options = struct('startedgeuniform', true,'printevery',100);
options.savetimes = [1e-4 1e-4+logspace(-3,-1)];
options.maxtime = max(options.savetimes)+0.01;
options.targetnodes = [];

s = rng();
[~,~,saveposdec,opt,curtime,nodepos,step] = simulateNetworkHopper(NTdec,nethopinfo,npart,options);
%%
nethopinfo0 = networkPropagatorRoots(NT0,struct('epsilon',1e-14));
nethopinfo0 = tabulateHopTimes(nethopinfo0,NT0);
s2 = rng();
[~,~,savepos0,opt,curtime,nodepos,step] = simulateNetworkHopper(NT0,nethopinfo0,npart,options);


%% convert to real-space positions
tracklist = savepos2tracklist(NTdec,saveposdec);
%
realpos = zeros(size(saveposdec));
for tc = 1:length(tracklist)
    realpos(tc,:,:) = tracklist{tc}';
end

realposdec = realpos;

%
tracklist = savepos2tracklist(NT0,savepos0);
%
realpos = zeros(size(savepos0));
for tc = 1:length(tracklist)
    realpos(tc,:,:) = tracklist{tc}';
end
realpos0 = realpos;

save('~/writeup/propagatormethods/workfig/MSDsims.mat')
%% get MSD
realpos = realposdec;
diffs = realpos - realpos(:,:,1);
sqdists = sum(diffs.^2,2);
MSDdec = squeeze(mean(sqdists,1));
stddec = squeeze(std(sqdists,1))/sqrt(npart);
tdiff = options.savetimes - options.savetimes(1);

realpos = realpos0;
diffs = realpos - realpos(:,:,1);
sqdists = sum(diffs.^2,2);
MSD0 = squeeze(mean(sqdists,1));
std0 = squeeze(std(sqdists,1))/sqrt(npart);
tdiff = options.savetimes - options.savetimes(1);


errorbar(tdiff,MSD0,std0,'b.-','LineWidth',2,'MarkerSize',10)
hold all
errorbar(tdiff,MSDdec,stddec,'m.-','LineWidth',2,'MarkerSize',10)
loglog(tdiff,2*tdiff,'k--','LineWidth',2)
hold off
set(gca,'XScale','log','YScale','log',...
    'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',14)
xlabel('dimensionless time $tD/L^2$')
ylabel('dimensionless MSD $\left<x^2\right>/L^2$')

%% plot networks
figure
plotopt = struct('plotedgepath',0);
plotopt.edgeplotopt = {'LineWidth',3, 'Color','k'};
plotopt.nodeplotopt = {'k','filled'};
plotopt.nodesize=10;
NT0.plotNetwork(plotopt);
set(gca,'Visible','off')
% ------------------------
%% Playing around with other network structures
% -------------------------

%% try on decimated hex network
%% make a hex network
NT = makeHexNetwork(36);
NT.nodepos = NT.nodepos*60+COM-[20,20];
%% constrain within circle
COM = mean(NT.nodepos,1);
NT.nodepos = NT.nodepos-COM;
nodedist = sqrt(sum((NT.nodepos).^2,2));
keepind = find(nodedist<10);
NT.keepNodes(keepind) ;  

NT.setupNetwork();
NT.interpolateEdgePaths(2);
NT.setCumEdgeLen(1:NT.nedge,1);


%% decimate edges
[NTdec,whichremoved]=decimateEdges(NT,0.5*NT.nedge);
%%
NTdec.setupNetwork();
NTdec.interpolateEdgePaths(2);
NTdec.setCumEdgeLen(1:NTdec.nedge,1);
%%
NTdec.plotNetwork()

% ----
%% Load ER network

dirname = '~/proj/propagatormethods/networks/';
fglob = 'WT_*.mat';
files = dir([dirname fglob]);
fnames = {files.name};

rc = 2;
load([dirname fnames{rc}],'NT')

NTER = copy(NT);

NTER.plotNetwork()
%% get a convex hull
ch = convhull(NTER.nodepos);

COM = mean(NTER.nodepos(ch,:)); % rough center

%% make a hex network
NT = makeHexNetwork(36);
NT.nodepos = NT.nodepos*60+COM-[30,30]

%% keep only nodes within convex hull
inch = inpolygon(NT.nodepos(:,1),NT.nodepos(:,2),NTER.nodepos(ch,1),NTER.nodepos(ch,2));
ind = find(inch);
NT.keepNodes(ind);
%%
NT.setupNetwork();
NT.interpolateEdgePaths(2);
NT.setCumEdgeLen(1:NT.nedge,1);

%%
NTER.plotNetwork()
hold all
NT.plotNetwork()
hold off