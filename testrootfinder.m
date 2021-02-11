%% test root finding for Laplace inversion of Pik
lens = [0.339,1,2.7,3.2,10.2]; % arbitrary edge lengths
checkfun = @(u) sin(lens(2)*u).*(sum(cot(lens'*u)))

ulist = linspace(0,10,10000);

plot(ulist,checkfun(ulist),'.-')
ylim([-5,5])

%% get roots
deg = 3;
umax = 10;

[foundroots,~,poles] = nodePropagatorRoots(lens,umax)

% plot roots as circles
hold all
plot(foundroots,zeros(size(foundroots)),'o', poles, zeros(size(poles)),'*')
hold off