alpha = 30; 
beta = 0.8; 
c = 3; 

th = 0.001:0.001:10;
r= @(t) (alpha - alpha./(1+beta*t));

f = @(t)(c*t);

% subplot(311)
plot(th, r(th),'Linewidth',1.5);

hold on 
% subplot(312)
plot(th, -1*f(th),'Linewidth',1.5)
xlabel('Time (s)');
ylabel('Cost of harvesting in Patch');

hold on
yline(0)

beautifyfig

figure
% subplot(313)
plot(th, r(th)-f(th),'Linewidth',1.5)
xlabel('Time (s)');
ylabel('Net Harvest in Patch');

beautifyfig