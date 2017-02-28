function [beta,learning_fun] = calc_learning_rate(timeToTarget)

%% inputs
narginchk(1,1)
if ~isvector(timeToTarget),
    error('timeToTarget should be a vector')
end

%% outputs
nargoutchk(2,2)

%% fit exponential learning curve to behavior
ls_opt = optimset('Display','off');
learning_fun = @(beta,x) (beta(1))*exp(-1*x*beta(2)) + beta(3);
beta = lsqcurvefit(learning_fun,[10,.1,2],(1:length(timeToTarget))',timeToTarget,[0,0,0],[],ls_opt);

x = 1:length(timeToTarget);
yhat = learning_fun(beta,x);

%% plot
hold on
plot(x,timeToTarget,'r^')
plot(x,yhat,'--k')
text(mean(get(gca,'XLim')),mean(get(gca,'YLim')),...
    sprintf('y = %.2f*exp^{-%.2f*x} + %.2f',beta(1),beta(2),beta(3)),...
    'FontSize',18,'HorizontalAlignment','Center')
xlabel('trial')
ylabel('reward time (s)')
title('BMI Learning Curve')
