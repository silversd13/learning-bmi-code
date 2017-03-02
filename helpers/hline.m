function hline(x,varargin)
%% hline(x)
% Plots horizontal line on current axes at X=x. (x can be a vector)
% Default is a black dashed line, but user can specify formatting

hold on
if nargin < 2,
    for i=1:length(x)
        plot(get(gca,'XLim'),[x(i),x(i)],'--k')
    end
else
    for i=1:length(x)
        plot(get(gca,'XLim'),[x(i),x(i)],varargin{:})
    end
end