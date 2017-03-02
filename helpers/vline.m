function vline(x,varargin)
%% vline(x)
% Plots vertical line on current axes at X=x. (x can be a vector)
% Default is a black dashed line, but user can specify formatting

hold on
if nargin < 2,
    for i=1:length(x)
        plot([x(i),x(i)],get(gca,'YLim'),'--k')
    end
else
    for i=1:length(x)
        plot([x(i),x(i)],get(gca,'YLim'),varargin{:})
    end
end