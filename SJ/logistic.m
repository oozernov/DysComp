function y = logistic(params,x)

% y = logistic(params,x)
%
% K  - capacity --> = params(1)
% r  - growth rate --> params(2)
% t0 - reflection point --> params(3)
% x  - x values
%
% Use the following to calculate threshold
% thresh_level = (log(-thoi./(thoi-K))+r*t0)./r;
% with thoi - threshold definition, e.g. 0.5
% ------------------------------------
% adapted from B. Herrmann, Email: bherrmann@cbs.mpg.de
% 
% Sung-Joo Lim (sungjoo@bu.edu)
%

K  = 1; % capacity, or limiting value
r  = params(1); % slope, or growth rate
t0 = params(2); % inflection point, threshold
t  = x;

y = K ./ (1 + exp(-r*(t-t0)));
