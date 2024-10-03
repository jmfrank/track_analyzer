function output = scattered1D(x,y,xq, varargin)
%% SCATTERED1D
% inerpolate 1D scattered data using simple binned averaging.  
% Inputs x, y must be column vectors of the same size. xq are the query
% points. You can optionally input an interpolation / extrapolation method:
%       scattered1D(x,y,xq,'pchip')
% Note that 'linear' is the default. You can also input an interpolation
% resolution (used for function fitting) as follows:
%         scattered1D(x,y,xq,200)
% Otherwise, the function will chooe a resolution for you. IF you want to
% see a plot of the data, switch the 'do_plotting' flag to true in the
% code.
%
% Here is some sample code to test the function. Just uncomment and press
% Ctrl + Enter:
%     var = 1; nargin = 3;
%     [X,Y] = meshgrid(-5:0.01:5, -1:0.05:1);
%     Z = 1.1779.*exp(-0.5.*X.^2./var).*exp(-0.5.*Y.^2./var) + rand(size(X)).*0.1;
%     x = X(1:end)'; y = Z(1:end)';
%     xq = diff([min(x), max(x)]);
%     xq = [min(x) - xq/50, rand([1,25]).*xq + min(x), max(x) + xq/50];


do_plotting = false

interp1_method = 'linear'; % the default
sugested_numpts = [];
if nargin > 3
    for n = 1:nargin-3
        if ischar(varargin{n})
            if (strcmp(varargin{n}, 'linear') || strcmp(varargin{n}, 'nearest') || strcmp(varargin{n}, 'next') || strcmp(varargin{n}, 'previous') || strcmp(varargin{n}, 'pchip') || strcmp(varargin{n}, 'cubic') || strcmp(varargin{n}, 'v5cubic') || strcmp(varargin{n}, 'makima') || strcmp(varargin{n}, 'spline'))
                interp1_method = varargin{n};
            else
                error('ERROR (scattered1D): interpolation method should match interp1 inputs'); 
            end
        elseif isnumeric(varargin{n}) && varargin{n} > 0
            sugested_numpts = varargin{n};
        else
            error('ERROR (scattered1D): unable to interpret extra inputs'); 
        end
    end
end

range = [min(x), max(x)];
if isempty(sugested_numpts)
    % extract the most common x differential value, to help with estimating resolution
    if numel(x) <2 || ~iscolumn(x) || ~isequal(size(x),size(y)), error('ERROR (scattered1D): x,y inputs should be column vectors of the same size'); end
    rd = sortrows([x,y]);
    mode_is = mean(diff(rd(:,1)));
    
    % figure out interpolation resolution and adjust for smoothing
    smoothing_rate = 100; % this should be adjusted in a more sophisticated way
    sugested_numpts = ceil(diff(range)/mode_is/smoothing_rate);
    if sugested_numpts < 10; sugested_numpts = 10; end
end

% smooth function
fitted_fcn = linspace(range(1), range(2), sugested_numpts)';
for n = 1:numel(fitted_fcn)-1
    % extract relevant values
    if n == 1
        relevant_values = (x <= fitted_fcn(2,1));
    else
        relevant_values = (x > fitted_fcn(n,1)) & (x <= fitted_fcn(n+1,1));
    end

    if any(relevant_values)
        fitted_fcn(n,2) = mean(y(relevant_values));
    else
        fitted_fcn(n,2) = nan;
    end
end

% adjust domain to centre of averaging intervals
fitted_fcn(:,1) = fitted_fcn(:,1) + diff(fitted_fcn(1:2,1));
fitted_fcn(end,:) = [];

% now produce the output, with extrapolation
output = interp1(fitted_fcn(:,1), fitted_fcn(:,2),xq,interp1_method, 'extrap');

% plotting
if do_plotting == true
    figure, 
    hold on, 
    plot(fitted_fcn(:,1), fitted_fcn(:,2))
    plot(x,y,'color', 'k', 'LineStyle','none', 'Marker','.','MarkerSize',0.5)
    plot(xq,output,'gx')
end
