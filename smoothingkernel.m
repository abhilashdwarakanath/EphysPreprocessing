function [x,y] = smoothingkernel(duration,fs,sigma,type);

% function [x,y] = smoothingkernel(duration,sample_rate,sigma,type);
%
% This function outputs a kernel to smooth any signal using either a
% Gaussian kernel or an Alpha function kernel. Duration is the duration of
% your signal and is in SECONDS. Sigma is the standard-deviation of your
% kernel and is also in SECONDS. Sampling rate is in HERTZ. Type takes in
% the argument "gaussian" or "alpha".
%
% DO NOT NORMALISE THE KERNEL! Do it only if you are sure how it changes
% the firing rate values.
% Abhilash Dwarakanath. MPI biological cybernetics. October 2015.

% Check for arguments

if nargin<4
    error('Please put in 4 arguments viz kernel duration, sampling rate, standard deviation and the type')
end

% Create the time-axis

l = duration*fs;
x = (floor(-l / 2):floor(l / 2)) ./ fs;

% Check for kernel type and output kernel

switch type
    
    case 'gaussian'
        
        y = (1 / (sqrt(2 * pi) * sigma)) .* exp(-(x.^2 ./ (2 * sigma^2)));
        
    case 'alpha'
        
        h = heaviside(x); % The heaviside function elegantly defines the decay and cutoff time-point
        %y = (x./(sigma*sigma)).*(exp(-x./sigma)).*h;
        y = (1 / (sqrt(2 * pi) * sigma)) .* exp(-(x.^2 ./ (2 * sigma^2))).*h;
end

end
