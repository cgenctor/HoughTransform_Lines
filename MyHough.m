   function [ accu, h, alpha ] = MyHough( edgeImage )
%MYHOUGH detects straight lines vie Hough transform
%   [ accu, h, alpha ] = MyHough( edgeImage ) calculates the Hough
%   transform of the binary edgeImage. The hough image is contained in
%   accu. h is a vector containing the distance values. With D being the
%   diagonal length of the edgeImage, the values in h range from -D to D.
%   alpha contains the angle values in degree. These range from -90 to 89.

edgeImage = double(edgeImage);  
[M,N] = size(edgeImage);    % get the size of edge image for diagonal calculation
dalpha = 1;                 % change in alpha
dh = 1;                     % change in rho
alpha = -90:dalpha:89;         % alpha changes from -90 to 89
n_alpha = length(alpha);    % number of alphas 
D = sqrt ((M-1)^2+(N-1)^2); % diagonal calculation
q = ceil(D/dh);             % rounding diagonal
h = -q*dh:1:q*dh; % linspace(-q*dh,q*dh,n_h);            % rho changes from -D to D
%h_flipped = flip(h);
n_h =  length(h);             % number of rhos
[x,y,~] = find(edgeImage); % find x,y coordinates  and corresponding value of edge pixels
% x = x-1;
% y = y-1;
numOfEdgePoints = length(x);
% Create the hough akkumulator.
accu = zeros (n_h,n_alpha); % initialize output with zeros
%% Implement the hough transform 
for i = 1 :  numOfEdgePoints   %For every edge point in edge image
    edgePx = x(i,1);            %Determine the coordinates of the edge point
    edgePy = y(i,1);
    rho = edgePx * cos(alpha(1,:)*pi/180) + (-edgePy)* sin(alpha(1,:)*pi/180) ;
    for ii=1:length(rho)
        [ ~, idx ] = min( abs( h - rho(ii)) );
        % accu(idx,181-ii) = accu(idx,181-ii) + 1;
        accu(idx,181-ii) = accu(idx,181-ii) + 1;
    end
end
end
