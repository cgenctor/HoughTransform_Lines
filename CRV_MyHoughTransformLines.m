%% CRV_11_MyHough
% name : Candas Genctor

%% clean up 
clear;
close all;
clc;
%% edge01
I = imread('edgetestimages/edge01.png');    % Load the edge image edge01.png
I = double(I);                              % Convert it to double
edgeImage = im2bw(I,1);                     % make it a binary image containing zeros and ones
[accu, h, alpha] = MyHough2(edgeImage);      % Perform a Hough transform
%% Visualization
figure();
subplot(131);
imshow(edgeImage);
title('edge image');
subplot(132);
imshow(accu,[]); % ,[] enables automatic scaling
axis square;
title('hough image');
subplot(133);
imshow(log(1+accu),[]);
axis square;
title('hough image (logarithmic scale)');
%% Determine the highest peak in the hough image and the corresponding values of h and alpha.
[x,y,~] = find(edgeImage);  % find edge coordinates
numOfEdgePoints = length(x);    % get the number of edge points
Nmax = numOfEdgePoints; % get Nmax entries. In this case all points are considered 
[ Avec, Ind ] = sort(accu(:),1,'descend');  % sort the accu values in descending order
max_values = Avec(1:Nmax);  % accu values are stored in descending order
[ ind_row, ind_col ] = ind2sub(size(accu),Ind(1:Nmax)); % fetch indices of accu values
rowcolPairs = [ind_row ind_col];    % to store the row col combinations of accu values

%% to suppress the values nearby peak points
A = rowcolPairs;
k = length(rowcolPairs);
nn = 1;
CC = zeros(5,5);
a = 1;
b = 1;
while k > 0 % this loop compares combinations and suppresses the ones near the peak
    for j = nn:length(rowcolPairs)-1
        CC(a,1:2) = rowcolPairs(b,:);
        CC(a,3:4) = rowcolPairs(j+1,:); %rowcolPairs(a+1,:); 
        CC(a,5) = sqrt((CC(a,4)-CC(a,2))^2+(CC(a,3)-CC(a,1))^2);
        if CC(a,5)< 5
            [~, locCC] = ismember(CC(a,3:4),rowcolPairs, 'rows');
            A(locCC,:) = 0;
        end
        a = a +1;
    end
    nn = nn+1;
    k = k -1;
    b = b +1;
end
% get rid of the nearby ones 
A( ~any(A,2), : ) = [];  %rows
% get the max desired number of hough peaks
N = 4;  % desired number of hough peaks 
houghPeaks = A(1:N,:);
% getting rho and alpha values for all corresponding lines in the figure
% and drawing all corresponding lines in the figure
figure()
imshow(edgeImage)
hold on
for mm = 1:N
    rhoVals(mm,1) =  h(houghPeaks(mm,1));
    if houghPeaks(mm,2) <= 90
        alphaVals(mm,1) = alpha(houghPeaks(mm,2) + 90);
    else
        alphaVals(mm,1) = alpha(houghPeaks(mm,2) - 90);
    end
    if alphaVals(mm,1) <= 0
        f = @(x,y) x * cos((alphaVals(mm,1)+1)*pi/180) + (-y)* sin((alphaVals(mm,1)+1)*pi/180) -rhoVals(mm,1);
        fimplicit(f,[0 150 0 150],'LineWidth',3)
        hold on
    else
        f = @(x,y) x * cos((alphaVals(mm,1)+1)*pi/180) + (-y)* sin((alphaVals(mm,1)+1)*pi/180) +rhoVals(mm,1);
        fimplicit(f,[0 150 0 150],'LineWidth',3)
        hold on
    end
end
%% clean up 
clear;
close all;
clc;
%% edge02
I = imread('edgetestimages/edge02.png');    % Load the edge image edge01.png
I = double(I);                              % Convert it to double
edgeImage = im2bw(I,1);                     % make it a binary image containing zeros and ones
[accu, h, alpha] = MyHough(edgeImage);      % Perform a Hough transform
%% Visualization
figure();
subplot(131);
imshow(edgeImage);
title('edge image');
subplot(132);
imshow(accu,[]); % ,[] enables automatic scaling
axis square;
title('hough image');
subplot(133);
imshow(log(1+accu),[]);
axis square;
title('hough image (logarithmic scale)');
%% Determine the highest peak in the hough image and the corresponding values of h and alpha.
[x,y,~] = find(edgeImage);  % find edge coordinates
numOfEdgePoints = length(x);    % get the number of edge points
Nmax = numOfEdgePoints; % get Nmax entries. In this case all points are considered 
[ Avec, Ind ] = sort(accu(:),1,'descend');  % sort the accu values in descending order
max_values = Avec(1:Nmax);  % accu values are stored in descending order
[ ind_row, ind_col ] = ind2sub(size(accu),Ind(1:Nmax)); % fetch indices of accu values
rowcolPairs = [ind_row ind_col];    % to store the row col combinations of accu values

%% to suppress the values nearby peak points
A = rowcolPairs;
k = length(rowcolPairs);
nn = 1;
CC = zeros(5,5);
a = 1;
b = 1;
while k > 0 % this loop compares combinations and suppresses the ones near the peak
    for j = nn:length(rowcolPairs)-1
        CC(a,1:2) = rowcolPairs(b,:);
        CC(a,3:4) = rowcolPairs(j+1,:); %rowcolPairs(a+1,:); 
        CC(a,5) = sqrt((CC(a,4)-CC(a,2))^2+(CC(a,3)-CC(a,1))^2);
        if CC(a,5)< 5
            [~, locCC] = ismember(CC(a,3:4),rowcolPairs, 'rows');
            A(locCC,:) = 0;
        end
        a = a +1;
    end
    nn = nn+1;
    k = k -1;
    b = b +1;
end
% get rid of the nearby ones 
A( ~any(A,2), : ) = [];  %rows
% get the max desired number of hough peaks
N = 4;  % desired number of hough peaks 
houghPeaks = A(1:N,:);
% getting rho and alpha values for all corresponding lines in the figure
% and drawing all corresponding lines in the figure
figure()
imshow(edgeImage)
hold on
for mm = 1:N
    rhoVals(mm,1) =  h(houghPeaks(mm,1));
    if houghPeaks(mm,2) <= 90
        alphaVals(mm,1) = alpha(houghPeaks(mm,2) + 90);
    else
        alphaVals(mm,1) = alpha(houghPeaks(mm,2) - 90);
    end
    if alphaVals(mm,1) <= 0
        f = @(x,y) x * cos((alphaVals(mm,1)+1)*pi/180) + (-y)* sin((alphaVals(mm,1)+1)*pi/180) -rhoVals(mm,1);
        fimplicit(f,[0 150 0 150],'LineWidth',3)
        hold on
    else
        f = @(x,y) x * cos((alphaVals(mm,1)+1)*pi/180) + (-y)* sin((alphaVals(mm,1)+1)*pi/180) +rhoVals(mm,1);
        fimplicit(f,[0 150 0 150],'LineWidth',3)
        hold on
    end
end
