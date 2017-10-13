% Amath 582 Homework 4
% Brian de Silva
% 142824

% This script can be used to extract position data from the bucket
% recordings

%% Extract locations of buckets for analysis
testno = 4;     % Specify which test we are using

% Test 1: Ideal signal
% Test 2: Noisy signal
% Test 3: Horizontal displacement
% Test 4: Horizontal displacement and rotation


% Load the video data if it's not already loaded
if ~exist(['vidFrames1_' num2str(testno)],'var')
    load(['cam1_' num2str(testno) '.mat'])
end
if ~exist(['vidFrames2_' num2str(testno)],'var')
    load(['cam2_' num2str(testno) '.mat'])
end
if ~exist(['vidFrames3_' num2str(testno)],'var')
    load(['cam3_' num2str(testno) '.mat'])
end

% Give all video matrices the same name
switch testno
    case 1
        vidFrames1 = vidFrames1_1;
        vidFrames2 = vidFrames2_1;
        vidFrames3 = vidFrames3_1;
    case 2
        vidFrames1 = vidFrames1_2;
        vidFrames2 = vidFrames2_2;
        vidFrames3 = vidFrames3_2;
    case 3
        vidFrames1 = vidFrames1_3;
        vidFrames2 = vidFrames2_3;
        vidFrames3 = vidFrames3_3;
    otherwise
        vidFrames1 = vidFrames1_4;
        vidFrames2 = vidFrames2_4;
        vidFrames3 = vidFrames3_4;
end

%% Get location of bucket at each time step:

% Prime first position by actual inspection then only check for new
% position within a window of pixels near previous position
indices1 = zeros(2,size(vidFrames1,4));
figure()
imshow(vidFrames1(:,:,:,1)), title('Please click the flashlight at the top of the bucket')
[x1, y1] = ginput(1);
indices1(:,1) = [y1; x1];

indices2 = zeros(2,size(vidFrames2,4));
figure()
imshow(vidFrames2(:,:,:,1)), title('Please click the flashlight at the top of the bucket')
[x2, y2] = ginput(1);
indices2(:,1) = [y2; x2];

indices3 = zeros(2,size(vidFrames3,4));
figure()
imshow(vidFrames3(:,:,:,1)), title('Please click the flashlight at the top of the bucket')
[x3, y3] = ginput(1);
indices3(:,1) = [y3; x3];

close all

% temp = vidFrames3_1(:,:,:,1);
% figure(), imshow(temp)
% temp_filtered = (sum(temp,3)>=760);
% figure(), imshow(255*temp_filtered)
% local = temp((indices3(1,1)-10):(indices3(1,1)+10),(indices3(2,1)-10):(indices3(2,1)+10),:);
% figure(), imshow(local)



% Set parameters for first camera
rwindsize1 = 15;                         % Size of row window about previously known bucket location about which to look for new location
cwindsize1 = 15;                          % Size of column window about previously known bucket location about which to look for new location
local1 = uint8(zeros(rwindsize1+1,cwindsize1+1,3));

% Set parameters for second camera
rwindsize2 = 15;                         % Size of row window about previously known bucket location about which to look for new location
cwindsize2 = 15;                          % Size of column window about previously known bucket location about which to look for new location
rgb2 = [170,255,255];                % RGB values near the flashlight in test 2.2
local2 = uint8(zeros(rwindsize2+1,cwindsize2+1,3));

% Set parameters for third camera
rwindsize3 = 15;                             % Size of row window about previously known bucket location about which to look for new location
cwindsize3 = 15;                          % Size of column window about previously known bucket location about which to look for new location
local3 = uint8(zeros(rwindsize3+1,cwindsize3+1,3));

% Temporary window size variables in case window extends outside the
% image
temprsize = 0;
tempcsize = 0;

    % --------------First camera-----------------
for k=2:size(vidFrames1,4)
    row1 = indices1(1,k-1);            % Previous row location of bucket
    col1 = indices1(2,k-1);             % Previous column location of bucket
    
    % Make sure window doesn't extend outside image
    if (row1 - rwindsize1) < 1
        temprsize = row1-1 + (row1==1);
    elseif (row1 + rwindsize1) > size(vidFrames1,1)
        temprsize = size(vidFrames1,1) - row1 + (row1 == size(vidFrames1,1));
    else
        temprsize = rwindsize1;
    end
    if (col1 - cwindsize1) < 1
        tempcsize = col1-1 + (col1==1);
    elseif (col1 + cwindsize1) > size(vidFrames1,1)
        tempcsize = size(vidFrames1,1) - col1 + (col1 == size(vidFrames1,1));
    else
        tempcsize = cwindsize1;
    end
        
    
    % Look in a small window about previous bucket location
    local1 = vidFrames1((row1-temprsize):(row1+temprsize),(col1-tempcsize):(col1+tempcsize),:,k);
    
    % Find pixel with maximum intensity
    local_filtered1 = (sum(local1,3)==max(max(sum(local1,3))));
    
    % Use top-left most pixel for location
    [r,c] = find(local_filtered1);
    indices1(1,k) = r(1) + (row1-temprsize)-1;
    indices1(2,k) = c(1) + (col1-tempcsize)-1;
end


% --------------Second camera-----------------
for k=2:size(vidFrames2,4)
    row2 = indices2(1,k-1);            % Previous row location of bucket
    col2 = indices2(2,k-1);             % Previous column location of bucket
    
    % Make sure window doesn't extend outside image
    if (row2 - rwindsize2) < 1
        temprsize = row2-1 + (row2==1);
    elseif (row2 + rwindsize2) > size(vidFrames2,1)
        temprsize = size(vidFrames2,1) - row2 + (row2 == size(vidFrames2,1));
    else
        temprsize = rwindsize2;
    end
    if (col2 - cwindsize2) < 1
        tempcsize = col2-1 + (col2==1);
    elseif (col2 + cwindsize2) > size(vidFrames2,1)
        tempcsize = size(vidFrames2,1) - col2 + (col2 == size(vidFrames2,1));
    else
        tempcsize = cwindsize2;
    end
    
    % Look in a small window about previous bucket location
    local2 = vidFrames2((row2-temprsize):(row2+temprsize),(col2-tempcsize):(col2+tempcsize),:,k);
    
    % Find pixel with maximum intensity (whitest pixel)
    local_filtered2 = (sum(local2,3)==max(max(sum(local2,3))));
    
    % Use top-left most pixel for location
    [r,c] = find(local_filtered2);
    indices2(1,k) = r(1) + (row2-temprsize)-1;
    indices2(2,k) = c(1) + (col2-tempcsize)-1;
end


% --------------Third camera-----------------
for k=2:size(vidFrames3,4)
    row3 = indices3(1,k-1);            % Previous row location of bucket
    col3 = indices3(2,k-1);             % Previous column location of bucket
    
    % Make sure window doesn't extend outside image
    if (row3 - rwindsize3) < 1
        temprsize = row3-1 + (row3==1);
    elseif (row3 + rwindsize3) > size(vidFrames3,1)
        temprsize = size(vidFrames3,1) - row3 + (row3 == size(vidFrames3,1));
    else
        temprsize = rwindsize3;
    end
    if (col3 - cwindsize3) < 1
        tempcsize = col3-1 + (col3==1);
    elseif (col3 + cwindsize3) > size(vidFrames3,1)
        tempcsize = size(vidFrames3,1) - col3 + (col3 == size(vidFrames3,1));
    else
        tempcsize = cwindsize3;
    end
    
    % Look in a small window about previous bucket location
    local3 = vidFrames3((row3-temprsize):(row3+temprsize),(col3-tempcsize):(col3+tempcsize),:,k);
    
    % Find pixel with maximum intensity (whitest pixel)
    local_filtered3 = (sum(local3,3)==max(max(sum(local3,3))));
    
    % Use top-left most pixel for location
    [r,c] = find(local_filtered3);
    indices3(1,k) = r(1) + (row3-temprsize)-1;
    indices3(2,k) = c(1) + (col3-tempcsize)-1;
end

% % Play the video
% for k=1:size(vidFrames1_1,4)
%     imshow(vidFrames1_1(:,:,:,k));
%     title(['t= ' num2str(k)])
%     drawnow
% end

% % Sanity check: plot locations the computer found against true video
% figure()
% for k=1:size(vidFrames2,4)
%     temp = vidFrames2(:,:,:,k);
%     temp((indices2(1,k)-3):(indices2(1,k)+3),(indices2(2,k)-3):(indices2(2,k)+3),:)=0;
%     imshow(temp), axis on
%     title(['t= ' num2str(k)])
%     drawnow
% end

% Phew, it worked!

% Plot indices corresponding to height
figure()
plot(indices1(1,:)), hold on
plot(indices2(1,:))
plot(indices3(2,:))
xlabel('Frame'), ylabel('Row index'), title('Row indices for cameras over time')
legend('Camera 1', 'Camera 2', 'Camera 3')

% Save the position measurements
save(['test' num2str(testno) '.mat'], 'indices1', 'indices2', 'indices3');




