%% Amath 582 Homework 4
% Brian de Silva
% 1422824

% This script is used to clean up and analyze the data extracted by
% extract_positions.m


%% Load and clean up data

testno = 4;         % Specify which test we are using

% Load the data produced by extract_positions
load(['test' num2str(testno) '.mat'])

% Max and min height should be the same for each camera so we normalize to
% [0, 1]
for i = 1:2
    indices1(i,:) = indices1(i,:) - min(indices1(i,:));
    indices1(i,:) = indices1(i,:) / max(indices1(i,:));
    
    indices2(i,:) = indices2(i,:) - min(indices2(i,:));
    indices2(i,:) = indices2(i,:) / max(indices2(i,:));
    
    indices3(i,:) = indices3(i,:) - min(indices3(i,:));
    indices3(i,:) = indices3(i,:) / max(indices3(i,:));
end



% % We want the bucket to start at the same position in each image (at the top of its trajectory) 
% figure()
% plot(indices1(1,2:end)), hold on
% plot(indices2(1,10:end))
% plot(indices3(2,1:end))
% xlabel('Frame'), ylabel('Row index'), title('Row indices for cameras over time (rescaled)')
% legend('Camera 1', 'Camera 2', 'Camera 3')

% Throw out some early data so that all three sources start at the same
% instant rather than at different times

switch testno
    case 1
        indices1 = indices1(:,9:end);
        indices2 = indices2(:,18:end);
        indices3 = indices3(:,8:end);
    case 2
        indices1 = indices1(:,12:end);
        indices2 = indices2(:,:);
        indices3 = indices3(:,12:end);
    case 3
        indices1 = indices1(:,12:end);
        indices2 = indices2(:,1:end);
        indices3 = indices3(:,4:end);
    otherwise
        indices1 = indices1(:,2:end);
        indices2 = indices2(:,10:end);
        indices3 = indices3(:,1:end);
end

% Normalize time to [0,1]
minsteps = min([size(indices1,2), size(indices2,2), size(indices3,2)]);
t = linspace(0,1,minsteps);

% We want the same number of time steps for each camera so we cut off all
% the sources as soon as one of the cameras stops recording
X = zeros(6,minsteps);
X(1,:) = indices1(1,1:minsteps);
X(2,:) = indices1(2,1:minsteps);
X(3,:) = indices2(1,1:minsteps);
X(4,:) = indices2(2,1:minsteps);
X(5,:) = indices3(1,1:minsteps);
X(6,:) = indices3(2,1:minsteps);

% Subtract off the mean for each row of data
% X = X - repmat(mean(X, 2), 1, nframes);
X = X - repmat(mean(X, 2), 1, minsteps);


% Compute SVD
% [U, S, V] = svd(X / sqrt(nframes-1));
[U, S, V] = svd(X / sqrt(minsteps-1));
lambda = diag(S);

% Plot percentage of variance captured by differnt modes
figure()
plot((lambda / sum(lambda)) * 100,'bo')
title('Percentage of variance captured by each mode')
xlabel('Mode')
ylabel('Percent variance captured')


% Project onto the left singular vectors
Y = U.' * X;

% Plot projections
figure()
plot(t,Y(1,:),t,Y(2,:),t,Y(3,:),t,Y(4,:),t,Y(5,:),t,Y(6,:))
legend('First','Second','Third','Fourth','Fifth','Sixth')
xlabel('Time')
ylabel('Normalized height')
title('Projections of the data onto the six left singular vectors')
