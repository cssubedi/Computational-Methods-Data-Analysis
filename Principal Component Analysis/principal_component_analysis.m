%% Clear Workspace
close all; clear variables; clc

%% Global Variables
Test = 1;
Verbose = 1;

%% Load position data
data = load(['./Data/Test' num2str(Test) '/position.mat']);

if Verbose
    figure(); hold on
    plot(data.position1(2,:)),
    plot(data.position2(2,:))
    plot(data.position3(1,:))
    xlabel('Frame'), ylabel('Y position')
    legend('Camera 1', 'Camera 2', 'Camera 3')
end
%% Generate same size position vector for each camera
% Extract region of clean data.
if Test == 1
    data.position1 = data.position1(:,10:210);
    data.position2 = data.position2(:,20:220);
    data.position3 = data.position3(:,10:210);
    size=201;
elseif Test == 2
    data.position1 = data.position1(:,10:210);
    data.position2 = data.position2(:,1:201);
    data.position3 = data.position3(:,15:215);
    size=201;
elseif Test == 3
    data.position1 = data.position1(:,18:218);
    data.position2 = data.position2(:,5:205);
    data.position3 = data.position3(:,10:210);
    size=201;
elseif Test == 4
    data.position1 = data.position1(:,2:352);
    data.position2 = data.position2(:,10:360);
    data.position3 = data.position3(:,1:351);
    size=351;
end

time = linspace(0,1,size);

if Verbose
    figure(); hold on
    plot(1:size,data.position1(2,:)),
    plot(1:size,data.position2(2,:))
    plot(1:size,data.position3(1,:))
    xlim([0 size])
    xlabel('Frame'), ylabel('X/Y position')
    legend('Camera 1', 'Camera 2', 'Camera 3')
end

if Verbose
    figure(); hold on
    plot(1:size,data.position1(1,:)),
    plot(1:size,data.position2(1,:))
    plot(1:size,data.position3(2,:))
    xlim([0 size])
    xlabel('Frame'), ylabel('Y position')
    legend('Camera 1', 'Camera 2', 'Camera 3')
end

%% Generate position data normalized to [0,1]
for i=[1,2]
    data.position1(i,:) = (data.position1(i,:) - min(data.position1(i,:)))/ ...
        (max(data.position1(i,:)) - min(data.position1(i,:)));
    data.position2(i,:) = (data.position2(i,:) - min(data.position2(i,:)))/ ...
        (max(data.position2(i,:)) - min(data.position2(i,:)));
    data.position3(i,:) = (data.position3(i,:) - min(data.position3(i,:)))/ ...
        (max(data.position3(i,:)) - min(data.position3(i,:)));
end

if Verbose
    figure(); hold on
    plot(1:size,data.position1(2,:))
    plot(1:size,data.position2(2,:))
    plot(1:size,data.position3(1,:))
    xlim([0 size])
    xlabel('Frame'), ylabel('X/Y position')
    legend('Camera 1', 'Camera 2', 'Camera 3')
end

%% Construct the matrix with zero-mean position data
A = [data.position1(1,:); data.position1(2,:); data.position2(1,:); ...
    data.position2(2,:); data.position3(1,:); data.position3(2,:)];
A = A-repmat(mean(A,2),1,size);
% A = A./repmat(std(A,0,2),1,size);


%% Perform the SVD on matrix A

[U, S, V] = svd(A); 
sigma = diag(S);

% Energy captured in each POD modes
figure();
plot(sigma/(sum(sigma))*100,'-bo')
xlabel('Mode'); ylabel('Sigma')

% Projection of modes onto the left singular vectors
B = U.'* A;
figure()
plot(time,B(1,:),time,B(2,:),time,B(3,:),time,B(4,:),time,B(5,:),time,B(6,:))
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6')
xlabel('Time'); ylabel('Normalized position')







