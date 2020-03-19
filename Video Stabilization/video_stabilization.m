%% Clear Workspace
close all; clear variables; clc

%% Global Variables
Verbose = 0;

%% Load the video and construct the matrix
vidReader = VideoReader('video.mp4','CurrentTime',1130);

% Using Horn-Schunck method
opticFlow = opticalFlowHS;

% Plotting axis for optical flow
h = figure;
movegui(h);
hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
hPlot = axes(hViewPanel);
set(h,'DefaultLineLineWidth',10)

% Initialize matrices
FX = [];
FY = [];

nFrames = 0;
while hasFrame(vidReader) && nFrames < 300
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);  
    frameResize = imresize(frameGray, 0.25);
    flow = estimateFlow(opticFlow,frameResize);
    
    f_x = flow.Vx;
    f_y = flow.Vy;
    FX = [FX double(reshape(f_x, [],1))];
    FY = [FY double(reshape(f_y, [], 1))];
    
    imshow(frameResize)
    hold on
    plot(flow,'DecimationFactor',[10 10],'ScaleFactor',100,'Parent',hPlot);
    hold off
    pause(0.001)
    nFrames = nFrames + 1;
end

%% Singular Value Decomposition
% Zero-mean the matrices before SVD
X_mean = repmat(mean(FX,1),129600,1);
FX = FX-X_mean;
Y_mean = repmat(mean(FY,1),129600,1);
FY = FY-Y_mean;

[UX,SX,VX] = svd(FX, 'econ');
[UY,SY,VY] = svd(FY, 'econ');

%% Plot principal modes
figure()
for i=1
    imagesc(reshape(UY(:,i),270,480)); colormap gray;
    set(gca,'xtick',[],'ytick',[])
    pause(2);
end

% Plot energy content of each singular modes
sigma_y = diag(SY);
figure()
plot((sigma_y/sum(sigma_y))*100, 'r*')

%% Extra plotting for the paper
if Verbose
    VX_p = UX(:,1:10)*SX(1:10,1:10)*VX(:,1:10)';
    VY_p = UX(:,1:10)*SX(1:10,1:10)*VX(:,1:10)';

    VX_p = VX_p + X_mean;
    VY_p = VY_p + Y_mean;

    [X,Y] = meshgrid(1:480,1:270);
    X_p = X-reshape(VX_p(:,200), 270, 480);
    Y_p = Y-reshape(VY_p(:,200), 270, 480);

    % Plot final flow velocity after warping
    figure()
    q = quiver(X_p,Y_p);
    q.AutoScale = 'on';
    q.AutoScaleFactor = 2;

    % Plot flow field of 200th frame.
    figure()
    q3 = quiver(reshape(FX(:,200), 270, 480),...
        reshape(FY(:,200), 270, 480), 'linewidth', 5);
    q3.AutoScale = 'on';
    q3.AutoScaleFactor = 2;
    axis off

    % Plot 10-rank approximation of the flow field of 200th frame
    figure()
    q4 = quiver(reshape(VX_p(:,200), 270, 480),...
        reshape(VY_p(:,200), 270, 480));
    q4.AutoScale = 'on';
    q4.AutoScaleFactor = 2;
    axis off
end

%% Reconstruct the video using low-rank approximation
% 20-rank approximation of flow vield matrices.
VX_p = UX(:,1:20)*SX(1:20,1:20)*VX(:,1:20)' + X_mean;
VY_p = UX(:,1:20)*SX(1:20,1:20)*VX(:,1:20)' + Y_mean;

[X,Y] = meshgrid(1:480,1:270);

vidReader = VideoReader('video.mp4','CurrentTime',1130);
nFrames = 1;

% Read each image frame
while hasFrame(vidReader) && nFrames <= 300
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);  
    frameResize = imresize(frameGray, 0.25);
    
    % Find previous pixel locations using flow velocity
    X_p = X-reshape(VX_p(:,nFrames), 270, 480);
    Y_p = Y-reshape(VY_p(:,nFrames), 270, 480);
    
    % Warp the image
    image = interp2(im2double(frameResize), X_p, Y_p);
    imshow(image)
    pause(0.1)
    
    % Save the image
    imwrite(image, ['./Result/image' num2str(nFrames) '.png'])
    nFrames = nFrames + 1;
end

%% Check reconstructed images
files = dir("./Result/image*");
for file=files'
    img = imread([file.folder '/' file.name]);
    imshow(img)
    pause(0.1)
end

%% Observe each column of the matrix
if Verbose
    figure()
    for i=1:300
        imagesc(reshape(FX(:,i),270,480)); colormap gray;
        set(gca,'xtick',[],'ytick',[])
        pause(2);
    end
end
