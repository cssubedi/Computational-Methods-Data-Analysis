%% Clear workspace
close all; clear variables; clc

%% Set up Global variables
Test = 4;   % 1, 2, 3 or 4
Save = 1;
Filter = 0;
Verbose = 1;

%% Load Test data
load(['./Data/Test' num2str(Test) '/cam1_' num2str(Test) '.mat']);
load(['./Data/Test' num2str(Test) '/cam2_' num2str(Test) '.mat']);
load(['./Data/Test' num2str(Test) '/cam3_' num2str(Test) '.mat']);

if Test == 1
    window = 15;
    camera1 = vidFrames1_1;
    camera2 = vidFrames2_1;
    camera3 = vidFrames3_1;
    X1=321; Y1=224;
    X2=275; Y2=275;
    X3=319; Y3=272;
elseif Test == 2
    window = 20;
    camera1 = vidFrames1_2;
    camera2 = vidFrames2_2;
    camera3 = vidFrames3_2;
    X1=325; Y1=308;
    X2=313; Y2=360;
    X3=348; Y3=245;
elseif Test == 3
    window = 20;
    camera1 = vidFrames1_3;
    camera2 = vidFrames2_3;
    camera3 = vidFrames3_3;
    X1=320; Y1=286;
    X2=238; Y2=292;
    X3=352; Y3=231;
elseif Test == 4
    window = 20;
    camera1 = vidFrames1_4;
    camera2 = vidFrames2_4;
    camera3 = vidFrames3_4;
    X1=423; Y1=325;
    X2=220; Y2=321;
    X3=412; Y3=182;
end

[m1,n1,~,t1] = size(camera1);
[m2,n2,~,t2] = size(camera2);
[m3,n3,~,t3] = size(camera3);

%% Check the data
if Verbose
    for j=1:t1
        imshow(camera1(:,:,:,j))
        drawnow;
    end
    for j=1:t2
        imshow(camera2(:,:,:,j))
        drawnow;
    end
    for j=1:t3
        imshow(camera3(:,:,:,j))
        drawnow;
    end
end

%% Find position of the bucket
% Camera 1
position1 = zeros(2, t1);
position1(:,1) = [X1; Y1];

for frame = 2:t1
    img = camera1(:,:,:,frame);
    X = position1(1,frame-1); Y = position1(2, frame-1);
    roi = crop(img, X, Y, window);
    mask = (sum(roi, 3) == max(max(sum(roi,3))));
    [r, c] = find(mask);
    position1(:,frame) = [X - window + c(1); Y - window + r(1)];
end 
figure(1)
plot(1:t1, position1(2,:))

% Camera 2
position2 = zeros(2, t2);
position2(:,1) = [X2; Y2];

for frame = 2:t2
    img = camera2(:,:,:,frame);
    X = position2(1,frame-1); Y = position2(2, frame-1);
    roi = crop(img, X, Y, window);
    mask = (sum(roi, 3) == max(max(sum(roi,3))));
    [r, c] = find(mask);
    position2(:,frame) = [X - window + c(1); Y - window + r(1)];
end 
figure(2)
plot(1:t2, position2(2,:))

% Camera 3
position3 = zeros(2, t3);
position3(:,1) = [X3; Y3];

for frame = 2:t3
    img = camera3(:,:,:,frame);
    X = position3(1,frame-1); Y = position3(2, frame-1);
    roi = crop(img, X, Y, window);
    mask = (sum(roi, 3) == max(max(sum(roi,3))));
    [r, c] = find(mask);
    position3(:,frame) = [X - window + c(1); Y - window + r(1)];
end 
figure(3)
plot(1:t3, position3(2,:))
figure(4)
plot3(position3(1,:), position3(2,:), 1:t3)

%% Plot position onto image
if Verbose
    figure()
    for frame = 1:t1
        imshow(camera1(:,:,:,frame))
        hold on
        plot(position1(1,frame), position1(2,frame), 'ro', 'Linewidth', 3)
        drawnow
    end

    figure()
    for frame = 1:t2
        imshow(camera2(:,:,:,frame))
        hold on
        plot(position2(1,frame), position2(2,frame), 'ro', 'Linewidth', 3)
        drawnow
    end

    figure()
    for frame = 1:t3
        imshow(camera3(:,:,:,frame))
        hold on
        plot(position3(1,frame), position3(2,frame), 'ro', 'Linewidth', 3)
        drawnow
    end
end

%% Shannon Filter the position data (Did not give much improvement!)
if Filter
    position1 = position1(:,10:210);
    position2 = position2(:,20:220);
    position3 = position3(:,10:210);

    xposition1_t = fft(position1(1,:));
    xposition1_t(10:end-10) = 0;
    position1(1,:) = abs(ifft(xposition1_t));

    yposition1_t = fft(position1(2,:));
    yposition1_t(10:end-10) = 0;
    position1(2,:) = abs(ifft(yposition1_t));

    xposition2_t = fft(position2(1,:));
    xposition1_t(10:end-10) = 0;
    position2(1,:) = abs(ifft(xposition2_t));

    yposition2_t = fft(position2(2,:));
    yposition1_t(10:end-10) = 0;
    position2(2,:) = abs(ifft(yposition2_t));

    xposition3_t = fft(position3(1,:));
    xposition3_t(10:end-10) = 0;
    position3(1,:) = abs(ifft(xposition3_t));

    yposition3_t = fft(position3(2,:));
    yposition3_t(10:end-10) = 0;
    position3(2,:) = abs(ifft(yposition3_t));

    if Verbose
        figure()
        for frame = 1:t1
            imshow(camera1(:,:,:,frame))
            hold on
            plot(position1(1,frame), position1(2,frame), 'ro', 'Linewidth', 3)
            drawnow
        end

        figure()
        for frame = 1:t2
            imshow(camera2(:,:,:,frame))
            hold on
            plot(position2(1,frame), position2(2,frame), 'ro', 'Linewidth', 3)
            drawnow
        end

        figure()
        for frame = 1:t3
            imshow(camera3(:,:,:,frame))
            hold on
            plot(position3(1,frame), position3(2,frame), 'ro', 'Linewidth', 3)
            drawnow
        end
    end
    figure(); hold on
    plot(position1(2,:))
    plot(position2(2,:))
    plot(position3(1,:))
end

%% Save the position data
if Save
    save(['./Data/Test' num2str(Test) '/position.mat'], 'position1', 'position2', 'position3');
end



