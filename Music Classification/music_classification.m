%% Clear Workspace
close all; clear variables; clc

%% Global Variables
Test = 3;
Save = 0;

%% Load all music 
if Test == 1
    datasets = dir("./Data/test1/*.wav");
elseif Test == 2
    datasets = dir("./Data/test2/*.wav");
elseif Test == 3
    datasets = dir("./Data/test3/*.wav");
end

spg_full = [];
for data=datasets'
    [Y, Fs] = audioread([data.folder '\' data.name]);
    t = (1:length(Y))/Fs;
    Y = downsample(Y,2)';
    t = t(1:2:end);
    N = length(Y);
    L = t(end);
    
    if mod(N,2) == 0
        k = (2*pi)/(L)*[0:N/2-1 -N/2:-1];     % Even N
    else
        k = (2*pi)/(L)*[0:(N-1)/2 -(N-1)/2:-1];
    end
    ks = fftshift(k);
    tslide = [0:0.1:L];
    spg = [];
    
    % Gabor Tranform
    for ind = 1:length(tslide)
        gauss = exp(-100*(t-tslide(ind)).^2);     % width~0.25sec
        Yf = gauss.*Y;
        Yft = fft(Yf);
        spg = [spg; abs(fftshift(Yft))];
    end
    spg_full = [spg_full; spg];
end

% SVD of the total spectrogram
[U,S,V] = svd(spg_full', 'econ');    % 45 principal modes
%% Save the spectrogram
% if Save
%     save(['./Data/test' num2str(Test) '/spectrogram.mat'], 'spg_full');
% end

%% Split the data into training and testing
V2 = V(:,1:5); % use first 5 modes
classification_errors = [];

for n=[1:100]
    group1_train = randsample(15,10);   % Training set of 30
    group2_train = randsample(15,10)+15;
    group3_train = randsample(15,10)+30;
    
    group1_test = setdiff(1:15,group1_train);   % Testing set of 15
    group2_test = setdiff(16:30,group2_train);
    group3_test = setdiff(31:45,group3_train);

    train = [V2(group1_train,:); V2(group2_train,:); V2(group3_train,:)];
    test = [V2(group1_test,:); V2(group2_test,:); V2(group3_test,:)];
    
    train_classification = [ones(10,1); 2*ones(10,1); 3*ones(10,1)];
    test_classification = [ones(5,1); 2*ones(5,1); 3*ones(5,1)];
    
    naive_base = fitcnb(train, train_classification);
    prediction = naive_base.predict(test);
    
    bar(prediction);
    set(gca, 'FontSize',14);        
    xlabel('Test music'); ylabel('Predicted Classification');
    drawnow
    pause(0.1)
    
    error = 0;
    for ind=1:length(prediction)
        if prediction(ind) ~= test_classification(ind)
            error = error + 1;
        end
    end
    error = (error/15)*100;
    classification_errors = [classification_errors error];
end

mean_error = mean(classification_errors);
std_error = std(classification_errors);

%% Generate plots for the paper
% Energy captured in each POD modes
figure();
sigma = diag(S);
plot(sigma/(sum(sigma))*100,'-bo')
semilogy(sigma/(sum(sigma))*100,'-bo')
xlabel('Mode'); ylabel('Percentage of energy in each mode.')

% True classification bar graph
figure()
bar([ones(5,1); 2*ones(5,1); 3*ones(5,1)])
set(gca, 'FontSize',14);        
xlabel('Test music'); ylabel('True Classification');

% First 3 principal components
figure();
plot3(V(1:15,1), V(1:15,2), V(1:15,3), 'o', 'MarkerFaceColor', 'r')
hold on
plot3(V(16:30,1), V(16:30,2), V(16:30,3), 'o', 'MarkerFaceColor', 'b')
hold on
plot3(V(31:45,1), V(31:45,2), V(31:45,3), 'o', 'MarkerFaceColor', 'k')
hold on
legend('Pear Jam', 'Soundgarden', 'Nirvana')
grid on

%% Spectrogram of Nirvana music
[Y, Fs] = audioread("./Data/test2/Nirvana.wav");
t = (1:length(Y))/Fs;
Y = downsample(Y,2)';
t = t(1:2:end);
N = length(Y);
L = t(end);

if mod(N,2) == 0
    k = (2*pi)/(L)*[0:N/2-1 -N/2:-1];     % Even N
else
    k = (2*pi)/(L)*[0:(N-1)/2 -(N-1)/2:-1];
end
ks = fftshift(k);
tslide = [0:0.1:L];
spg = [];

% Gabor Tranform
for ind = 1:length(tslide)
    gauss = exp(-100*(t-tslide(ind)).^2);     % width~0.25sec
    Yf = gauss.*Y;
    Yft = fft(Yf);
    spg = [spg; abs(fftshift(Yft))];
end

figure();
pcolor(tslide, ks, spg.')
shading interp
colormap(hot)
xlabel("Time [sec]")
ylabel("Frequency [Hz]")
