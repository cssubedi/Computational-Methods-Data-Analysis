%% Clear Workspace
close all; clear variables; clc

%% Analysis on cropped-images
folders = dir("./Data/yalefaces_cropped/CroppedYale/yaleB*");
X = [];
for ind=1:length(folders)-1    % use last folder for test images
    files = dir([folders(ind,1).folder '/' folders(ind,1).name '/*.pgm']);
    
    for file=files'
        img = imread([file.folder '/' file.name]);
        img_vec = double(reshape(img,[],1));
        X = [X img_vec];
    end
end

X = X-repmat(mean(X,1),192*168,1);
[U,S,V] = svd(X, 'econ');
%% Show Eigenfaces
for i=1:5
    pause(0.1);
    imagesc(reshape(U(:,i),192,168)); colormap gray;
    set(gca,'xtick',[],'ytick',[])
    saveas(gcf, ['./Figures/cropped' num2str(i) '.png'])
end

%% Reconstruct an image
image = imread("./Data/yalefaces_cropped/CroppedYale/yaleB39/yaleB39_P00A+000E+45.pgm");
img_vec = double(reshape(image, [],1));
mean_img = mean(img_vec);
img_vec = img_vec - mean_img;   % Subtract mean

ranks = [25 100 200 400 1000];
for r=ranks
    img_gen = U(:,1:r)*U(:,1:r)'*img_vec + mean_img;
    imagesc(reshape(img_gen,192,168)); colormap gray;
    set(gca,'xtick',[],'ytick',[])
    title(['Rank - ' num2str(r) ' approximation.'])
    saveas(gcf, ['./Figures/r' num2str(r) '.png'])
    pause(0.2)
end
imagesc(reshape(img_vec,192,168)); colormap gray;
set(gca,'xtick',[],'ytick',[])
title(['Test Image'])
saveas(gcf, './Figures/test1.png')

%% Number of components needed to capture 50% of energy.
sigma = diag(S);
plot((sigma/sum(sigma))*100, 'b*')
energy = (sigma./sum(sigma)).*100;

e = 0;
r = 1;
while true
    if e < 50
        e = e + energy(r);
        r = r + 1;
    else
        break
    end
end

