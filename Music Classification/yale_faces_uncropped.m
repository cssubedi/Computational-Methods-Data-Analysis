%% Clear Workspace
close all; clear variables; clc

%% Analysis on Uncropped-images
files = dir("./Data/yalefaces/subject*");
X = [];
for file=files(1:end-1,1)'
    img = imread([file.folder '/' file.name]);
    imshow(img)
    img_vec = double(reshape(img,[],1));
    X = [X img_vec];
end

X = X-repmat(mean(X,1),243*320,1);
[U,S,V] = svd(X, 'econ');

%% Show Eigenfaces
for i=1:5
    pause(0.1);
    imagesc(reshape(U(:,i),243,320)); colormap gray;
    set(gca,'xtick',[],'ytick',[])
    saveas(gcf, ['./Figures/uncropped' num2str(i) '.png'])
end

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

