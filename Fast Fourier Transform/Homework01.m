%% Clear Workspace
clear all; close all; clc;
%% Load the ultrasound data
load Testdata
%% Setup global variable
AnalysisPlotting = true;
%% Averaging the data in frequency domain to find center frequency.
L=15;
n=64;                                 % fourier modes

x2 = linspace(-L, L, n+1); x=x2(1:n); y=x; z=x;
[X,Y,Z]=meshgrid(x,y,z);

k=(2*pi)/(2*L)*[0:(n/2-1) -n/2:-1];   % scaled wavenumber
ks=fftshift(k);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

K=max(abs(k));

Uave=zeros(n,n,n);
for t=1:20
    Un(:,:,:)=reshape(Undata(t,:),n,n,n);
    Utn=fftn(Un);
    Uave=Uave+Utn;
end
Uave=fftshift(Uave)/20;               % shifted average spectrum

[M,I]=max(Uave(:));
[j,i,k]=ind2sub(size(Uave), I);
kc=[ks(i), ks(j), ks(k)];             % center frequencies
%% 3D Gaussian Filter to filter the data around center frequency.
% filt(:,1)=exp(-(ks-kc(1)).^2);
% filt(:,2)=exp(-(ks-kc(2)).^2);
% filt(:,3)=exp(-(ks-kc(3)).^2);
% 
% gauss_filter=meshgrid(filt(:,1),filt(:,2),filt(:,3));
gauss_filter= exp(-((Kx-kc(1)).^2 + (Ky-kc(2)).^2 + (Kz-kc(3)).^2));

%% Find the path of the marble
path=[];
for t=1:20
    Un(:,:,:)=reshape(Undata(t,:), n,n,n);
    Utn=fftn(Un);                     % fourier transform
    Utf=gauss_filter.*fftshift(Utn);  % filter around center freq
    Uf=ifftn(fftshift(Utf));          % inverse fourier transform
    
    [M,I]=max(abs(Uf(:)));
    [j,i,k]=ind2sub(size(Uf), I);
    location=[x(i), y(j), z(k)];
    path=[path; location];
end
%% Plot the path of the marble
figure(1)
plot3(path(:,1), path(:,2), path(:,3), '-o', 'Linewidth', [3])
grid on, xlabel("x"), ylabel("y"), zlabel("z")
text(path(20,1),path(20,2),path(20,3),...
    "  \leftarrow Location of marble at 20th instance.",...
    'HorizontalAlignment','left','FontSize',10);
fprintf('The location of the marble at 20th instance is %f, %f, %f \n', path(20,:))

%% Additional plot for analysis
if AnalysisPlotting
    points=[[0 0 64]; [64 0 64]; [64 0 0]; [0 0 0]];
    for jj=1:64
        offset(:,:,jj)=[zeros(1,4);jj*ones(1,4);zeros(1,4)]';
    end

    points = bsxfun(@plus, points, offset);
    size(points)

    Xr = reshape(points(:,1,:),4,64);
    Yr = reshape(points(:,2,:),4,64);
    Zr = reshape(points(:,3,:),4,64);

    % Plot the representation of data
    figure(2); grid on;
    fill3(Xr, Yr, Zr, 'r');
    hold on;
    [Xs,Ys,Zs] = meshgrid(0:0.1:20,0:0.1:20,0:0.1:20);
    isosurface(Xs, Ys, Zs, (Xs-10).^2+(Ys-10).^2+(Zs-10).^2, 5)
    alpha(0.3);                         % transparency

    % Plot the isosurface of averaged spectrum
    figure(3);
    set(0, 'defaultTextFontSize',15);
    isosurface(Kx,Ky,Kz,abs(Uave)/max(abs(Uave(:))),0.6)
    axis([-K K -K K -K K]), grid on, drawnow
    xlabel("Wave number (k)")
    ylabel("Wave number (k)")
    zlabel("Wave number (k)")
    
    % Plot the isosurface for 4 stages of tranformation
    figure(4);
    set(0, 'defaultTextFontSize',15);
    hold on;
    
    Un(:,:,:)=reshape(Undata(10,:), n,n,n);
    subplot(2,2,1); view(3); camlight; lighting gouraud
    isosurface(X,Y,Z,abs(Un)/max(abs(Un(:))),0.5)
    axis([-20 20 -20 20 -20 20]), grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title("Spatial raw data")

    Utn=fftn(Un);
    Utn=fftshift(Utn);
    subplot(2,2,2); view(3); camlight; lighting gouraud
    isosurface(Kx,Ky,Kz,abs(Utn)/max(abs(Utn(:))),0.5)
    axis([-K K -K K -K K]), grid on
    xlabel("Wave number (k)")
    ylabel("Wave number (k)")
    xlabel("Wave number (k)")
    title("FFT: Frequency spectrum")

    Utf=gauss_filter.*Utn; 
    subplot(2,2,4); view(3); camlight; lighting gouraud
    isosurface(Kx,Ky,Kz,abs(Utf)/max(abs(Utf(:))),0.3)
    axis([-K K -K K -K K]), grid on
    xlabel("Wave number (k)")
    ylabel("Wave number (k)")
    xlabel("Wave number (k)")
    title("Filtered Frequency spectrum")

    Uf=ifftn(fftshift(Utf));
    subplot(2,2,3);view(3); camlight; lighting gouraud
    isosurface(X,Y,Z,abs(Uf)/max(abs(Uf(:))),0.6)
    axis([-20 20 -20 20 -20 20]), grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title("Denoised spatial data")
end


