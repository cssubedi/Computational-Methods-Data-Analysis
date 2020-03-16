function Af = Smoothing(image_grayscale)
A = double(image_grayscale);
At = fft2(A);
Ats = fftshift(At);

[r,c] = size(image_grayscale);
kx=1:c; ky=1:r; 
[Kx,Ky]=meshgrid(kx,ky); 
F=exp(-0.0001*(Kx-321).^2-0.0001*(Ky-241).^2);
Atsf=Ats.*F; 
Atf=ifftshift(Atsf); 
Af=ifft2(Atf);

figure(1), 
subplot(1,2,1), pcolor(log(abs(Atsf)))
colormap(gray), shading interp
set(gca,'Xtick',[],'Ytick',[]) 
subplot(1,2,2), imshow(uint8(Af)), colormap(gray)
end
