function img = crop(image, X, Y, a)
% Crop image using a square window of side 2a at (X,Y).
[M,N,~] = size(image);
% Check if the window is inside the image
inside = (X < (N-a) && X > (a)) && ...
    (Y < (M-a) && Y > (a));
if inside
    img = image(Y-a:Y+a, X-a:X+a,:);
else
    error('Window lies outside the image.');
end