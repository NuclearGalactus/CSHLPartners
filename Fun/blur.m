img = im2double(imread('dog.jpg'));
W = size(img,2);
H = size(img,1);
radius = 1;
width = (2 * radius) + 1;
kernel = zeros(width,width);
kernel(2,1:3) = 1 / 3;
%img2 = zeros(H - (width - 1), W - (width - 1), 3);
%iterate through each RGB channel
for i = 1:3
    img2(:,:,i) = deconvblind(img(:,:,i),kernel);
end
g = fspecial('gaussian',15,2);
imagesc(g)
surfl(g)
gclown = conv2(clown,g,'same');
imagesc(conv2(clown,[-1 1],'same'));
imagesc(conv2(gclown,[-1 1],'same'));
dx = conv2(g,[-1 1],'same');
imagesc(conv2(clown,dx,'same');
lg = fspecial('log',15,2);
lclown = conv2(clown,lg,'same');
imagesc(lclown)
imagesc(clown + .2*lclown)
figure;
subplot(1,2,1);
imshow(img);
subplot(1,2,2);
imshow(img2);