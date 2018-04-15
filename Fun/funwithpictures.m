img = im2double(imread('dog.jpg'));
W = size(img,2);
H = size(img,1);
radius = 20;
Ws = 1:radius + 1:W;
Hs = 1:radius + 1:H;
newW = length(Ws);
newH = length(Hs);
img2 = zeros(newH,newW,3);
for w = 1:newW
    for h = 1:newH
        %for i = Ws(w) - 1:Ws(w) + 1
          %  for j = Hs(h) - 1:Hs(h) + 1
         img2(h,w,:) = img(Hs(h),Ws(w),:);
         %   end
        %end
    end
end

figure;
subplot(1,2,1);
imshow(img);
subplot(1,2,2);
imshow(img2);