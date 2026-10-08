clear all, close all

load vid2.mat
recording = double(recording1);
recording = imresize(recording, [640, 480]);

figure
subplot(1,2,1)
imagesc(recording(:,:,1,1))
colormap('gray')
grid on

subplot(1,2,2)
imagesc(recording(:,:,1,3))
colormap('gray')
grid on

figure
imagesc(abs(recording(:,:,1,1)-recording(:,:,1,3)))
colormap('hsv')
grid on

