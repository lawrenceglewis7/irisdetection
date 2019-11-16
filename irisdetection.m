
clear all
close all 
clc

eye = imread('C:\Users\Lawrence\Documents\MATLAB\imageProcessing\eyes\eye2.JPG'); %unfiltered image

greyEye = rgb2gray(eye);
[rows, cols] = size(greyEye);
imshow(greyEye);

%%
%gausian filter
lowpass = GausFilter(greyEye,rows/32, 0,1);

minlowpas = double(min(lowpass(:)));
maxlowpas = double(max(lowpass(:)));

lowpass = lowpass * (255/maxlowpas);




%%
%derivatives
dxKernel = [1 1 1;0 0 0; -1 -1 -1];
dyKernel = [1 0 -1; 1 0 -1; 1 0 -1];

partialx = conv2(lowpass, dxKernel,'same');
partialy = conv2(lowpass, dyKernel,'same');

partialxx = conv2(partialx, dxKernel,'same');
partialyy = conv2(partialy, dyKernel,'same');

partialx = uint8(partialx);
partialy = uint8(partialy);

partialxx = uint8(partialxx);
partialyy = uint8(partialyy);

%%
% try to binerize the image. 
%partialxx = imbinarize(partialxx,'adaptive','Sensitivity',0);
%partialyy = imbinarize(partialyy,'adaptive','Sensitivity',.000000001);

addxxandyy = partialxx + partialyy;

binaryimage = imbinarize(addxxandyy,'adaptive','Sensitivity',.4);

%%
%pickout each and every individual blob
allgroups = DifBlobs(binaryimage);

%%
minpartialxx = double(min(partialxx(:)));
maxpartialxx = double(max(partialxx(:)));
figure
imshow(partialxx,[minpartialxx,maxpartialxx])

minpartialyy = double(min(partialyy(:)));
maxpartialyy = double(max(partialyy(:)));
figure
imshow(partialyy,[minpartialyy,maxpartialyy])

figure
histogram(partialxx)

figure
histogram(partialyy)


minLowp = double(min(lowpass(:)));
maxLowp = double(max(lowpass(:)));
figure
imshow(lowpass,[minLowp,maxLowp])

figure
imshow(addxxandyy)

figure
imshow(binaryimage)

figure
imshow(allgroups{1})

figure
imshow(allgroups{2})

figure
imshow(allgroups{3})

figure
imshow(allgroups{4})

figure
imshow(allgroups{5})

figure
imshow(allgroups{6})

figure
imshow(allgroups{7})
figure
imshow(allgroups{8})
figure
imshow(allgroups{9})
figure
imshow(allgroups{10})
figure
imshow(allgroups{11})
figure
imshow(allgroups{12})
figure
imshow(allgroups{13})
%%
%functions
%if loworhigh = 0 lowpass else highpass
function gausianFilter = GausFilter(imgInput,standardDeviation, mean,loworHigh)
    
    [ro, col] = size(imgInput);
    gausianLight = zeros(ro,col);
    centerY = ro/2;
    centerX = col/2;

    for i=1:ro
        for j=1:col
            r = sqrt((i-centerY).^2 +(j-centerX).^2);
            gausianLight(i,j) = 1./(standardDeviation*sqrt(2*pi))*exp(-(1/2)*((r-mean)./standardDeviation).^2);
        end
    end

    invGL = zeros(ro,col);
    for i=1:ro
        for j=1:col
            invGL(i,j) =  1 - gausianLight(i,j);
        end
    end
    
    fftshifted = fftshift(fft2(imgInput));
    lowpassImage  = abs(ifft2(ifftshift(fftshifted  .* gausianLight)));
    highpassImage = abs(ifft2(ifftshift(fftshifted  .* gausianLight)));
    
    if loworHigh == 0
        gausianFilter = lowpassImage;
    else
        gausianFilter = highpassImage;
    end
end
function differentiatbinaryBlobswithintensity = DifBlobs(binaryinputImg)
    [ro, col] = size(binaryinputImg);
    allBlobs{20} = zeros(ro,col,20);
    previousBlob = zeros(ro,col);
    SE = strel('square',3);
    binaryinputImg = logical(binaryinputImg);
    
    blobscounter = 1;
    for i = 1:ro
        for j= 1:col
            if binaryinputImg(i,j) == 1
                currentBlobs = zeros(ro,col);
                currentBlobs (i,j) = 1;
                while(~isequal(currentBlobs, previousBlob))
                    previousBlob = logical(currentBlobs) ;
                    currentBlobs = logical(imdilate(currentBlobs,SE));
                    currentBlobs = and(currentBlobs, binaryinputImg);
%                     figure
%                     imshow(currentBlobs)
%                     title('blobx')
                end
                
                binaryinputImg = binaryinputImg - currentBlobs;
                allBlobs{blobscounter} = currentBlobs;
                blobscounter = 1 + blobscounter;
            end
            
        end
    end
    differentiatbinaryBlobswithintensity = allBlobs;
end