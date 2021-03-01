load('mnist_all_converted2.mat');

for i=1:10
    OCRtrain{i}=im2double(OCRtrain{i});
    OCRtest{i}=im2double(OCRtest{i});
end

rk = 32*ones(length(OCRtrain),1);
[U_r] = train(OCRtrain,rk);

accuracy = zeros(length(OCRtrain),1);
num = 5;
correct = 0;
i = 1;
for j = 1:length(OCRtest{1,num+1}(1,:))
    digit = classify(OCRtest{1,num+1}(:,j), U_r);
    
    if digit == num
        correct = correct+1;
    else
        wrong(i) = digit;
        ind_wrong(i) = j;
        i = i + 1;
    end
    
end
%%
im=reshape(OCRtest{num+1}(:,ind_wrong(39)),[28,28]);
imshow(im,[0,1]);

% 1 2 4 7 8 9 11 17 19 22 25 26 29 32 33 34 44
% 5    12    13    18    20    23    27    28    30    36    38    39    40    42    46    47    51


