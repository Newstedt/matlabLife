%MAIN_OCR - Performs tests relevant in Matrix Computations and
%Applications, Assignment 2, by calling numerous functions and conducting
%a number of plots.
%
%   MINIMAL WORKING EXAMPEL: 
%   ">> main_OCR" to see results in relevant tests.
%
% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-26: Initial version .
%
% Function code starts here...

%load train/test-data
load('mnist_all_converted2.mat');

%change data from uint8 to double
for i=1:10
    OCRtrain{i}=im2double(OCRtrain{i});
    OCRtest{i}=im2double(OCRtest{i});
end

%call function for finding r_k for each digit
rk = find_rk(OCRtrain);
%define r sweep values
r_val = [1 2 4 8 16 32 64 128];
%pre-allocate cell array for storing the accuracy
acc_Cell = cell(length(r_val),1);

%loop over all r
for i = 1:8
    %prepare r for calling train.m : must be in vector form
    r = r_val(i)*ones(length(OCRtrain),1);
    %call function train to compress training data
    [U_r] = train(OCRtrain,r);
    
    %pre-allocate accuracy vector for current r -> will be overwritten
    accuracy = zeros(length(OCRtrain),1);
    %loop over all digits
    for num = 0:9
        %reset number of correct classifications
        correct = 0;
        %loop over all images of a digit
        for j = 1:length(OCRtest{1,num+1}(1,:))
            %call classify.m to classify an image of a digit
            digit = classify(OCRtest{1,num+1}(:,j), U_r);
            
            %check if correctly classified
            if digit == num
                %increase correctly classified digits by one
                correct = correct+1;
            end
            
        end
        %compute accuracy for the current digit of the current r
        accuracy(num+1) = correct/length(OCRtest{1,num+1}(1,:));
    end
    %store vector of accuracies for current r
    acc_Cell{i} = accuracy;
end

%calculate average accuracy for each r in r_val over the digits
for i = 1:8
   mean_r(i) = mean(acc_Cell{i,1});
end

%%
%Calculate average accuracy for each digit over r_val
for i = 1:10
    for j = 1:8
        digit(j) = acc_Cell{j,1}(i);
    end
    mean_dig(i) = mean(digit);
end

%visualise the r-dependence of the mean accuracy
figure
hold on
plot(r_val, mean_r)
scatter(r_val, mean_r, 50, 'filled', 'r')
xlabel('r')
ylabel('avg accuracy')

%visualise which digits is easiest/hardest to classify
figure
hold on
plot(0:9,mean_dig)
scatter(0:9,mean_dig,50,'filled','r')
axis([-0.5 9.5 0.8 1])
xlabel('digit')
ylabel('mean accuracy')

%% This part is to look closer at the digit 5 (hardest to classify)
%set r = 32 since this gives highest accuracy
r = 32*ones(length(OCRtrain),1);
%"train" model
[U_r] = train(OCRtrain,r);

%choose digit 5
num = 5;
%reset number of correctly classified digits
correct = 0;
i = 1;
%loop over all images of 5 in the training set
for j = 1:length(OCRtest{1,num+1}(1,:))
    %classify current image
    digit = classify(OCRtest{1,num+1}(:,j), U_r);
    %check if correctly classified
    if digit == num
        %%increase correctly classified
        correct = correct+1;
    else
        %store wrongly classified
        wrong(i) = digit;
        %store index of wrongly classified
        ind_wrong(i) = j;
        i = i + 1;
    end
    
end
%visualise random wrongly classified image 
im=reshape(OCRtest{num+1}(:,...
           ind_wrong(randsample(length(ind_wrong),1))),[28,28]);
imshow(im,[0,1]);

