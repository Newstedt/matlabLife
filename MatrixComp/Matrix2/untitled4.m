clear all; clc;
load('mnist_all_converted2.mat');

for i=1:10
    OCRtrain{i}=im2double(OCRtrain{i});
    OCRtest{i}=im2double(OCRtest{i});
end

k=0;

im=reshape(OCRtrain{k+1}(:,1),[28,28]);
imshow(im,[0,1]);
%%
k=0;
A = OCRtrain{1,k+1};
[U,S,V] = svd(A,'econ');
Ak = U*S*V';
norm_Ak = norm(Ak);
U = U(:,1:9);
S = S(1:9,1:9);
V = V(:,1:9);
comp = norm(Ak-U*S*V')/norm_Ak


figure
semilogy(S,'.')
xlabel('Basis vectors')
ylabel('Singular values  (\sigma_i)')
title('Singular value-plot for digit: 1')
%%
r_k = zeros(length(OCRtrain),1);
for j = 1:length(OCRtrain)
    A = OCRtrain{1,j};
    [U,S,V] = svd(A,'econ');
    Ak = U*S*V';
    norm_Ak = norm(Ak);
    c = 1;
    i = 20;
    while c == 1 && i >= 1
        U = U(:,1:i);
        S = S(1:i,1:i);
        V = V(:,1:i);
        comp = norm(Ak-U*S*V')/norm_Ak;
        
        if comp >= 0.12
            r_k(j) = i+1;
            c = 0;
        end
        i = i-1;
    end
end

rk = [10 9 15 12 15 16 12 12 13 14];


