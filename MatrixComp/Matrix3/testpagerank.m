%TESTPAGERANK - Tests pagerank as asked for in "Reporting"
%
%   MINIMAL WORKING EXAMPLE: Perform the test
%
%   >> testpagerank

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-10-12: Initial version .
%
% Function code starts here...

%find number of pages with no outgoing links
nonOut = find(sum(HT,1) == 0);
num_nonOut = length(nonOut);

%save urls for pages with no outgoing links
nonOut_urls = cell(num_nonOut,1);
for i = 1:num_nonOut
    nonOut_urls{i,1} = urls{nonOut(i),1};
end

%find which url that has most num of outgoing links and how many they are
[max_out,i_max] = max(sum(HT,1));
max_out_url = urls{i_max,1};
%find which url that has most num of incoming links and how many they are
[max_in,i_maxIn] = max(sum(HT,2));
max_in_url = urls{i_maxIn,1};

alpha = 0.1:0.1:0.9; %define vector for all alphas we want to test
n = length(alpha); %save number number of test-alphas
eps = 1e-6; %define epsilon
%construct empty cell-array of length n for storing ranks
ranks = cell(n,1);
%pre-define zero-matrix for storing number of iterations for each alpha
iter = zeros(n,1);
%scale HT using scalemat enable calling of pagerank(HT,...)
HT_scale = scalemat(HT);

%iterate over each alpha
for i = 1:length(alpha)
    %perform pagerank on HT for each alpha
    [ranks{i}, iter(i)] = pagerank(HT_scale,alpha(i),eps);
end

%find top-ten pages with highest score for alpha = 0.65
rank_top = pagerank(HT_scale,0.65,eps);
[tTen_val,tTen_ind] = maxk(rank_top,10);

%save urls for top-ten-score pages
tTen_urls = cell(10,1);
for i = 1:10
    tTen_urls{i,1} = urls{tTen_ind(i),1};
end


