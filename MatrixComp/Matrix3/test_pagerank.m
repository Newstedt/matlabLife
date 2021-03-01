%TEST_PAGERANK - Tests pagerank on example from book and randomised ST
%
%   MINIMAL WORKING EXAMPLE: Perform the test
%
%   >> test_pagerank

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-10-12: Initial version .
%
% Function code starts here...
P = [0 1/2 1/2 0 0 0; 0 0 0 0 0 0; 1/3 1/3 0 0 1/3 0;...
     0 0 0 0 1/2 1/2; 0 0 0 1/2 0 1/2; 0 0 0 1 0 0];
 
ranks = pagerank(P',0.9,1e-8)

rand('seed',28);
n=1e6;
ST=sprand(n,n,1e-8);
tic,pagerank(ST,0.1,1e-8);toc