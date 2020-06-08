function H=histeq(H);

%HISTEQ rescales the values in the matrix H uniformly to the interval [0,1]
%   H=histeq(H);
%   In histogram equalization, a monotonous mapping is applied such that
%   the original values in H are mapped to values in [0,1] such that all
%   values are evenly spread as in a uniform distribution.
%   The function is used in ltvplot and tfdplot for instance for improved visibility.
%
%   Example:
%     m=exltv('ar2',1000);
%     M=ltv2tfd(m);
%     subplot(2,1,1),  tfdplot(M,'view','surf','histeq','on')
%     subplot(2,1,2),  tfdplot(M,'view','surf','histeq','off')
%
%   See also: ltvplot, tfdplot

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

[n,m]=size(H);
[h,ind]=sort(H(:));
h(ind)=(1:length(ind))/length(ind);
H=reshape(h,n,m);
