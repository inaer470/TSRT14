function m=exltv(ex,N);

%EXLTI generates standard LTV objects of length N
%   m=exltv(ex,N);
%   Number of samples of the time-varying model is optional (default 500)
%
%   Ex  Type Time   Description
%   -----------------------------
%   mean1  Three changes in the mean (1,2,4) of length 200
%   mean2  As mean1, with softer transitions
%   mean3  As mean1, with random means N(0,10)
%   ar1    An abrupt switch between two AR(2) models
%   ar2    A soft switch between two AR(2) models
%
%   Example:
%      m=exltv('ar2',1000);
%      y=ltv2sig(m);
%      mhat=sig2ltv(y,'ar',[2 0 0 0],'adg',0.95);
%      ltvplot(m,mhat)
%
%   See also: exlti, getsignal, randltv

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<2, N=500; end

switch lower(ex)
    case 'mean1'
        TH=[0 1 3 7];
        j=[1 2 3 4]/4*N;
        nnn=[0 0 0 0];
        m=th2ltv('mean',nnn,TH,j);
    case 'mean2'
        TH=[0 1 3 7];
        j=[1 2 3 4]/4*N;
        nnn=[0 0 0 0];
        m=th2ltv('mean',nnn,TH,j,'ip','on','ipn',2);
    case 'mean3'
        TH=10*randn(1,4);
        j=[1 2 3 4]/4*N;
        nnn=[0 0 0 0];
        m=th2ltv('mean',nnn,TH,j);
    case 'ar1'
        TH=[-1.2 1.2;0.72 0.72];
        j=[1 2]/2*N;
        nnn=[2 0 0 0];
        m=th2ltv('ar',nnn,TH,j);
    case 'ar2'
        TH=[-1.2 1.2;0.72 0.72]';
        j=[1 2]/2*N;
        nnn=[2 0 0 0];
        m=expand(rarx(2),j,TH,[],1,'ip','on','ipn',2,'ipL',40);
end
