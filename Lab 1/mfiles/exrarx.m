function m=exrarx(ex,N);

%EXRARX generates standard RARX objects of length N
%   m=exrarx(ex,N);
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
%     mt=exrarx('rar2',1000);
%     y=simulate(m);
%     mthat=estimate(rarx(2),y,'adg',0.95);
%     plot(mt,mthat,'Ylim',[-2 2])
%
%     mt=exrarx('mean1',1000);
%     u=getsignal('unit',1000);
%     y=simulate(mt,u);
%     plot(y)
%
%   See also: getsignal, exlti, randrarx

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<2, N=500; end

switch lower(ex)
    case 'mean1'
        TH=[0 1 3 7]';
        j=[1 2 3 4]/4*N;
        nnn=[0 1 0];
        m=expand(rarx(nnn),j,TH);
    case 'mean2'
        TH=[0 1 3 7]';
        j=[1 2 3 4]/4*N;
        nnn=[0 1 0];
        m=expand(rarx(nnn),j,TH);
    case 'mean3'
        TH=10*randn(4,1);
        j=[1 2 3 4]/4*N;
        nnn=[0 1 0];
        m=expand(rarx(nnn),j,TH);
    case 'rar1'
        TH=[-1.2 1.2;0.72 0.72]';
        j=[1 2]/2*N;
        nnn=[2 0 0];
        m=expand(rarx(nnn),j,TH);
    case 'rar2'
        TH=[-1.2 1.2;0.72 0.72]';
        j=[1 2]/2*N;
        nnn=[2 0 0];
        m=expand(rarx(nnn),j,TH,[],[],'ip','on','ipn',2,'ipL',40);
    case 'rarx1'
        TH=[-1.2 0.72 1 -0.5;1.2 0.72 1 -0.8];
        j=[1 2]/2*N;
        nnn=[2 2 1];
        m=expand(rarx(nnn),j,TH);
    case 'rarx2'
        TH=[-1.2 0.72 1 -0.5;1.2 0.72 1 -0.8];
        j=[1 2]/2*N;
        nnn=[2 2 1];
        m=expand(rarx(nnn),j,TH,[],[],'ip','on','ipn',2,'ipL',40);
    otherwise
        error(['Unknown string option: ',ex])
end
