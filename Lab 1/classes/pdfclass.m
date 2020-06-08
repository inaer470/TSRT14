classdef pdfclass < handle %%% handle or not!? /ps
%PDFCLASS
%   All PDF objects inherit the PDF class
%
%   help pdfclass.pdfclass      gives help for the constructor
%   methods pdfclass            lists all methods for this class

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

  properties % (Static)
    clslist = {};
  end

  methods (Access = public)

function P = pdfclass
    % No constructor for this class


end

    function l=list(X)
      %LIST returns a list of available PDF objects
      %   l=list(pdfclass)
      if isempty(X.clslist) %%%
        files = dir([fileparts(which('pdfclass')) '/*.m']);
        for n={files.name};
          try
            cls = strrep(n{1}, '.m', '');
            obj = eval(cls);
            if isa(obj, 'pdfclass') & strcmp(cls, 'pdfclass')==0
              X.clslist{end+1} = cls;
            end
          catch ME %%%
          end
        end
      end
      l = X.clslist; %%%
    end %%%

    function Y=horzcat(varargin)
      %HORZCAT
      error('PDFCLASS: Stochastic vectors are column vectors by convention. Use vertcat, [X1;X2], instead')
    end %%%

    function Y=vertcat(varargin)
      %VERTCAT [X1;X2;...] is used to create multivariate statistical variables
      MC=100; %Minimum
      X1=varargin{1};
      if nargin==1
        Y=X1;
      elseif nargin>1
        X2=varargin{2};
        if isa(X1,'pdfclass')
          MC=max([MC,X1.MC]);
        end
        if isa(X2,'pdfclass')
          MC=max([MC,X2.MC]);
        end
        if isa(X1,'pdfclass')
          x1=rand(X1,MC,1);
        else
          x1=X1*ones(MC,1);
        end
        if isa(X2,'pdfclass')
          x2=rand(X2,MC,1);
        else
          x2=X2*ones(MC,1);
        end
        y=[x1 x2];
        Y=empdist(y);
      end
      if nargin>2;  % Recursion
        Y=vertcat(Y,varargin{3:end});
      end
    end

    function str=char(X,format)
    %char, see symbolic
       str=symbolic(X,format)
    end

    %==================
    %----Stats---------
    %==================

    function x=rand(X,n,m)
      %RAND generates random numbers in an (n,m) matrix
      %   x=rand(X,n,m)  for stochastic variables X
      %   x=rand(X,n)    for stochastic vectors X, where x is an (n,length(X)) matrix
      if nargin<3; m=1; end
      if nargin<2; n=1; end
      [P,x]=cdf(X);
      rand('state',X.state);
      u=rand(1,n*m);
      tmp=interp(sig(x,P),u,'method','hold','degree',1);
      x=tmp.y;
      x=reshape(x,n,m);
    end

    function I=erf(X,x)
      %ERF evaluates the error function I(x)=P(X<x) numerically
      %   I=erf(X,x)
      %   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz
      if nargin<2
        error('PDFCLASS.ERF: Syntax erf(X,x)')
      end
      [P,xg]=cdf(X);
      tmp=interp(sig(P,xg),x,'degree',1);
      I=tmp.y;
    end %%%

    function x=erfinv(X,I)
      %ERFINV evaluates the inverse error function I(x)=P(X<x) numerically
      %   x=erfinv(X,I)
      %   The function generalizes the standard erfinv function to non-Gaussian distributions.
      %   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz
      %   That is, I(x)=P(X<x)
      if nargin<2
        error('PDFCLASS.ERFINV: Syntax erfinv(X,I)')
      end
      if I>1 | I<0
        error('PDFCLASS.ERFINV: I must be between 0 and 1')
      end
      [P,xg]=cdf(X);
      tmp=interp(sig(xg,P),I,'degree',1);
      x=tmp.y;
    end %%%

    function [P,x]=cdf(X,x)
      %CDF is the cumulative density function
      %   P=cdf(X,x) for given (vector) x
      %   [P,x]=cdf(X) automatically chooses an interval for x
      %   p is (length(x),1)

      if nargin<2, x=[]; end
      [p,xsort]=pdf(X,x);   % ,x added
      P=cumsum(p.*diff([xsort(1);xsort(:)])); % Riemann approximation of CDF
      P=P-P(1);
      P=P/P(end);
      if isempty(x)
        x=xsort;
      else
        %   if any(x<xsort(1)) | any(x>xsort(end))
        %      error('PDFCLASS.CDF: Specified x value outside observation range in CDF')
        %   end
        ind1=find(x<xsort(1));
        ind2=find(x>=xsort(1) & x<xsort(end));
        ind3=find(x>=xsort(end));
        Porg=P;
        P=zeros(size(x));
        P(ind2)=interp(Porg,xsort,x(ind2),'degree',1);
        P(ind3)=ones(size(ind3));
      end
    end

      function [I,Istd]=ia(X)
        %IA computes the intrisic accuracy of a distribution
        %   I=ia(X)
        %   IA is defined as E[grad(log(p(x|mu)))grad(log(p(x|mu)))'],
        %   where the gradient is with respect to the mean mu
        %'
        x=rand(X,X.MC);
        epsi=std(X)*1e-6;
        p2=pdf(X,x+epsi);
        p1=pdf(X,x);
        g=(log(p2)-log(p1))/epsi;
        I=mean(g.^2);
        Istd=[];
    end

    function [psi,psistd]=ra(X)
        %RA computes the relative accuracy of a distribution
        %   r=ra(X)
        %   RA is defined as IA(x)*var(X) for scalar distributions.
        %   Application: RA says how many times smaller the CRLB is
        %   using X compared with a Gaussian noise with the same variance.
        [I,Istd]=ia(X);
        psi=I*cov(X);
        psistd=[];
    end

    %==================
    %----Operators-----
    %==================

    % Default operators with one argument

    function Y=uminus(X)
      %uminus       - Compute the unary negation of a matrix.
      Y=evalfun('uminus',X);
      static function Y=uplus(X)
      %uplus        - Compute the unary plus of a matrix.
      Y=evalfun('uplus',X);
    end %%%

    % Default operators with two arguments

    %%% static function Y=ldivide(X1,X2)
    function Y=ldivide(X1,X2)
      %ldivide    Divide matrices pointwise.
      Y=evalfun2('ldivide',X1,X2);
    end %%%

    %%% static function Y=minus(X1,X2)
    function Y=minus(X1,X2)
      %minus      Subtract matrices pointwise.
      Y=evalfun2('minus',X1,X2);
    end %%%


    %%% static function Y=mldivide(A,X)
    function Y=mldivide(A,X)  %%% OBS! kommer inte funkar!!!! /ps
      %MLDIVIDE, solves equation system AY=X
      if ~isa(A,'double')  % Y=A*X
        error('First argument must be numeric to MLDIVIDE')
      end
      if size(A,1)~=length(X)
        error('Incorrect size of A in AY=X (size(A,1)=length(X) required)')
      end
      MC=max([100,X.MC]);
      x=rand(X,MC);
      for k=1:MC
        y(k,:)=(A\x(k,:).').';
      end
      Y=empdist(y);
    end

    function Y=mpower(X1,X2)
      %mpower     Matrix power.
      Y=evalfun2('mpower',X1,X2);
    end

    function Y=mrdivide(X1,X2)
      %mrdivide   Solve a linear system of equations.
      Y=evalfun2('mrdivide',X1,X2);
    end

    function Y=mtimes(X1,X2)
      %mtimes       - Compute a matrix product.
      Y=evalfun2('times',X1,X2);
    end

    function Y=plus(X1,X2)
      %plus       Add matrices pointwise.
      Y=evalfun2('plus',X1,X2);
    end

    function Y=power(X1,X2)
      %power      Compute a matrix power pointwise.
      Y=evalfun2('power',X1,X2);
    end

    function Y=rdivide(X1,X2)
      %rdivide    Divide matrices pointwise.
      Y=evalfun2('rdivide',X1,X2);
    end

    function Y=times(X1,X2)
      %times      Multiply matrices pointwise.
      Y=evalfun2('times',X1,X2);
    end



    function Y=fusion(X1,X2)
      %FUSION combines two independent unbiased estimates in an optimal way
      %   X=fusion(X1,X2)
      %
      %   Fusion implements the sensor fusion formula
      %   P = inv( pinv(P1) + pinv(P2) )
      %   x = P * (pinv(P1) * x1 + pinv(P2) * x2)
      %
      %   Example:
      %    X1=ndist([0;0],[2 1;1 1]);
      %    X2=ndist([0;0],[1 1;1 2]);
      %    X=fusion(X1,X2)
      %    plot2(X1,X2,X)

      if ~ (isa(X1,'ndist') & isa(X2,'ndist'))
        error('pdfclass.fusion only implemented for the normal distribution')
      end
      if length(X1)~=length(X2)
        error('pdfclass.fusion: Length of X1 and X2 must be the same')
      end
      I1=pinv(X1.P);
      I2=pinv(X2.P);
      P=pinv(I1+I2);
      mu=P*(I1*X1.mu+I2*X2.mu);
      Y=ndist(mu,P);
    end

    function Y=safefusion(X1,X2)
      %SAFEFUSION combines two dependent unbiased estimates in a worst case approach
      %   X=safefusion(X1,X2)
      %
      %   Safe sensor fusion avoids double-counting information
      %
      %   Example:
      %    X1=ndist([0;0],[2 1;1 1]);
      %    X2=ndist([0;0],[1 1;1 2]);
      %    X=safefusion(X1,X2)
      %    plot2(X1,X2,X)

      if ~ (isa(X1,'ndist') & isa(X2,'ndist'))
        error('pdfclass.safefusion only implemented for the normal distribution')
      end
      if length(X1)~=length(X2)
        error('pdfclass.safefusion: Length of X1 and X2 must be the same')
      end
      I1=pinv(X1.P);
      I2=pinv(X2.P);
      % Step 1
      [U1,D1]=svd(I1);
      % Step 2
      Dhalf=inv(sqrtm(D1));
      [U2,D2]=svd(Dhalf*U1'*I2*U1*Dhalf');
      % Step 3
      T=U2'*sqrtm(D1)*U1';
      %T=U2'*Dhalf*U1'
      % Step 4
      mub1=T*X1.mu;
      mub2=T*X2.mu;
      n=length(X1);
      mub=zeros(n,1);
      Ib=zeros(n);
      % Step 5
      for i=1:n;
        if D2(i,i)>1
          mub(i)=mub2(i);
          Ib(i,i)=D2(i,i);
        else
          mub(i)=mub1(i);
          Ib(i,i)=1;
        end
      end
      % Step 6
      Ti=pinv(T);
      mu=Ti*mub;
      %I=Ti*Ib*Ti';
      I=T'*Ib*T;
      Y=ndist(mu,pinv(I));
    end

    %==================
    %----Elementary----
    %==================

    % Default elementary functions with one argument

    % Default functions with one argument
    function Y=abs(X)
      %abs        Absolute value.
      Y=evalfun('abs',X);
    end

    function Y=acos(X)
      %acos       Inverse cosine.
      Y=evalfun('acos',X);
    end

    function Y=acosh(X)
      %acosh      Inverse cosine
      Y=evalfun('acosh',X);
    end

    function Y=acot(X)
      %acot       Inverse cotangent.
      Y=evalfun('acot',X);
    end

    function Y=acoth(X)
      %acoth      Inverse hyperbolic cotangent.
      Y=evalfun('acoth',X);
    end

    function Y=acsc(X)
      %acsc       Inverse cosecant.
      Y=evalfun('acsc',X);
    end

    function Y=acsch(X)
      %acsch      Inverse hyperbolic cosecant.
      Y=evalfun('acsch',X);
    end

    function Y=angle(X)
      %angle      Polar angle.
      Y=evalfun('angle',X);
    end

    function Y=asec(X)
      %asec       Inverse secant.
      Y=evalfun('asec',X);
    end

    function Y=asech(X)
      %asech      Inverse hyperbolic secant.
      Y=evalfun('asech',X);
    end

    function Y=asin(X)
      %asin       Inverse sine.
      Y=evalfun('asin',X);
    end

    function Y=asinh(X)
      %asinh      Inverse hyperbolic sine.
      Y=evalfun('asinh',X);
    end

    function Y=atan(X)
      %atan       Inverse tangent.
      Y=evalfun('atan',X);
    end

    function Y=atanh(X)
      %atanh      Inverse hyperbolic tangent.
      Y=evalfun('atanh',X);
    end

    function Y=ceil(X)
      %ceil       Round to nearest larger integer.
      Y=evalfun('ceil',X);
    end

    function Y=complex(X)
      %complex    Create a complex matrix.
      Y=evalfun('complex',X);
    end

    function Y=conj(X)
      %conj       Complex conjugate.
      Y=evalfun('conj',X);
    end

    function Y=cos(X)
      %cos        Cosine.
      Y=evalfun('cos',X);
    end

    function Y=cosh(X)
      %cosh       Hyperbolic cosine.
      Y=evalfun('cosh',X);
    end

    function Y=cot(X)
      %cot        Cotangent.
      Y=evalfun('cot',X);
    end

    function Y=coth(X)
      %coth       Hyperbolic cotangent.
      Y=evalfun('coth',X);
    end

    function Y=csc(X)
      %csc        Cosecant.
      Y=evalfun('csc',X);
    end

    function Y=csch(X)
      %csch       Hyperbolic cosecant.
      Y=evalfun('csch',X);
    end

    function Y=exp(X)
      %exp        Exponential.
      Y=evalfun('exp',X);
    end

    function Y=fix(X)
      %fix        Round to nearest integer towards zero.
      Y=evalfun('fix',X);
    end

    function Y=floor(X)
      %floor      Round to nearest smaller integer.
      Y=evalfun('floor',X);
    end

    function Y=imag(X)
      %imag       Imaginary part.
      Y=evalfun('imag',X);
    end

    function Y=isreal(X)
      %isreal     Test if a value is a real matrix.
      Y=evalfun('isreal',X);
    end

    function Y=log(X)
      %log        Natural logarithm.
      Y=evalfun('log',X);
    end

    function Y=log10(X)
      %log10      Base-10 logarithm.
      Y=evalfun('log10',X);
    end

    function Y=log2(X)
      %log2       Base-2 logarithm.
      Y=evalfun('log2',X);
    end

    function Y=mod(X)
      %mod        Compute the modulus of matrices.
      Y=evalfun('mod',X);
    end

    function Y=pow2(X)
      %pow2       Compute or multiply by power of 2.
      Y=evalfun('pow2',X);
    end

    function Y=real(X)
      %real       Real part.
      Y=evalfun('real',X);
    end

    function Y=reallog(X)
      %reallog    Natural logarithm of nonnegative real number.
      Y=evalfun('reallog',X);
    end

    function Y=realpow(X)
      %realpow    Compute power of real matrix.
      Y=evalfun('realpow',X);
    end

    function Y=realsqrt(X)
      %realsqrt   Square root of nonnegative real number.
      Y=evalfun('realsqrt',X);
    end

    function Y=rem(X)
      %rem        Remainder.
      Y=evalfun('rem',X);
    end

    function Y=round(X)
      %round      Round to nearest integer.
      Y=evalfun('round',X);
    end

    function Y=sec(X)
      %sec        Secant.
      Y=evalfun('sec',X);
    end

    function Y=sech(X)
      %sech       Hyperbolic secant.
      Y=evalfun('sech',X);
    end

    function Y=sign(X)
      %sign       Sign of argument.
      Y=evalfun('sign',X);
    end

    function Y=sin(X)
      %sin        Sine.
      Y=evalfun('sin',X);
    end

    function Y=sinh(X)
      %sinh       Hyperbolic sine.
      Y=evalfun('sinh',X);
    end

    function Y=sqrt(X)
      %sqrt       Square root.
      Y=evalfun('sqrt',X);
    end

    function Y=tan(X)
      %tan        Tangent.
      Y=evalfun('tan',X);
    end

    function Y=tanh(X)
      %tanh       Hyperbolic tangent.
      Y=evalfun('tanh',X);
    end

    function Y=unwrap(X)
      %unwrap     Remove phase jumps.
      Y=evalfun('unwrap',X);
    end

    % Default operators with two arguments
    function Y=atan2(X1,X2)
      %atan2      Binary atan.
      Y=evalfun2('atan2',X1,X2);
    end

    %==================
    %----Plots---------
    %==================

    function plot(varargin)
      %PLOT illustrates the PDF of X
      %   plot calls pdfplot
      if length(varargin{1})>1
        error('PDFCLASS.PLOT: Specify the dimension of the stochastic vector for plot by using the call plot(X(dim))')
      end
      pdfplot(varargin{:},'view','pdf')
    end

    function cdfplot(varargin)
      %CDFPLOT illustrates the CDF of X
      %   plot calls pdfplot
      if length(varargin{1})>1
        error('PDFCLASS.CDFPLOT: Specify the dimension of the stochast vector for plot by using the call plot(X(dim))')
      end
      pdfplot(varargin{:},'view','cdf')
    end

    function erfplot(varargin)
      %ERFPLOT illustrates the error function ERF of X
      %   plot calls pdfplot
      if length(varargin{1})>1
        error('PDFCLASS.ERFPLOT: Specify the dimension of the stochast vector for plot by using the call plot(X(dim))')
      end
      pdfplot(varargin{:},'view','erf')
    end

    function varargout=plot2(varargin)
      %PLOT2 illustrates the PDF of X in a two dimensional plot
      %   h=plot2(X1,X2,...,dim,Property1,Value1,...)
      %   Several distributions and/or data sets can be plotted at the same time
      %
      %   X1,X2,...  PDF object
      %   dim is a vector of two integers [i j].
      %   Default dim=[1 2] and the first two dimensions are plotted
      %
      %
      %   Property   Value        Description
      %   ---------------------------------------------------------
      %   conf      {[90]}        Confidence levels in % (scalar or vector)
      %   axis      {gca}         Axis handle where plot is added
      %   col       {'bgrmyk'}    Colors in order of appearance
      %   fontsize  14            Font size
      %   linewidth 2             Line width
      %   markersize10            Markersize
      %   Xlim      {}            Limits on x axis
      %   Ylim      {}            Limits on y axis
      %   legend    {'auto'}      Automatic legend

      if nargin==0
        error('PDFCLASS.PLOT2: No inputs given')
      end

      N=0; K=0; optvar='';
      k=0;
      dim=[1 2];
      while k<length(varargin)
        k=k+1;
        if isa(varargin{k},'pdfclass')  %dist as an object
          N=N+1;
          M{N}=varargin{k};
        elseif isstr(varargin{k})
          K=K+1;
          optvar{K}=varargin{k};
          K=K+1;
          k=k+1;
          optvar{K}=varargin{k};
        elseif isvector(varargin{k})   % Data vector
          dim=varargin{k};
        else
          error('PDFCLASS.PLOT2: Unknown syntax: pdf or property value pairs expected')
        end
      end
      if length(dim)~=2;
        error('PDFCLASS.PLOT2: Dimension dim must be two-dimensional')
      end
opt=struct('conf',[90],'view','pdf','axis',0,'col','bgrmyk','fontsize',[],'linewidth',[],'Xlim',[],'Ylim',[],'legend','auto','markersize',10);
      opt=optset(opt,optvar);
      if opt.axis==0; opt.axis=gca; end
      oldhold = ishold(opt.axis);
      h=[];
      for k=1:N;
        X=M{k};
        leg{k}=symbolic(X);
        if length(X)==1
          error('PDFCLASS.PLOT2 applies only to PDFCLASS objects representing stochastic vectors, not stochastic numbers')
        end
        if max(dim)>length(X) | any(round(dim)~=dim) | min(dim)<0
          error('PDFCLASS.PLOT2: Incorrect dimension in dim')
        end
        if isa(X,'empdist')
          %x=rand(X(dim),X.MC);  %YYYY
          x=rand(X,X.MC);
          kcol=rem(k-1,length(opt.col))+1; %XXX
          htmp=plot(x(:,1),x(:,2),['x' opt.col(kcol)],'parent', ...
                    opt.axis,'linewidth',opt.linewidth,'markersize',opt.markersize);
          hh(k)=htmp(1);
          hold(opt.axis,'on')
        elseif isa(X,'ndist')
          mu=E(X);
          P=cov(X);
          [U,S,V]=svd(P);
          Sroot=sqrt(S);
          Ph=U*Sroot;
          N=100;
          phi=linspace(0,2*pi,N);
          kcol=rem(k-1,length(opt.col))+1; %XXX
	  hh(k)=plot(mu(1),mu(2),['+' opt.col(kcol)],'parent',opt.axis);
          hold(opt.axis,'on')
          for c=opt.conf
	    l=erfinv(expdist(1),c/100);
            plot(mu(1)+l*Ph(1,:)*[cos(phi);sin(phi)],...
              mu(2)+l*Ph(2,:)*[cos(phi);sin(phi)],...
              opt.col(kcol),'parent',opt.axis,'linewidth',opt.linewidth);
          end
        else % generic solution
             %error(['PDFCLASS.PLOT2 does not support pdf of class ',class(X)])
          x=rand(X,X.MC);
          p=pdf(X,x);
          %          size(p),size(x)
          contour2(x(1,:),x(2,:),log(p))
        end
        if isa(X.xlabel,'cell') & length(X.xlabel)>=dim(1) & ~isempty(X.xlabel{dim(1)})
          xlabel(X.xlabel{dim(1)})
        else
          xlabel(['X',num2str(dim(1))])
        end
        if isa(X.xlabel,'cell') & length(X.xlabel)>=dim(2) & ~isempty(X.xlabel{dim(2)})
          ylabel(X.xlabel{dim(2)})
        else
          ylabel(['X',num2str(dim(2))])
        end
      end

      if ~isempty(opt.Xlim)
        set(opt.axis,'Xlim',opt.Xlim);
      end
      if ~isempty(opt.Ylim)
        set(opt.axis,'Ylim',opt.Ylim);
      end
      if strcmp(opt.legend,'auto')
          %legend(gca,leg{:}) %XXX
        legend(hh,leg{:})
      elseif strcmp(opt.legend,'off')
      elseif ~isempty(opt.legend)
        legend(opt.legend)
      end
      if oldhold==1
        hold(opt.axis,'on')
      else
        hold(opt.axis,'off')
      end
    if nargout>0
        varargout{1}=h;
    end

    end

    function surf(X,dim,varargin)
      %SURF illustrates the PDF of X in a two dimensionsal surf plot
      %   surf(X,dim,Property1,Value1,...)
      %   Only one distribution per axis can be illustrated
      %   Empirical distributions are smoothed with a Kernel with width w
      %   which is normalized to the empirical principle axles
      %
      %   Only the multi-variate normal distribution is computed analytically
      %
      %   X      PDF object
      %   dim    Dimensions [i,j] of X. Default dim=[1 2]
      %
      %   Property   Value        Description
      %   ---------------------------------------------------------
      %   w         0.1           Kernel size
      %   axis      {gca}         Axis handle where plot is added
      %   col       {'bgrmyk'}    Colors in order of appearance
      %   fontsize  14            Font size
      %   linewidth 2             Line width
      %   Xlim      {}            Limits on x axis and grid on X(dim(1))
      %   Ylim      {}            Limits on y axis and grid on X(dim(2))

      if nargin==0
        error('PDFCLASS.SURF: No inputs given')
      end

      if nargin<2, dim=[1 2]; end
      if length(dim)~=2;
        error('PDFCLASS.SURF: Dimension dim must be two-dimensional')
      end
      if isstr(dim)
        error('PDFCLASS.SURF: Second argument is the dimension [i,j], not a string')
      end

      opt=struct('w',1,'view','pdf','axis',0,'col','bgrmyk','fontsize',[],'linewidth',[],'Xlim',[],'Ylim',[]);
      opt=optset(opt,varargin);
      if opt.axis==0; opt.axis=gca; end
      if length(X)==1 | ~isa(X,'pdfclass')
        error('PDFCLASS.SURF: applies only for PDFCLASS objects for stochastic vectors, not stochastic numbers')
      end
      if max(dim)>length(X) | any(round(dim)~=dim) | min(dim)<0
        error('PDFCLASS.SURF: Incorrect dimension in dim')
      end

      X.MC=200;
      x=rand(X(dim),X.MC);
      P=cov(X(dim));
      S=inv(opt.w*P);
      xmin=min(x);
      xmax=max(x);
      if ~isempty(opt.Xlim)
        xmin(1)=opt.Xlim(1);
        xmax(1)=opt.Xlim(2);
      end
      if ~isempty(opt.Ylim)
        xmin(2)=opt.Ylim(1);
        xmax(2)=opt.Ylim(2);
      end
      N=20;
      x1=linspace(xmin(1),xmax(1),N)';
      x2=linspace(xmin(2),xmax(2),N)';
      p=zeros(N);
      for k=1:X.MC
        for m=1:N
          for n=1:N
            xx=[x1(m) x2(n)];
            p(m,n)=p(m,n)+exp(-0.5*(xx-x(k,:))*S*(xx-x(k,:))');
          end
        end
      end
      surf(x1,x2,p)
      set(opt.axis,'box','on')
      if ~isempty(opt.Xlim)
        set(opt.axis,'Xlim',opt.Xlim);
      end
      if ~isempty(opt.Ylim)
        set(opt.axis,'Ylim',opt.Ylim);
      end
    end


    function pdfplot(varargin)
      %PDFPLOT plots and compares univariate probabiity density functions (pdf)s
      %   pdfplot(dist1,dist2,...,Property1,Value1,...)
      %   Several distributions and/or data sets can be plotted at the same time
      %   If dist1 is a multivariate distribution, use plot(dist1(dim)), etc.
      %
      %   dist      PDF object
      %
      %   Property   Value        Description
      %   ---------------------------------------------------------
      %   view      {'pdf'}       Plot of probability density function
      %             {'cdf'}       Plot of cumulative density function
      %                                 P(X<x)=int(p(x),-Inf,x)
      %             {'erf'}       Plot of error function (inverse of cdf)
      %   axis      {gca}         Axis handle where plot is added
      %   col       {'bgrmyk'}    Colors in order of appearance
      %   s         []            Smoothing parameter for empdist
      %   x         []            Vector of x values
      %   fontsize  14            Font size
      %   linewidth 2             Line width
      %   Xlim      {}            Limits on x axis
      %   Ylim      {}            Limits on y axis
      %   legend    'auto'        Automatic text for legend using the symbolic method
      %
      %   Examples:
      %      X=ndist(0,1)
      %      x=rand(X,100);
      %      plot(X,empdist(x))
      %
      %   See also: pdfclass.plot2, pdfclass.surf

      %   Fredrik Gustafsson 08-May-2006
      %   Copyright (c) 1994-2006 by COMSOL AB
      %   $ Revision: 28-Oct-2019 $

      if nargin==0
        error('PDFCLASS.PLOT: No inputs given')
      end

      N=0; K=0; optvar='';
      k=0;
      while k<length(varargin)
        k=k+1;
        if isa(varargin{k},'pdfclass')  %dist as an object
          N=N+1;
          M{N}=varargin{k};
        elseif isstr(varargin{k})
          K=K+1;
          optvar{K}=varargin{k};
          K=K+1;
          k=k+1;
          optvar{K}=varargin{k};
        else   % Data vector
          error('PDFCLASS.PLOT: Unknown syntax: pdf or property value pairs expected')
        end
      end
      opt=struct('view','pdf','axis',0,'s',[],'x',[],'col','bgrmyk','fontsize',[],'linewidth',[],'Xlim',[],'Ylim',[],'legend','auto');
      opt=optset(opt,optvar);
      if opt.axis==0; opt.axis=gca; end
      for k=1:N;
        dist=M{k};
        if length(dist)>1
          error('pdfclass.plot cannot plot multivariate distributions') %XXX
        end
        leg{k}=symbolic(dist);
        kcol=rem(k-1,length(opt.col))+1; %XXX
        if strcmpi(opt.view,'pdf')
          if isa(dist,'empdist')
            [p,X]=pdf(dist,[],opt.s);
          else
            [p,X]=pdf(dist,opt.x(:));
          end
          h(k)=plot(X,p,opt.col(kcol),'parent',opt.axis);
          title(opt.axis,'Probability Density Function')
          xlabel(opt.axis,'x')
          ylabel(opt.axis,'p(x)')
        elseif strcmpi(opt.view,'cdf')
          [P,X]=cdf(dist);
	   h(k)=plot(X,P,opt.col(kcol),'parent',opt.axis);
          title(opt.axis,'Cumulative Density Function')
          xlabel(opt.axis,'x')
          ylabel(opt.axis,'P(x)')
        elseif strcmpi(opt.view,'erf')
          [P,X]=cdf(dist);
          h(k)=plot(P,X,opt.col(kcol),'parent',opt.axis);
          title(opt.axis,'Error Function')
          xlabel(opt.axis,'P(x)')
          ylabel(opt.axis,'x')
        end
        hold(opt.axis,'on')
      end
      if ~isempty(opt.Xlim)
        set(opt.axis,'Xlim',opt.Xlim);
      end
      if ~isempty(opt.Ylim)
        set(opt.axis,'Ylim',opt.Ylim);
      end
      if strcmp(opt.legend,'auto')
        legend(h,leg{:})
      elseif ~isempty(opt.legend)
        legend(opt.axis,opt.legend)
      end
      hold(opt.axis,'off')

    end

  end
end


function Y=evalfun(fun,X)
%EVALFUN uses feval to evaluate functions of a s.v.
MC=100; %Minimum
MC=max([MC,X.MC]);
x=rand(X,MC,1);
y=feval(fun,x);
Y=empdist(y);
Y.MC=MC;
end

function Y=evalfun2(fun,X1,X2)
%EVALFUN2 uses feval to evaluate functions of two arguments
if isa(X2,'sig')  % Exception case, the SIG method should be used
  %  Y=sig.evalfun2(fun,X2,X1);
  Y=feval(fun,X2,X1);
  return
end
MC=100; %Minimum
if isa(X1,'pdfclass')
  MC=max([MC,X1.MC]);
end
if isa(X2,'pdfclass')
  MC=max([MC,X2.MC]);
end
if isa(X1,'pdfclass')
  x1=rand(X1,MC,1);
else
  x1=X1;
end
if isa(X2,'pdfclass')
  x2=rand(X2,MC,1);
else
  x2=X2;
end
y=feval(fun,x1,x2);
Y=empdist(y);
Y.MC=MC;
end
