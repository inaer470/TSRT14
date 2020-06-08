classdef logdist < pdfclass
%LOGDIST defines the logistic distribution log(mu,P)
%
%   p=exp(-(x-mu)/s)/s/(1+exp(-(x-mu)/s))
%   P=1/(1+exp(-(x-mu)/s)) = sigmoid function
%
%   help logdist.logdist      gives help for the constructor
%   methods logdist         lists all methods for this class


%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

  properties (SetAccess = public)
    mu;
    s;
  end

  properties (SetAccess = public)
    MC;
    state;
    xlabel;
  end

  methods
    function X = logdist(mu,s)
      % Constructor
      %   X=logdist(mu,s)  defines the distribution with variance/covariance s>0
      %   X=logdist        defines an exmpty distribution
      %   X=logdist(X)     makes a copy of the distribution with a new state

%#
%#

      X.MC=1000; % Default number
      X.state=round(1e6*rand); %Random state
      if nargin==0
        X.mu=[]; X.s=[];
        X.xlabel='';
      elseif isa(mu,'logdist') % Copy PDF with new state
        X.mu=mu.mu;
        X.s=mu.s;
        X.MC=mu.MC;
        X.xlabel=mu.xlabel;
        X.state=round(1e6*rand); %Random state
      elseif nargin==2
        if size(mu,2)>1 | ~isreal(mu)
          error('mu must be a real column vector')
        end
        n=length(mu);
        for i=1:n
          X.xlabel{i}=['x',num2str(i)]; %%%
        end
        if ~isscalar(s)
          error('LOGDIST: s must be scalar')
        end
        if  ~isreal(s)
          error('LOGDIST: s must be a real')
        end
        if s<=0
          disp('LOGDIST warning: s must be positive')
        end
        X.s = s;   %%%
        X.mu = mu; %%%
      else
        error('Syntax: logdist(mu,s), or just logdist')
      end

    end

    function X=set.mu(X,mu)
    if ~isempty(X.mu) & length(mu)~=length(X.mu)
        error(['Cannot change the dimensions of mu'])
      end
      X.mu=mu;
    end
    function X=set.s(X,s)
      if ~isempty(X.s)
          if length(s)~=length(X.s)
             error(['Cannot change the dimensions of s'])
          elseif ~isscalar(s) | s<=0
              error(['s is not a valid scaling'])
          end
      end
      X.s=s;
    end


    function out=fieldread(X, arg)
      % Reads the field named arg, ex X.arg
      out=eval(X.arg);
    end


    function out=fieldwrite(X, arg1,arg2)
      % Writes the field named arg1, ex X.arg1=arg2
      if ~all(size(eval(X.arg1))==size(arg2))
        error(['Cannot change the dimensions of ',arg1])
      end
      eval([X.arg1,'=arg2;'])
    end

function Xi=Xsubsref(X,S)
%function Xi=subsref(X,S)
%SUBSREF is used to pick out parts of stochastic vector
%   Xi=X(i)
%   i is the row index/indices
      if S.type~='()'
        error('EMPDIST: Use parenthesis for indexing')
      end
      if length(S.subs)>1
         error('EMPDIST: object is one-dimensional')
      end
      i=S.subs{1};
      if any(i<1) | any(i>length(X.mu))
        error('Index out of range')
      end
      Xi=logdist(X.mu(i),X.s(i,i));
      for k=1:length(i)
         Xi.xlabel{k}=X.xlabel{i(k)};
      end
    end

    function disp(X)
      format='%11.3g';
      if isempty(X.mu)
        mustr='mu';
      else
        if length(X.mu)==1
          mustr=num2str(X.mu,format);
        else
          mustr=mat2strformat(X.mu,format);
        end
      end
      if isempty(X.s)
        sstr='s';
      else
        if length(X.mu)==1
          sstr=num2str(X.s,format);
        else
          sstr=mat2strformat(X.s,format);
        end
      end
      disp(['logdist(',mustr,',',sstr,')']);
    end

    function out=desc(X)
      %DESC Returns a description of the distribution defined by the class.
      out='Logistic distribution';
    end


    function str=symbolic(X,format)
      %SYMBOLIC returns a symbolic expression logdist(mu,s)
      %   str=symbolic(X,format)
      %   format is the numeric format used in mat2str
      %   Default format is '%11.3g'
      if nargin<2, format='%11.3g'; end
      if isempty(X.mu)
        mustr='mu';
      else
        if length(X)==1
          mustr=num2str(X.mu,format);
        else
          mustr=mat2strformat(X.mu,format);
        end
      end
      if isempty(X.s)
        sstr='s';
      else
        if length(X)==1
          sstr=num2str(X.s,format);
        else
          sstr=mat2strformat(X.s,format);
        end
      end
      str=['logdist(',mustr,',',sstr,')'];
    end

    function texcode=tex(X,varargin)
      %TEX returns latex code for the LOGDIST object logdist(mu,s)
      %   texcode=tex(X,Property1,Value1,...)
      %
      %   Property   Value/{Default}  Description
      %   --------------------------------------------------------------
      %   filename   {''}             Name of the .tex file (none for '')
      %   decimals   {1}              Number of decimals
      %   env        {'eqnarray*'}    Tex environment, '' means no env

      opt=struct('filename','','decimals',1,'env','eqnarray*');
      opt=optset(opt,varargin);

      if isempty(X.mu)
        mustr='mu';
      else
        mustr=texmatrix(X.mu,'decimals',opt.decimals,'env',opt.env);
      end
      if isempty(X.s)
        sstr='s';
      else
        sstr=texmatrix(X.s,'decimals',opt.decimals,'env',opt.env);
      end
      texcode=char({'N\left(',mustr,',',sstr,'\right)'});

      if ~isempty(opt.filename)
        eval(['fid=fopen(''',opt.filename,'.tex'',''w'')']);
        for j=1:size(texcode,1);
          rowj=texcode(j,:);
          rowj=strrep(rowj,'\','\\');
          fprintf(fid,[rowj,' \n']);
        end
        fclose(fid);
      end
    end

    function n=length(X)
      %LENGTH is the length of the stochastic vector X
      n=size(X.mu,1);
    end



    %==================
    %----Moments-------
    %==================

    function mu=median(X)
      %MEDIAN is the median operator
      mu=X.mu;
    end

    function mu=mode(X)
      %MODE is the mode operator
      mu=X.mu;
    end

    function mu=E(X)
      %E is the expectation operator
      mu=X.mu;
    end

    function mu=mean(X)
      %MEAN is the expectation operator
      mu=X.mu;
    end

    function v=var(X)
      %VAR is the variance operator
      b=pi^2/3*X.s^2;
    end

    function s=std(X)
      %STD is the standard deviation operator
      s=pi/sqrt(3)*X.s;
    end

    function out=skew(X)
      %SKEW is the skewness operator
      if length(X.s)>1
        error('skew(X) not defined for multi-dimensional variables')
      end
      out=0;
    end

    function out=kurt(X)
      %KURT is the kurtosis operator
      if length(X.s)>1
        error('kurt(X) not defined for multi-dimensional variables')
      end
      out=6/5;
    end

    function s=cov(X)
      %COV is the covariance operator
      s=pi^2/3*X.s^2;
    end

    %==================
    %----Estimators----
    %==================

    function [Y,msg]=estimate(X,Z)
      %ESTIMATE computes a moment based estimate of Z in the Gaussian class
      %   X=estimate(logdist,Z)
      m=E(Z);
      v=cov(Z);
      Y=logdist(m,sqrt(v*3/pi^2));
      msg=[];
    end

    %==================
    %----Stats---------
    %==================

    function x=rand(X,n,m)
      %RAND generates random numbers in an (n,m) matrix
      %   x=rand(X,n,m)  for stochastic variables X
      %   x=rand(X,n)    for stochastic vectors X, where x is an (n,length(X)) matrix
      randn('state',X.state);
      if nargin<3,  m=1; end
      if nargin<2,  n=1; end
      if isempty(n) | ~isreal(n) | n<0 | round(n)~=n
        error('logdist.rand: n must be a non-negative integer'),
      end
      if isempty(X.mu) | isempty(X.s)
          error(['logdist.rand: distribution parameters not ' ...
                 'specified'])
      end

      u=rand(n,m);
      x=X.mu-X.s.*log(1./u-1);
    end

    function I=cdf(X,x)
      %CDF evaluates the error function I(x)=P(X<x)
      %   I=cdf(X,x)
      %   The error function is defined as I(x)=int_{-Inf}^{x} p(z) dz
      if nargin<2
        error('Syntax: cdf(X,x)')
      end
      I=1./(1+exp(-(x-X.mu)./X.s));
    end


    function x=cdfinv(X,I)
      %CDFINV evaluates the inverse error function I(x)=P(X<x) numerically
      %   x=cdfinv(X,I)
      %   The function generalizes the standard cdfinv function to non-Gaussian distributions.
      %   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz
      %   That is, I(x)=P(X<x)
      if nargin<2
        error('Syntax: cdfinv(X,I)')
      end
      if length(X)>1
        error('CDF only defined for scalar stochastic variables')
      end
      if I>1 | I<0
        error('I must be between 0 and 1')
      end
      x=X.mu-X.s.*log(1./I-1);
    end


    function [p,x]=pdf(X,x)
      %PDF is the probability density function
      %   p=pdf(X,x) for given x (dim (N,nx))
      %   [p,x]=pdf(X) automatically chooses a (N,nx) grid for x
      %   If X is multivariate, then the pdf is evaluated for each row of x, and
      %   length(p) = size(x,1)
      m=X.mu;
      s=X.s;
      nx=length(X);
      if nargin<2 | isempty(x)
          x=linspace(m-6*s,m+6*s,1000)';
      else
          if (size(x,2)~=length(m)) & length(m)~=1
            error('LOGDIST.PDF: x must have the same number of columns as length(mu)')
          end
      end
      e=exp(-(x-m)./s);
      p=e./s./(1+e).^2;
    end

    %==================
    %----Operators-----
    %==================

    function Y=vertcat(X1,varargin)
      %VERTCAT [X1;X2] is used to create multivariate normal distributions
      if nargin==1
        Y=X1; %%% copy or not??? /ps
        %error('LOGDIST.VERTCAT: at least two input arguments required')
      else
        X2=varargin{1};
        if isa(X2,'logdist')
          Y=logdist([X1.mu;X2.mu],blkdiag(X1.P,X2.P));
        else
          MC=X1.MC;
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
      end
      if nargin>2;  % Recursion
        Y=vertcat(Y,varargin{2:end});
      end
    end

    function Y=uplus(X)
      %UPLUS unitary plus, Y=+X
      Y=logdist(X.mu,X.P); %%% return a copy or not? alt. just return X; /ps
    end

    function Y=uminus(X)
      %UMINUS unitary minus, Y=-X
      Y=logdist(-X.mu,X.P);
    end

    function Y=plus(X1,X2)
      %PLUS, Y=X1+X2
      if isa(X2,'sig')  % Exception for calls like n+y
        Y=X2+X1;
        return
      end
      if isa(X1,'double')  % Y=AX
        mu2=X2.mu;
        P2=X2.P;
        if length(X1)==1
          mu1=X1*ones(length(mu2),1);
        else
          mu1=X1;
        end
        P1=zeros(length(mu1));
        Y=logdist(mu1+mu2,P1+P2);
      elseif isa(X2,'double')  % Y=Xa
        mu1=X1.mu;
        P1=X1.P;
        if length(X2)==1
          mu2=X2*ones(length(mu1),1);
        else
          mu2=X2;
        end
        P2=zeros(length(mu2));
        Y=logdist(mu1+mu2,P1+P2);
      elseif isa(X2,'pdfclass')
        if isa(X2,'gmdist')
	  Y=X2+X1;   % Use the code in gmdist
        elseif isa(X2,'logdist')
          mu1=X1.mu; P1=X1.P;
          mu2=X2.mu; P2=X2.P;
          if size(mu1,1)~=size(mu2,1)
            error('X1 and X2 not compatible for PLUS')
          end
          Y=logdist(mu1+mu2,P1+P2);
        else
          MC=max([100 X1.MC X2.MC]);
          x1=rand(X1,MC,1);
          x2=rand(X2,MC,1);
          Y=empdist(x1+x2);
        end
      else
        error('Incorrect argument to PLUS in LOGDIST class')
      end
    end


    function Y=minus(X1,X2)
      %MINUS, Y=X1-X2
      Y=X1+(-X2);
    end

    function Y=mrdivide(X1,X2,MC)
      %MRDIVIDE, division Y=X1/X2
      if isa(X2,'sig')  % Exception for calls like n.*y
        Y=X2/X1;
      elseif isa(X2,'double')
        Y=X1*(1/X2);
      end
    end


    function Y=rdivide(X1,X2,MC)
      %RDIVIDE, elementwise division Y=X1./X2
      if isa(X2,'sig')  % Exception for calls like n.*y
        Y=X2./X1;
      elseif isa(X2,'double')
        Y=X1.*(1./X2);
      end
    end


    function Y=times(X1,X2,MC)
      %TIMES, elementwise multiplication Y=X1.*X2
      if isa(X2,'sig')  % Exception for calls like n.*y
        Y=X2.*X1;
        return
      end
      if nargin<3, MC=1000; end
      if isa(X1,'double')  % Y=A.*X
        if size(X1,1)~=size(X2.mu,1);
          error('X1 and X2 not compatible for TIMES')
        end
        Y=logdist(X1.*X2.mu,diag(X1)*X2.P*diag(X1));
      elseif isa(X2,'double') %Y=XA
        if size(X1.mu,1)~=size(X2,1);
          error('X1 and X2 not compatible for TIMES')
        end
        Y=logdist(X1.mu.*X2,diag(X2)*X1.P*diag(X2));
      else  % Both Gaussian
        x1=rand(X1,MC);
        x2=rand(X2,MC);
        Y=empdist(x1.*x2);
      end
    end

    function Y=power(X,n,MC)
      %POWER, Y=X.^n
      if nargin<3; MC=1000; end
      Y=mpower(X,n,MC);
    end

    function Y=mpower(X,n,MC)
      %MPOWER, Y=X^n
      if X.mu==0 & X.P==1 & n==2
        Y=chi2dist(1);
      else
        if nargin<3; MC=1000; end
        x=rand(X,MC,1);
        Y=empdist(x.^n);
      end
    end


    function Y=mtimes(X1,X2,MC)
      %MTIMES, matrix multiplication Y=X1*X2
      if nargin<3, MC=1000; end
      if isa(X1,'double')  % Y=A*X
        if ~isscalar(X1) & size(X1,2)~=size(X2.mu,1);
          error('X1 and X2 not compatible for MTIMES')
        end
        Y=ndist(X1*X2.mu,X1*X2.P*X1');
      elseif isa(X2,'double') %Y=XA
        if ~isscalar(X2)
          error('X1 and X2 not compatible for MTIMES')
        end
        Y=logdist(X1.mu.*X2,X2.^2*X1.P);
        %elseif isa(X1,'logdist') & isa(X2,'logdist')   % Both Gaussian
      else  % Both distributions
        if length(X1)>1 | length(X2)>1
          error('X1 and X2 not compatible for MTIMES')
        end
        x1=rand(X1,MC,1);
        if length(X1)==1
          x1=x1*ones(1,length(X2));
        end
        x2=rand(X2,MC,1);
        if length(X2)==1
          x2=x2*ones(1,length(X1));
        end
        Y=empdist(x1.*x2);
      end
    end

  end
end
