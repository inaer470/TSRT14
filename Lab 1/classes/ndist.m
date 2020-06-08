classdef ndist < pdfclass
%NDIST defines the Gaussian, or normal, distribution N(mu,P)
%
%   help ndist.ndist      gives help for the constructor
%   methods ndist         lists all methods for this class


%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $


  properties (SetAccess = public)
    mu;
    P;
  end

  properties (SetAccess = public)
    MC;
    state;
    xlabel;
  end

  methods
    function X = ndist(mu,P)
      % Constructor
      %   X=ndist(mu,P)  defines the distribution with variance/covariance P>0
      %   X=ndist        defines an exmpty distribution
      %   X=ndist(X)     makes a copy of the distribution with a new state

%#
%#

      X.MC=1000; % Default number
      X.state=round(1e6*rand); %Random state
      if nargin==0
        X.mu=[]; X.P=[];
        X.xlabel='';
      elseif isa(mu,'ndist') % Copy PDF with new state
        X.mu=mu.mu;
        X.P=mu.P;
        X.MC=mu.MC;
        X.xlabel=mu.xlabel;
        X.state=round(1e6*rand); %Random state
      elseif isa(mu,'gmdist') %
        X.mu=mean(mu);
        X.P=cov(mu);
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
        if size(P,1)~=n | size(P,2)~=n
          error('NDIST: P must be square with same dimension as mu')
        end
        if  ~isreal(P)
          error('NDIST: P must be a real matrix')
        end
        if norm(P-P')>1e-10*norm(P)
          error('NDIST: P must be symmetric matrix')
        end
        if any(abs(eig(P))<-eps)
          disp('NDIST warning: P must be positive semidefinite matrix')
        end
        X.P = P;   %%%
        X.mu = mu; %%%
      else
        error('Syntax: ndist(mu,P), or just ndist')
      end

    end

    function X=set.mu(X,mu)
    if ~isempty(X.mu) & length(mu)~=length(X.mu)
        error(['Cannot change the dimensions of mu'])
      end
      X.mu=mu;
    end

    function X=set.P(X,P)
      if ~isempty(X.P) & length(P)~=length(X.P)
        error(['Cannot change the dimensions of P'])
      end
      if ~iscov(P)
        error(['P is not a valid covariance'])
      end
      X.P=P;
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
      Xi=ndist(X.mu(i),X.P(i,i));
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
      if isempty(X.P)
        Pstr='P';
      else
        if length(X.mu)==1
          Pstr=num2str(X.P,format);
        else
          Pstr=mat2strformat(X.P,format);
        end
      end
      disp(['N(',mustr,',',Pstr,')']);
    end

    function out=desc(X)
      %DESC Returns a description of the distribution defined by the class.
      out='Normal distribution';
    end


    function str=symbolic(X,format)
      %SYMBOLIC returns a symbolic expression N(mu,P)
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
      if isempty(X.P)
        Pstr='P';
      else
        if length(X)==1
          Pstr=num2str(X.P,format);
        else
          Pstr=mat2strformat(X.P,format);
        end
      end
      str=['N(',mustr,',',Pstr,')'];
    end

    function texcode=tex(X,varargin)
      %TEX returns latex code for the NDIST object N(mu,P)
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
      if isempty(X.P)
        Pstr='P';
      else
        Pstr=texmatrix(X.P,'decimals',opt.decimals,'env',opt.env);
      end
      texcode=char({'N\left(',mustr,',',Pstr,'\right)'});

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

    function P=var(X)
      %VAR is the variance operator
      P=X.P;
      if length(P)>1
        error('X is multi-dimensional, use cov instead')
      end
    end

    function s=std(X)
      %STD is the standard deviation operator
      if length(X.P)>1
        error('X is multi-dimensional, use cov instead')
      end
      s=sqrt(X.P);
    end

    function out=skew(X)
      %SKEW is the skewness operator
      if length(X.P)>1
        error('skew(X) not defined for multi-dimensional variables')
      end
      out=0;
    end

    function out=kurt(X)
      %KURT is the kurtosis operator
      if length(X.P)>1
        error('kurt(X) not defined for multi-dimensional variables')
      end
      out=0;
    end

    function P=cov(X)
      %COV is the covariance operator
      P=X.P;
    end

    %==================
    %----Estimators----
    %==================

    function [Y,msg]=estimate(X,Z)
      %ESTIMATE computes a moment based estimate of Z in the Gaussian class
      %   X=estimate(ndist,Z)
      m=E(Z);
      P=cov(Z);
      Y=ndist(m,P);
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
        error('ndist.rand: n must be a non-negative integer'),
      end
      if length(X.mu)==0
        error('ndist.rand: X must have a specified mean and variance'),
      elseif length(X.mu)==1
        x=X.mu+sqrt(X.P)*randn(n,m);
      else
        [U,D]=svd(X.P);
        x=ones(n,1)*X.mu.'+randn(n,length(X.mu))*sqrt(D)*U';
      end
    end

    function I=erf(X,x)
      %ERF evaluates the error function I(x)=P(X<x) numerically
      %   I=erf(X,x)
      %   The error function is defined as I(x)=int_{-Inf}^{x} p(z) dz
      if nargin<2
        error('Syntax: erf(X,x)')
      end
      if length(X)>1
        error('ERF only defined for scalar stochastic variables')
      end
      I=(erf( (x-mean(X))/std(X)/sqrt(2))+1)/2;
    end


    function x=erfinv(X,I)
      %ERFINV evaluates the inverse error function I(x)=P(X<x) numerically
      %   x=erfinv(X,I)
      %   The function generalizes the standard erfinv function to non-Gaussian distributions.
      %   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz
      %   That is, I(x)=P(X<x)
      if nargin<2
        error('Syntax: erfinv(X,I)')
      end
      if length(X)>1
        error('ERF only defined for scalar stochastic variables')
      end
      if I>1 | I<0
        error('I must be between 0 and 1')
      end
      x=mean(X)+std(X)*erfinv(I*2-1)*sqrt(2);
    end


    function [P,x]=cdf(X,x)
      %CDF is the cumulative density function
      %   P=cdf(X,x) for given (vector) x
      %   [P,x]=cdf(X) automatically chooses an interval for x
      %   p is (length(x),1)

      if length(X)>1
        error('CDF only defined for scalar stochastic variables')
      end
      if nargin<2 | isempty(x)
        m=E(X); s=std(X);
        x=linspace(m-6*s,m+6*s,1000)';
      end
      P=(erf( (x-mean(X))/std(X)/sqrt(2))+1)/2;
    end


    function [p,x]=pdf(X,x)
      %PDF is the probability density function
      %   p=pdf(X,x) for given x (dim (N,nx))
      %   [p,x]=pdf(X) automatically chooses a (N,nx) grid for x
      %   If X is multivariate, then the pdf is evaluated for each row of x, and
      %   length(p) = size(x,1)
      m=X.mu;
      P=X.P;
      nx=length(X);
      if nx==1;
        var=X.P; s=sqrt(var);
        if nargin<2 | isempty(x)
          x=linspace(m-6*s,m+6*s,1000)';
        else
          if (size(x,2)~=length(m)) & length(m)~=1
            error('NDIST.PDF: x must have the same number of columns as length(mu)')
          end
        end
        if isreal(x)
          scale=0.5;
        else
          scale=1;
        end
        p=1/(2*pi*var)^scale*exp(-scale*(x-m).^2/var);
      elseif nx>1
        if nargin>1
          [N,m]=size(x);
          if nx~=m
            error('NDIST.PDF: x must have the same number of columns as length(mu)')
          end
          N=size(x,1);
          if isreal(x)
            scale=0.5;
          else
            scale=1;
          end
          if isequal(X.P,diag(diag(X.P)))
            Pinv=diag(1./diag(X.P));
            const=1/(2*pi*prod(diag((X.P))))^(nx*scale);
            epsi=x.'-repmat(X.mu,1,N);
            %           V=sum(conj(Pinv*epsi).*epsi,1);  % Slower option
            V=sum(conj(repmat(diag(Pinv),1,N).*epsi).*epsi,1);
            p=const*exp(-scale*V)';
          else
            const=1/(2*pi*det(X.P))^(nx/2);
            epsi=x.'-repmat(X.mu,1,N);
            V=sum(conj(pinv(X.P)*epsi).*epsi,1);
            %V=sum(conj(X.P\epsi).*epsi,1);
            %V=sum(epsi'*pinv(X.P)*epsi,1);
            p=const*exp(-scale*V)';
          end
        else
          % Grid 1000 points in nx-dimensional space
          x=[];
          for k=1:nx
            x =[x m(k)+linspace(-4*sqrt(P(k,k)),4*sqrt(P(k,k)),round(10000^(1/nx)))'];
          end
          if isreal(x)
            scale=0.5;
          else
            scale=1;
          end
          N=size(x,1);
          Pinv=inv(X.P);
          const=1/sqrt(2*pi*det(X.P));
          for k=1:N
            for l=1:N
              epsi=[x(k,1);x(l,2)]-X.mu;
              V=epsi'*Pinv*epsi;
              p(k,l)=const*exp(-scale*V);
            end
          end
        end
      end
    end

    %==================
    %----Operators-----
    %==================

    function Y=vertcat(X1,varargin)
      %VERTCAT [X1;X2] is used to create multivariate normal distributions
      if nargin==1
        Y=X1; %%% copy or not??? /ps
        %error('NDIST.VERTCAT: at least two input arguments required')
      else
        X2=varargin{1};
        if isa(X2,'ndist')
          Y=ndist([X1.mu;X2.mu],blkdiag(X1.P,X2.P));
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
      Y=ndist(X.mu,X.P); %%% return a copy or not? alt. just return X; /ps
    end

    function Y=uminus(X)
      %UMINUS unitary minus, Y=-X
      Y=ndist(-X.mu,X.P);
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
        Y=ndist(mu1+mu2,P1+P2);
      elseif isa(X2,'double')  % Y=Xa
        mu1=X1.mu;
        P1=X1.P;
        if length(X2)==1
          mu2=X2*ones(length(mu1),1);
        else
          mu2=X2;
        end
        P2=zeros(length(mu2));
        Y=ndist(mu1+mu2,P1+P2);
      elseif isa(X2,'pdfclass')
        if isa(X2,'gmdist')
	  Y=X2+X1;   % Use the code in gmdist
        elseif isa(X2,'ndist')
          mu1=X1.mu; P1=X1.P;
          mu2=X2.mu; P2=X2.P;
          if size(mu1,1)~=size(mu2,1)
            error('X1 and X2 not compatible for PLUS')
          end
          Y=ndist(mu1+mu2,P1+P2);
        else
          MC=max([100 X1.MC X2.MC]);
          x1=rand(X1,MC,1);
          x2=rand(X2,MC,1);
          Y=empdist(x1+x2);
        end
      else
        error('Incorrect argument to PLUS in NDIST class')
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
        Y=ndist(X1.*X2.mu,diag(X1)*X2.P*diag(X1));
      elseif isa(X2,'double') %Y=XA
        if size(X1.mu,1)~=size(X2,1);
          error('X1 and X2 not compatible for TIMES')
        end
        Y=ndist(X1.mu.*X2,diag(X2)*X1.P*diag(X2));
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
        Y=ndist(X1.mu.*X2,X2.^2*X1.P);
        %elseif isa(X1,'ndist') & isa(X2,'ndist')   % Both Gaussian
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

    function Y=tt1eval(X,f,varargin)
      %TT1EVAL implements Gauss approximation of Y=f(X) based on Taylor expansion
      %   Y=tt1eval(X,f,varargin)
      %
      %   Based on a first order Taylor expansion around mx=E(X):
      %        f(X) = f(mx) + A (X-mx),
      %           A = df/dx(mx)
      %   If X=N(mx,Px), then the approximation is Y=N(f(mx),A*Px*A')
      %
      %   f is an inline function or m-file with syntax f(x,varargin)
      %
      %   Example:
      %     X=ndist(0,1)
      %     Y=tt1eval(X,inline('(x+1).^2'))   % True E(Y)=2, Var(Y)=6
      %     R=90+ndist(0,5);
      %     Phi=pi/4+ndist(0,0.1);
      %     N1=tt1eval([R;Phi],inline('[x(1,:).*cos(x(2,:)); x(1,:).*sin(x(2,:))]'))
      %     plot2(N1,[R*cos(Phi);R*sin(Phi)])

      if isa(f,'inline') || isa(f, 'function_handle')
        %ok
      elseif exist(f,'file')==2
        %ok
      elseif isstr(f)
        try
          f=inline(f);
        catch
          error('ndist.tt1eval: inline(f) failed')
        end
      else
        error('ndist.tt1eval: f must be either a file name, an inline object or a string')
      end

      n=length(X);
      mu=mean(X);
      P=cov(X);
      F=numgrad(f,1,mu,varargin{:});
      Y=ndist(f(mu,varargin{:}),F*P*F');
    end


    function Y=tt2eval(X,f,varargin)
      %TT2EVAL approximates Y=f(X) based on a second order Taylor expansion
      %   Y=tt2eval(X,f,varargin)
      %
      %   Based on a second order Taylor expansion around mx=E(X):
      %        f(X) = f(mx) + A*(X-mx) + (X-mx)'*B*(X-mx)
      %           A = df/dx(mx)
      %           B = d2f/dx2(mx)
      %   If X=N(mx,Px), then the approximation is
      %         Y = N(f(mx)-tr(B*Px), A*Px*A', 0.5*tr(Px*B*Px*B))
      %   These relations have to be interpreted elementwise in f for vector-valued
      %   functions f.
      %
      %   f is an inline function or m-file with syntax f(x,varargin)
      %
      %   Example:
      %     X=ndist(0,1)
      %     Y=tt2eval(X,inline('(x+1).^2'))   % True E(Y)=2, Var(Y)=6
      %     R=90+ndist(0,5);
      %     Phi=pi/4+ndist(0,0.1);
      %     N2=tt2eval([R;Phi],inline('[x(1,:).*cos(x(2,:)); x(1,:).*sin(x(2,:))]'))
      %     plot2(N2,[R*cos(Phi);R*sin(Phi)])

      if isa(f,'inline') || isa(f, 'function_handle')
        %ok
      elseif exist(f,'file')==2
        %ok
      elseif isstr(f)
        try
          f=inline(f);
        catch
          error('ndist.tt2eval: inline(f) failed')
        end
      else
        error('ndist.tt2eval: f must be either a file name, an inline object or a string')
      end

      n=length(X);
      mu=mean(X);
      P=cov(X);
      J=numgrad(f,1,mu,varargin{:});
      H=numhess(f,1,mu,varargin{:});
      muy=f(mu,varargin{:});
      Py=J*P*J';
      nf=size(H,1);
      for i=1:nf
         muy(i)=muy(i)+ 0.5*trace(squeeze(H(i,:,:))*P);
      end
      for i=1:nf
        for j=1:nf
           Py(i,j)=Py(i,j) + 0.5*trace(squeeze(H(i,:,:))*P*squeeze(H(j,:,:))*P);
        end
      end
      Y=ndist(muy,Py);
    end


    function Y=mceval(X,f,varargin)
      %MCEVAL implements Monte Carlo based approximation of Y=f(X)
      %   Y=mceval(X,f,varargin)
      %
      %
      %   f is an inline function or m-file with syntax f(x,varargin)
      %   The number of MC samples is given by X.MC
      %
      %   Example:
      %     X=ndist(0,1)
      %     Y=mceval(X,inline('(x+1).^2'))   % True E(Y)=2, Var(Y)=6
      %     R=90+ndist(0,5);
      %     Phi=pi/4+ndist(0,0.1);
      %     Nut=mceval([R;Phi],inline('[x(1,:).*cos(x(2,:)); x(1,:).*sin(x(2,:))]'))
      %     plot2(Nut,[R*cos(Phi);R*sin(Phi)])

      if isa(f,'inline') || isa(f, 'function_handle')
        %ok
      elseif exist(f,'file')==2
        %ok
      elseif isstr(f)
        try
          f=inline(f);
        catch
          error('ndist.mceval: inline(f) failed')
        end
      else
        error('ndist.mceval: f must be either a file name, an inline object or a string')
      end

% $$$       if nargin<3 | isempty(NMC)
% $$$         NMC=100;
% $$$       end
      x=rand(X,X.MC);
      try
        y = f(x',varargin{:});
        if size(y,2)~=X.MC
              for k=1:X.MC;
                  y(:,k) = f(x(k,:).',varargin{:});
              end
        end
      catch
          try
              for k=1:X.MC;
                  y(:,k) = f(x(k,:).',varargin{:});
              end
          catch
              error(['ndist.mceval: function f cannot be evaluated ' ...
                     'for X'])
          end
      end
      Y=empdist(y.');
      Y=ndist(mean(Y),cov(Y));
    end


    function Y=qteval(X,f,m,varargin)
      %QTEVAL implements quadrature approximation of Y=f(X)
      %   Y=qteval(X,f,m,varargin)
      %
      %
      %   f is an inline function or m-file with syntax f(x,varargin)
      %   m hermite order in the quadrature rule
      %
      %   Example:
      %     X=ndist(0,1)
      %     Y=qteval(X,inline('(x+1).^2'))   % True E(Y)=2, Var(Y)=6
      %     R=90+ndist(0,5);
      %     Phi=pi/4+ndist(0,0.1);
      %     Nqt=qteval([R;Phi],inline('[x(1,:).*cos(x(2,:)); x(1,:).*sin(x(2,:))]'))
      %     plot2(Nqt,[R*cos(Phi);R*sin(Phi)])

      if isa(f,'inline') || isa(f, 'function_handle')
        %ok
      elseif exist(f,'file')==2
        %ok
      elseif isstr(f)
        try
          f=inline(f);
        catch
          error('ndist.qteval: inline(f) failed')
        end
      else
        error('ndist.qteval: f must be either a file name, an inline object or a string')
      end

      if nargin<3 | isempty(m)
        m=10;
      end
      [p,x,w]=hermite(m);
      x=x';
      try
        for k=1:m;
          y(:,k) = f(x(:,k),varargin{:});
        end
      catch
        error('ndist.qteval: function f cannot be evaluated for X')
      end
      ybar=y*w;
      xbar=x*w;
      nx=1; ny=1; %tmp
      Pyy=zeros(ny,ny);
      Pxy=zeros(nx,ny);
      Pxx=zeros(nx,nx);
      for k=1:m
        Pyy=Pyy+w(k)*(y(:,k)-ybar)*(y(:,k)-ybar)';
        Pxy=Pxy+w(k)*(x(:,k)-xbar)*(y(:,k)-ybar)';
        Pxx=Pxx+w(k)*(x(:,k)-xbar)*(x(:,k)-xbar)';
      end
      A=Pxy'*inv(Pxx);
      b=ybar-A*xbar;
      Y=ndist(b,A*Pxx*A');
    end


    function [Y,S,fS]=uteval(X,f,type,par,varargin)
      %UTVEAL implements the unscented transformation of Y=f(X)
      %   [Y,S,fS]=uteval(X,f,type,par,varargin)
      %
      %   The unscented transformation works as follows:
      %   1. The sigma points S are computed. These are the mean and symmetric
      %      deviations around the mean computed from the covariance matrix of X
      %   2. The sigma points are mapped by the function f to fS
      %   3. The Gaussian distribution Y is fitted to the mapped sigma points
      %
      %   f is an inline function or m-file with syntax f(x,varargin)
      %
      %   type  'ut1'|'ut2'|'ct'
      %         For ut1, par=w0 with default w0=1-n/3
      %         For ut2, par=[beta,alpha,kappa] with default [2 1e-3 0]
      %         For ct,  par=[a] with default [1]
      %
      %   Example:
      %     X=ndist(0,1)
      %     Y=uteval(X,inline('(x+1).^2'))    % True E(Y)=2, Var(Y)=6
      %     R=90+ndist(0,5);
      %     Phi=pi/4+ndist(0,0.1);
      %     Nut=uteval([R;Phi],inline('[x(1,:).*cos(x(2,:)); x(1,:).*sin(x(2,:))]'))
      %     plot2(Nut,[R*cos(Phi);R*sin(Phi)])

      %   Based on code by Gustaf Hendeby

      if nargin<3 | isempty(type)
        type='ut1';
      end
      if isa(f,'inline') || isa(f, 'function_handle')
        %ok
      elseif exist(f,'file')==2
        %ok
      elseif isstr(f)
        try
          f=inline(f);
        catch
          error('ndist.uteval: inline(f) failed')
        end
      else
        error('ndist.uteval: f must be either a file name, an inline object or a string')
      end

      n=length(X);
      mu=mean(X);
      Sigma=cov(X);
      % compute sigma points
      N = 2*n + 1;
      switch type
        case 'ut1'
          if nargin<4 | isempty(par)
            w0=1-n/3;     % Optimal for Gaussian
          else
            w0 = par{1};
          end
          wm = zeros(1, N);
          wm(1) = w0;
          wm(2:end) = (1-w0)/(2*n);
          wc = wm;
          scaling = n/(1-w0);
        case  'ut2'
          beta = 2;          % Optimal for Gaussian
          alpha = 1e-3;
          kappa=0;
          if nargin<4 | isempty(par)
            % Default values above
          else
            beta = par(1);
            if length(par)>1
              alpha = par(2);
            end
            if length(par)>2
              kappa = par(3);
            end
          end
          lambda =  alpha^2*(n+kappa) - n;
          w0 = lambda/(alpha^2*(n+kappa));
          wm = zeros(1, N);
          wm(1) = w0;
          wm(2:end) = 1/(2*alpha^2*(n+kappa));
          wc = wm;
          wc(1) = wc(1)+(1-alpha^2+beta);
          scaling = n + lambda;
        case  'ct'
          if nargin<4 | isempty(par)
            a = 1;
          else
            a = par(1);
          end
          wm = zeros(1, N);
          wm(1) = 0;
          wm(2:end) = 1/(2*n);
          wc = wm;
          scaling = n * a^2;
        otherwise
          error(['Type ',type,' of sigma points are not supported. Use ut1, ut2 or ct']) ;
      end
      S = repmat(mu, 1, N);
      %  dev = chol(scaling*Sigma)';
      [Usvd,Dsvd]= svd(scaling*Sigma);
      dev=Usvd*diag(sqrt(diag(Dsvd)));
      S(:, 2:1+n) = S(:, 2:1+n) + dev;
      S(:, 2+n:end) = S(:, 2+n:end) - dev;

      % Transform
      try
        for k=1:size(S,2);
          fS(:,k) = f(S(:,k),varargin{:});
        end
      catch
        error('ndist.uteval: function evaluation f(X,varargin) failed')
      end
      nz = size(fS, 1);
      muz = zeros(nz, 1);
      for i=1:nz
        muz(i) = wm*fS(i, :)';
        fS2(i, :) = fS(i, :) - muz(i);
      end
      Sigmaz = zeros(nz, nz);
      for i=1:size(wc, 2)
        Sigmaz = Sigmaz + wc(i)*fS2(:, i)*fS2(:, i)';
      end
      Y=ndist(muz,Sigmaz);
    end

    %==================
    %----Plots---------
    %==================


    function plot(varargin)
      %PLOT illustrates one or more PDF objects
      %   plot calls pdfplot
      pdfplot(varargin{:},'view','pdf')
    end


    function cdfplot(varargin)
      %CDFPLOT illustrates the CDF of X
      %   plot calls pdfplot
      pdfplot(varargin{:},'view','cdf')
    end


    function erfplot(varargin)
      %CDFPLOT illustrates the error function ERF of X
      %   plot calls pdfplot
      pdfplot(varargin{:},'view','erf')
    end


    function contour(p,ind)
      %CONTOUR illustrates the PDF object p in 2D
      %   contour(p,ind)
      %   ind is the two indices X(ind) to plot (default the first two ones)

      if nargin<2, ind=[1 2]; end
      nx=length(p);
      if nx==1
        error('contour requires the s.v. X to be at least two-dimensional')
      end
      if length(ind)~=2;
        error('index vector ind must be two dimensional')
      end
      %[p,x]=pdf(X(ind));
      %surf(x(:,1),x(:,2),histeq(p))
      %mesh(x(:,1),x(:,2),histeq(p))
      mu=E(p);
      P=cov(p);
      [U,S,V]=svd(P);
      Sroot=sqrt(S);
      Ph=U*Sroot;
      N=100;
      phi=linspace(0,2*pi,N);
      plot(mu(1),mu(2),'+')
      hold on
      plot(mu(1)+Ph(1,:)*[cos(phi);sin(phi)],...
        mu(2)+Ph(2,:)*[cos(phi);sin(phi)],'-')
      plot(mu(1)+2*Ph(1,:)*[cos(phi);sin(phi)],...
        mu(2)+2*Ph(2,:)*[cos(phi);sin(phi)],'-')
      plot(mu(1)+3*Ph(1,:)*[cos(phi);sin(phi)],...
        mu(2)+3*Ph(2,:)*[cos(phi);sin(phi)],'-')
      hold off
      axis('equal')
      legend('mean','1 std','2 std','3 std')
    end


    function surf(p,ind)
      %SURF illustrates the PDF object p in 2D using a surf plot
      %   surf(p,ind)
      %   ind is the two indices X(ind) to plot (default the first two ones)

      if nargin<2, ind=[1 2]; end
      nx=length(p);
      if nx==1
        error('contour requires the s.v. X to be at least two-dimensional')
      end
      if length(ind)~=2;
        error('index vector ind must be two dimensional')
      end
      [p,x]=pdf(p(ind));
      %surf(x(:,1),x(:,2),p)
      surf(x(:,1),x(:,2),p)
    end

  end
end
