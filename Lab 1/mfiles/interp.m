function y2=interp(varargin);

%INTERP interpolates y1(t1) to y2(t2)
%   y2=interp(y1,t1,t2,Property1,Value1,...) with y1 vector valued
%   y2=interp(y1,t2,Property1,Value1,...) with y1 as a signal struct
%
%   Interpolation is based on either a band-limited assumption, where
%   perfect reconstruction and re-sampling can be done, a spline
%   interpolation, or using an assumption of intersample behaviour. This
%   can be zero-order hold for piece-wise constant signal, first-order hold
%   for piece-wise linear signal.
%
%   Property   Value       Description
%   ---------------------------------------------------------
%   method    'BL' | {'hold'} | 'spline'
%   degree    {0} | 1      Degree in hold function
%
%   Examples:
%   t1=sort(rand(30,1));  % Non-uniformly sampled data
%   y1=sin(2*pi*t1); t2=0.1:0.01:0.9;
%   y2=interp(y1,t1,t2,'method','hold','degree',1);

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

optvar=[];
if isstruct(varargin{1})
     out='struct';
     sig=varargin{1};
     t2=varargin{2};
     for k=3:length(varargin)
         optvar{k-2}=varargin{k};
     end
     y1=sig.y;
     if isfield(sig,'t')
         t1=sig.t;
     else
         t1=(0:length(y1)-1);
     end
else
     out='vec';
     y1=varargin{1};
     t1=varargin{2};
     t2=varargin{3};
     for k=4:length(varargin)
         optvar{k-3}=varargin{k};
     end
     y1=y1(:);
     t1=t1(:);
     t2=t2(:);
end

opt=struct('method','hold','degree',0);
opt=optset(opt,optvar);

if strcmpi(opt.method,'bl')
    if any(abs(diff(diff(t1)))>1e-12) % Non-uniform sampling
        % Use the relation y1(t1)=sum_k y(kT)sinc((t1-kT)/T)
        %                  y2(t2)=sum_k y(kT)sinc((t2-kT)/T)
        N=length(y1);
        kT=t1(1):(t1(end)-t1(1))/(N-1):t1(end);
        T=kT(2)-kT(1);
        tm=repmat(t1,1,N)-repmat(kT,N,1);
        ind0=find(tm==0);
        ind=find(tm~=0);
        tm(ind0)=tm(ind0)+eps;
        A=sin(pi*tm/T)./(pi*tm/T);
        %ind=find(isnan(A));
        A(ind0)=1;
        ykT=A\y1;
        tm=repmat(resize(t2,length(t2),1),1,N)-repmat(kT,length(t2),1);
        A=sin(pi*tm/T)./(pi*tm/T);
        ind=find(isnan(A));
        A(ind)=1;
        y2=A*ykT;
    else     % Dedicated algorithms for uniform sampling
        % Use the relation y2(t2)=sum_k y1(kT)sinc((t2-kT)/T)
        N=length(y1);
        kT=t1(:)';
        T=kT(2)-kT(1);
        tm=repmat(resize(t2,length(t2),1),1,N)-repmat(kT,length(t2),1);
        A=sin(pi*tm/T)./(pi*tm/T);
        ind=find(isnan(A));
        A(ind)=1;
        y2=A*y1;
    end
elseif strcmpi(opt.method,'hold')
    n=opt.degree;
    if 1%any(abs(diff(diff(t1)))>1e-12) % Non-uniform sampling
        if n==0; % Zero order hold
            for i=1:length(t2);
                [dum,ind]=min(abs(t2(i)-t1));
                y2(i)=y1(ind);
            end
        elseif n==1  % First order hold
            for i=1:length(t2);
                ind=find(t2(i)<t1);
                if ~isempty(ind)
                   ind2=ind(1);
                   ind1=ind2-1;
                   tau=t1(ind2)-t1(ind1);
                   y2(i)=y1(ind1)+(y1(ind2)-y1(ind1))*(t2(i)-t1(ind1))/tau;
                else
                   error('interp: some t2 is larger than t1')
                end
            end
        end
    else     % Dedicated algorithms for uniform sampling
    end
elseif strcmpi(opt.method,'spline')
    y2=interp1(t1,y1,x2,'spline','extrap');
end
y2=y2(:);
if strcmp(out,'struct')
    y2=struct('y',y2,'t',t2);
end
