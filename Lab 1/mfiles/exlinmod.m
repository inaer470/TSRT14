function m=exlinmod(ex,opt1);
%EXLINMOD generates standard LTI objects
%   m=exlinmod(ex,opt1);
%   Standard examples of LTI objects of different structures are returned
%   depending on the integer ex:
%
%   Ex       Type Time   Description
%   ----------------------------------------------------------
%   tf1c     tf   cont   DC motor
%   tf1d     tf   disc   DC motor
%   tf2c     tf   cont   Slightly undamped second order system
%   tf2d     tf   disc   Slightly undamped second order system
%   motion1D ss   disc   Disrete time motion model in 1D
%   cp1D     ss   disc   Constant position model in 1D. opt1=T {0.5}
%   cv1D     ss   disc   Constant velocity model in 1D. opt1=T {0.5}
%   ca1D     ss   disc   Constant acceleration model in 1D. opt1=T {0.5}
%   cj1D     ss   disc   Constant jerk model in 1D. opt1=T {0.5}
%   cv2D     ss   disc   Constant velocity model in 2D. opt1=T {0.5}
%   ca2D     ss   disc   Constant acceleration model in 2D. opt1=T {0.5}
%   cv3D     ss   disc   Constant velocity model in 3D. opt1=T {0.5}
%   ca3D     ss   disc   Constant acceleration model in 3D. opt1=T {0.5}
%
%   Example:
%      % step response
%      m=exlinmod('tf2c')
%      step(m)
%      Kalman filter target tracking simulation and estimation
%      m=exlinmod('ca3D');
%      z=simulate(m,10);
%      xhat=kalman(m,z);
%      xplot2(z,xhat,'conf',90,[2 3]);
%
%   See also: exltv, getsignal

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

switch lower(ex)
    case 'tf1c'
        m=tf([0 0 1],[1 1 0]);
        m.name='DC motor';
    case 'tf1d'
        m=tf([0 0 1],[1 1 0]);
        m=c2d(m,1);
        m.name='DC motor';
    case 'tf2c'
        m=tf(0.68*[0 1 -0.5],[1 -1.4 0.74],1);
        m=d2c(m);
        m.name='Slightly undamped second order system';
    case 'tf2d'
        m=tf(0.68*[0 1 -0.5],[1 -1.4 0.74],1);
        m.name='Slightly undamped second order system';
    case 'tf3c'
        m=tf(1,[1 1.4 1]);
        m.name='Spring model';
    case 'tf3d'
        m=tf(1,[1 1.4 1]);
        m=c2d(m,1);
        m.name='Spring model';
    case 'tf2duncertain'
        m=tf([0 1 -0.5],[1 -1.4 0.74],2);
    case {'astrom', 'ÿstrÿm'}
        m=tf([1 -0.5],[1 -1.8 0.8],2);
        m.name='The ''ÿstrÿm'' system';
    case 'gripen'
            Ac=[-0.292 8.13 -201 9.77 0 -12.5 17.1;...
                -0.152 -2.54 0.561 -0.0004 0 107 7.68;...
                0.0364 -0.0678 -0.481 0.0012 0 4.67 -7.98;...
                0 1 0.0401 0 0 0 0;...
                0 0 1 0 0 0 0;...
                0 0 0 0 0 -20 0;...
                0 0 0 0 0 0 -20];
            Bc=[0 -2.15; -31.7 0.0274; 0 1.48; 0 0; 0 0; 20 0; 0 20];
            Cc=[0 0 0 1 0 0 0; 0 0 0 0 1 0 0];
            Dc=zeros(2,2);
            m=lss(Ac,Bc,Cc,Dc);
            m.name='JAS 39 Gripen';
            m.ulabel={'aileron command','rudder command'};
            m.ylabel={'roll angle','yaw angle'};
            m.xlabel={'lateral speed','roll rate','yaw rate','roll angle','yaw angle','aileron angle','rudder angle'};
            m.desc='Gripen model from Control Theory bu Glad and Ljung';
    case 'airc'
            Ac=[0 0 1.1320 0 -1; 0 -0.0538 -0.1712 0 0.0705; 0 0 0 1 0;...
                0 0.0485 0 -0.8556 -1.013; 0 -0.2909 0 1.0532 -0.06859];
            Bc=[0 0 0;-0.12 1 0;0 0 0; 4.419 0 -1.665; 1.575 0 -0.0732];
            Cc=[eye(3) zeros(3,2)];
            Dc=zeros(3,3);
            m=lss(Ac,Bc,Cc,Dc);
            m.name='AIRC';
            m.ulabel={'spoiler angle','forward acceleration','elevator angle'};
            m.ylabel={'altitude','forward speed','pitch angle'};
            m.xlabel={'altitude','forward speed','pitch angle','pitch rate','vertical speed'};
            m.desc='Case study in the book of Maciejowski, see Appendix A.1';
     case 'kairc'
            Ac=[zeros(3,5); 1.0968 0.0016 0.2793 -10.7225 0;...
               -0.6950 -0.0008 0.3246 0 -10.7225];
            Bc=[eye(3); 2.1937 0.0031 0.5587; -1.3900 -0.0017 0.6491];
            Cc=[-29.3803 0.0167 1.8067 164.6078 -94.9601;...
                -3.5182 5.0006 -0.2555 20.8859 -11.1922;...
                -78.0162 -0.0028 -33.74 743.9029 232.0489];
            Dc=[-58.7607 0.0333 -3.6134; -7.0365 10.001 -0.5110;...
                -156.0324 -0.0056 -67.4800];
            m=lss(Ac,Bc,Cc,Dc);
            m.name='AIRC controller';
            m.ulabel={'altitude','forward speed','pitch angle'};
            m.ylabel={'spoiler angle','forward acceleration','elevator angle'};
            m.xlabel='';
            m.desc='Design by Maciejowski, see eq (4.171)-(4.174) in his book';
       case 'tgen'
            Ac=[-18.4456 4.2263 -2.2830 0.2260 0.4220 -0.0951;...
                 -4.0977 -6.0706 5.6825 -0.6966 -1.2246 0.2873;...
                 1.4449 1.4336 -2.6477 0.6092 0.8979 -0.2300;...
                 -0.0093 0.2302 -0.5002 -0.1764 -6.3152 0.1350;...
                 -0.0464 -0.3489 0.7238 6.3117 -0.6886 0.3645;...
                 -0.0602 -0.2361 0.2300 0.0915 -0.3214 -0.2087];
            Bc=[-0.2748 3.1463; -0.0501 -9.3737; -0.1550 7.4296;...
                0.0716 -4.9176; -0.0814 -10.2648; 0.0244 13.7943];
            Cc=[0.5971 -0.7697 4.8850 4.8608 -9.8177 -8.8610;...
                3.1013 9.3422 -5.6000 -0.7490 2.9974 10.5719];
            Dc=zeros(2,2);
            m=lss(Ac,Bc,Cc,Dc);
            m.name='TGEN';
            m.ulabel={'throttle valve position','excitation control'};
            m.ylabel={'generator terminal voltage','generator load angle'};
            m.xlabel='';
            m.desc='Case study in the book of Maciejowski, see Appendix A.2';
       case 'motion1d'
            A=[1 1;0 1];
            B=[];
            Bv=[0.5;1];
            C=[1 0];
            D=[];
            pv=1; pe=1;
            m=lss(A,B,C,D,pv*Bv*Bv',pe,1); %'
            m.name='Motion model in 1D';
            m.xlabel={'position','velocity'};
            m.ylabel='position';
       case 'cp1d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
            A=1;
            B=[];
            Bv=T;
            C=1;
            D=[];
            pv=1; pe=1;
            m=lss(A,B,C,D,pv*Bv*Bv',pe,1); %'
            m.name='Motion model in 1D';
            m.xlabel={'position'};
            m.ylabel='position';
       case 'cv1d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
            A=[1 T;0 1];
            B=[];
            Bv=[0.5*T^2;T];
            C=[1 0];
            D=[];
            pv=1; pe=1;
            m=lss(A,B,C,D,pv*Bv*Bv',pe,1); %'
            m.name='Motion model in 1D';
            m.xlabel={'position','velocity'};
            m.ylabel='position';
       case 'ca1d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
	    A=[1 T T^2/2;0 1 T;0 0 1];
            B=[];
	    Bv=[T^3/6;T^2/2;T];
            C=[1 0 0];
            D=[];
            pv=1; pe=1;
            m=lss(A,B,C,D,pv*Bv*Bv',pe,1); %'
            m.name='Motion model in 1D';
	    m.xlabel={'position','velocity','acceleration'};
            m.ylabel='position';
       case 'cj1d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
            A=[1 T T^2/2 T^3/6;0 1 T T^2/2;0 0 1 T;0 0 0 1];
            B=[];
	    Bv=[T^4/24;T^3/6;T^2/2;T];
            C=[1 0 0 0];
            D=[];
            pv=1; pe=1;
            m=lss(A,B,C,D,pv*Bv*Bv',pe,1); %'
            m.name='Motion model in 1D';
	    m.xlabel={'position','velocity','acceleration','jerk'};
            m.ylabel='position';
       case 'cp2d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
            A=[1 0; 0 1];
            B=[T 0; 0 T];
            C=[1 0; 0 1];
            R=0.01*eye(2);
            m=lss(A,[],C,[],B*B',R,1/T);
            m.xlabel={'X','Y'};
            m.ylabel={'X','Y'};
            m.name='Constant position motion model';
       case 'cv2d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
            A=[1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
            B=[T^2/2 0; 0 T^2/2; T 0; 0 T];
            C=[1 0 0 0; 0 1 0 0];
            R=0.01*eye(2);
            m=lss(A,[],C,[],B*B',R,1/T);
            m.xlabel={'X','Y','vX','vY'};
            m.ylabel={'X','Y'};
            m.name='Constant velocity motion model';
       case 'ca2d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
            A=[1 0 T 0 0 0; 0 1 0 T 0 0; 0 0 1 0 T 0; 0 0 0 1 0 T; 0 0 0 0 1 0;0 0 0 0 0 1];
            B=[T^3/6 0; 0 T^3/6; T^2/2 0; 0 T^2/2; T 0; 0 T];
            C=[1 0 0 0 0 0; 0 1 0 0 0 0];
            R=0.01*eye(2);
            m=lss(A,[],C,[],B*B',R,1/T);
            m.xlabel={'X','Y','vX','vY','aX','aY'};
            m.ylabel={'X','Y'};
            m.name='Constant acceleration motion model';
       case 'cv3d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
            A=[1 0 0 T 0 0;...
               0 1 0 0 T 0;...
               0 0 1 0 0 T;...
               0 0 0 1 0 0;...
               0 0 0 0 1 0;...
               0 0 0 0 0 1];
            B=[T^2/2 0 0; 0 T^2/2 0; 0 0 T^2/2; T 0 0; 0 T 0; 0 0 T];
            C=[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0];
            R=0.01*eye(3);
            m=lss(A,[],C,[],B*B',R,1/T);
            m.xlabel={'X','Y','Z','vX','vY','Z'};
            m.ylabel={'X','Y','Z'};
            m.name='Constant velocity 3D motion model';
       case 'ca3d'
            if nargin<2
                T=0.5;
            else
                T=opt1;
            end
            A=[1 0 0 T 0 0 0 0 0;...
               0 1 0 0 T 0 0 0 0;...
               0 0 1 0 0 T 0 0 0;...
               0 0 0 1 0 0 T 0 0;...
               0 0 0 0 1 0 0 T 0;...
               0 0 0 0 0 1 0 0 T;...
               0 0 0 0 0 0 1 0 0;...
               0 0 0 0 0 0 0 1 0;...
               0 0 0 0 0 0 0 0 1];
            B=[T^3/6 0 0; 0 T^3/6 0; 0 0 T^3/6; T^2/2 0 0; 0 T^2/2 0; 0 0 T^2/2; T 0 0; 0 T 0; 0 0 T];
            C=[1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0];
            R=0.01*eye(3);
            m=lss(A,[],C,[],B*B',R,1/T);
            m.xlabel={'X','Y','Z','vX','vY','vZ','aX','aY','aZ'};
            m.ylabel={'X','Y','Z'};
            m.name='Constant acceleration motion model';

       case 'auv'
       A=[.1273, 1.1069; -0.4665, -2.6960];
       B=[1.2888; -2.5942];
       C=eye(2);
       D=[0;0];
       m=lss(A,B,C,D);
       m.name='AUV steering dynamics';
       m.xlabel={'v','r'};
       m.ylabel=m.xlabel;
       m.ulabel='deltar';

end
