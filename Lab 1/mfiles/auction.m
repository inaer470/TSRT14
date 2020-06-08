function [ind1,ind2,m]=auction(D)
%AUCTION implements the auction algorithm
%   [ind1,ind2]=auction(D)
%
%   The auction algorithm minimizes the sum of D(i,ind1(i))
%
%   D     is an (ny,nt) matrix, where ny is the number of observations from
%         the sensor, and nt the number of targets.
%         Entry D(i,j) is the distance from obs y(i) to target x(j)
%   ind1  is the association index vector,
%         so y(i) is associtated to x(ind1(i))
%   ind2  is the reverse association index vector,
%         so y(ind2(j)) is associtated to x(j)
%   m     The maximum
%
%
%   Example:
%   D=rand(4,4)
%   ind=auction(D)

% Copyright: Fredrik Gustafsson and Umut Orguner
%$ Revision: 28-Oct-2019 $


epsilon=0.1; %amount of deviation from the optimal reward

[ny,nt]=size(D);
if (ny>nt)
    error('Number of columns must be greater than or equal to the number of rows');
end

ind2=zeros(1,nt);
ind1=zeros(1,ny);
Dorg=D;
while ~isempty(find(ind1==0)),
    if (nt==1) %if there is only one item
        [maxval ind2]=max(D);%Assign the item to the best customer
        ind1(ind2)=1;%Assign the corresponding customer to the item
    else
        for i=1:ny,
            if ~ind1(i),
                [maxval,maxind]=max(D(i,:));%find maximum element value and its index
                D(i,maxind)=min(D(i,:))-1;%make the maximum minimum to find second maximum
                [secondmaxval,secondmaxind]=max(D(i,:));%find the second maximum value and its index
                D(i,maxind)=maxval; %restore the maximum value

                ind1(i)=maxind; %Assign the customer the item
                if ind2(maxind),%if item is already assigned
                    ind1(ind2(maxind))=0;%unassign the corresponding customer
                end
                ind2(maxind)=i; %Assign the item to the customer
                D(:,maxind)=D(:,maxind)-(maxval-secondmaxval+epsilon);%reduce the item's value
            end
        end
    end
end
m=0;
for i=1:ny
   m=m+Dorg(i,ind1(i));
end
