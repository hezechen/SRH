function [X, Y] = Eliminate(L,tolerance)
% To delete the double nodes

X = L(:,1) ;
Y = L(:,2) ;
delx = diff(X) ; dely = diff(Y) ;
dist = sqrt(delx.^2+dely.^2) ;
count = 0 ; doop = [];
for i = 1:length(dist)
    if dist(i)<=tolerance
        count = count+1 ;
        doop(count) = i ;
    end
end
geom = [X Y] ;
geom(doop,:) = [] ;
%fprintf('deleted %d points\n',count)
X = geom(:,1) ;
Y = geom(:,2) ;


