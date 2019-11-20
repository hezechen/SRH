% Torsion of a circular cross-section
clearvars; clc; close all

%% Input
C  = [0,0]; % Coordinate of center of circle
R  = 0.1;   % Radius of the circle (m)
% Meshing
M  = 50;    % Number of points along the circumference
N  = 20;    % Number of points along the radius
% Material for layers
r1=2/3*R;   % Radius of 1st circle (m)
r2=3/4*R;   % Radius of 2nd circle (m)
r3=R;       % Radius of 3rd circle (m)
G1=1e7;     % Shear modulus of material in 1st circle (Pa)
G2=1e7;     % Shear modulus of material between 1st circle and 2nd circle (Pa)
G3=1e7;     % Shear modulus of material between 2st circle and 3nd circle (Pa)
% Loading
theta=0.001;% rad/m

%% Mesh square at the center of circle
% Square at the centre
S = R/4;
P1 = [C(1)-S C(2)+S];
P2 = [C(1)+S C(2)+S];
P3 = [C(1)-S C(2)-S];
P4 = [C(1)+S C(2)-S];
% Draw lines
NS = round(M/4) ;
t = linspace(0,1,NS) ;
L1 = zeros(NS,2) ;
L2 = zeros(NS,2) ;
for i = 1:NS
    L1(i,:) = P1+t(i)*(P2-P1) ;
    L2(i,:) = P3+t(i)*(P4-P3) ;
end
% Generate mesh along square
X1 = zeros(NS,NS) ;
Y1 = zeros(NS,NS) ;
for i = 1:NS
    for j = 1:NS
        X1(i,j) = L1(i,1)*(1-t(j))+L2(i,1)*t(j) ;
        Y1(i,j) = L1(i,2)*(1-t(j))+L2(i,2)*t(j) ;
    end
end

%% Mesh region between square and circle boundary
% Get boundaries of square
x1 = [X1(1,:)' ; X1(:,end) ; flipud(X1(end,:)') ;flipud(X1(:,1))] ;
y1 = [Y1(1,:)' ; Y1(:,end) ; flipud(Y1(end,:)') ;flipud(Y1(:,1))] ;
% Delete double points if any
[x1, y1] = Eliminate([x1 y1],10^-3) ;
% Circle coordinates
M = length(x1) ;
th0 = linspace(0,2*pi,M) ;
% Arrange th
idx = find(th0>=3*pi/4) ;
th = [th0(idx) th0(2:idx(1))] ;
x2 = C(1)+R*cos(th') ;
y2 = C(2)+R*sin(th') ;
% Mesh the square region
t = linspace(0,1,N) ;
X2 = zeros(M,N) ;
Y2 = zeros(M,N) ;
for i = 1:M
    for j = 1:N
        X2(i,j) = x1(i)*(1-t(j))+x2(i)*t(j) ;
        Y2(i,j) = y1(i)*(1-t(j))+y2(i)*t(j) ;
    end
end

X1=X1';
Y1=Y1';
X2=X2';
Y2=Y2';

%% Plot meshing
figure('Name','Mesh grid')
plotgrid(X1,Y1) ;
plotgrid(X2,Y2) ;

%% Local system
Ntotal=(NS-1)*(NS-1)+(M-1)*(N-1);
V(1:Ntotal)=struct('G',[],'Part',[],'cont',[],'Vertices',[],'Area',[],'Cent',[],'A',[],'coJi',[],'Ji',[],'Length',[],'n',[],'c5c6',[],'ab',[],'L',[],'Lc',[],'u',[],'W',[],'Jc',[],'strainc',[],'stressc',[]);

% Inner square
for i=1:NS-1
    for j=1:NS-1
        order=(i-1)*(NS-1)+j;
        V(order).Part=1;
        V(order).cont=[1 1 1 1];
        if i==1
            V(order).cont(3)=2;
        end
        if i==NS-1
            V(order).cont(1)=2;
        end
        if j==1
            V(order).cont(4)=2;
        end
        if j==NS-1
            V(order).cont(2)=2;
        end
        V(order).Vertices=[X1(i+1,j) Y1(i+1,j);X1(i+1,j+1) Y1(i+1,j+1);X1(i,j+1) Y1(i,j+1);X1(i,j) Y1(i,j)];
        
    end
end

% Outer the square
for i=1:N-1
    for j=1:M-1
        order=(NS-1)*(NS-1)+(i-1)*(M-1)+j;
        V(order).Part=2;
        V(order).cont=[4 4 4 4];
        if i==1
            V(order).cont(1)=3;
        end
        if i==N-1
            V(order).cont(3)=0;
        end
        if j==1
            V(order).cont(2)=5;
        end
        if j==M-1
            V(order).cont(4)=5;
        end
        V(order).Vertices=[X2(i,j+1) Y2(i,j+1); X2(i,j) Y2(i,j);X2(i+1,j) Y2(i+1,j);X2(i+1,j+1) Y2(i+1,j+1)];
    end
end

% Give information to each element
for order=1:Ntotal
    V(order).Area=1/2*(V(order).Vertices(1,1)*V(order).Vertices(2,2)+V(order).Vertices(2,1)*V(order).Vertices(3,2)+V(order).Vertices(3,1)*V(order).Vertices(4,2)+V(order).Vertices(4,1)*V(order).Vertices(1,2)-V(order).Vertices(2,1)*V(order).Vertices(1,2)-V(order).Vertices(3,1)*V(order).Vertices(2,2)-V(order).Vertices(4,1)*V(order).Vertices(3,2)-V(order).Vertices(1,1)*V(order).Vertices(4,2));
    
    V(order).Cent=[(V(order).Vertices(1,1)+V(order).Vertices(2,1)+V(order).Vertices(3,1)+V(order).Vertices(4,1))/4;
        (V(order).Vertices(1,2)+V(order).Vertices(2,2)+V(order).Vertices(3,2)+V(order).Vertices(4,2))/4];
    
    if V(order).Cent(1)^2+V(order).Cent(2)^2<=r1^2
        V(order).G=[G1,G1];   % Material 1
    elseif V(order).Cent(1)^2+V(order).Cent(2)^2<=r2^2
        V(order).G=[G2,G2];   % Material 2
    elseif V(order).Cent(1)^2+V(order).Cent(2)^2<=r3^2
        V(order).G=[G3,G3];   % Material 3
    end
    
    V(order).A=1/4*[-V(order).Vertices(1,1)+V(order).Vertices(2,1)+V(order).Vertices(3,1)-V(order).Vertices(4,1);
        V(order).Vertices(1,1)-V(order).Vertices(2,1)+V(order).Vertices(3,1)-V(order).Vertices(4,1);
        -V(order).Vertices(1,1)-V(order).Vertices(2,1)+V(order).Vertices(3,1)+V(order).Vertices(4,1);
        -V(order).Vertices(1,2)+V(order).Vertices(2,2)+V(order).Vertices(3,2)-V(order).Vertices(4,2);
        V(order).Vertices(1,2)-V(order).Vertices(2,2)+V(order).Vertices(3,2)-V(order).Vertices(4,2);
        -V(order).Vertices(1,2)-V(order).Vertices(2,2)+V(order).Vertices(3,2)+V(order).Vertices(4,2);];
    
    V(order).coJi=1/(V(order).A(1)*V(order).A(6)-V(order).A(3)*V(order).A(4));
    V(order).Ji=V(order).coJi*[V(order).A(6) -V(order).A(4); -V(order).A(3) V(order).A(1)];
    
    V(order).Length=[sqrt((V(order).Vertices(2,1)-V(order).Vertices(1,1))^2+(V(order).Vertices(2,2)-V(order).Vertices(1,2))^2);
        sqrt((V(order).Vertices(3,1)-V(order).Vertices(2,1))^2+(V(order).Vertices(3,2)-V(order).Vertices(2,2))^2);
        sqrt((V(order).Vertices(4,1)-V(order).Vertices(3,1))^2+(V(order).Vertices(4,2)-V(order).Vertices(3,2))^2);
        sqrt((V(order).Vertices(1,1)-V(order).Vertices(4,1))^2+(V(order).Vertices(1,2)-V(order).Vertices(4,2))^2)];
    
    V(order).n(1,1)=(V(order).Vertices(2,2)-V(order).Vertices(1,2))/V(order).Length(1);
    V(order).n(1,2)=(-V(order).Vertices(2,1)+V(order).Vertices(1,1))/V(order).Length(1);
    V(order).n(2,1)=(V(order).Vertices(3,2)-V(order).Vertices(2,2))/V(order).Length(2);
    V(order).n(2,2)=(-V(order).Vertices(3,1)+V(order).Vertices(2,1))/V(order).Length(2);
    V(order).n(3,1)=(V(order).Vertices(4,2)-V(order).Vertices(3,2))/V(order).Length(3);
    V(order).n(3,2)=(-V(order).Vertices(4,1)+V(order).Vertices(3,1))/V(order).Length(3);
    V(order).n(4,1)=(V(order).Vertices(1,2)-V(order).Vertices(4,2))/V(order).Length(4);
    V(order).n(4,2)=(-V(order).Vertices(1,1)+V(order).Vertices(4,1))/V(order).Length(4);
    
    % Self-defined quantities
    V(order).c5c6(:,1)=V(order).G(1)*theta*V(order).Cent(2)*[V(order).n(1,1); V(order).n(2,1); V(order).n(3,1); V(order).n(4,1)];
    V(order).c5c6(:,2)=V(order).G(2)*theta*V(order).Cent(1)*[V(order).n(1,2); V(order).n(2,2); V(order).n(3,2); V(order).n(4,2)];
    
    V(order).ab(:,1)=V(order).G(1)*V(order).A(4)*[V(order).n(1,1); V(order).n(2,1); V(order).n(3,1); V(order).n(4,1)]-V(order).G(2)*V(order).A(1)*[V(order).n(1,2); V(order).n(2,2); V(order).n(3,2); V(order).n(4,2)];
    V(order).ab(:,2)=V(order).G(1)*V(order).A(6)*[V(order).n(1,1); V(order).n(2,1); V(order).n(3,1); V(order).n(4,1)]-V(order).G(2)*V(order).A(3)*[V(order).n(1,2); V(order).n(2,2); V(order).n(3,2); V(order).n(4,2)];
    
    % Local stiffness matrix
    V(order).L=V(order).coJi*[2*V(order).ab(1,1) 1/2*V(order).ab(1,2) V(order).ab(1,1) -1/2*V(order).ab(1,2);
        1/2*V(order).ab(2,1) 2*V(order).ab(2,2) -1/2*V(order).ab(2,1) V(order).ab(2,2);
        -V(order).ab(3,1) 1/2*V(order).ab(3,2) -2*V(order).ab(3,1) -1/2*V(order).ab(3,2);
        1/2*V(order).ab(4,1) -V(order).ab(4,2) -1/2*V(order).ab(4,1) -2*V(order).ab(4,2);]...
        -V(order).coJi/(V(order).ab(1,1)*V(order).Length(1)+V(order).ab(2,2)*V(order).Length(2)-V(order).ab(3,1)*V(order).Length(3)-V(order).ab(4,2)*V(order).Length(4))*[V(order).ab(1,1);V(order).ab(2,2);V(order).ab(3,1);V(order).ab(4,2)]...
        *[2*V(order).ab(1,1)*V(order).Length(1)+1/2*V(order).ab(2,1)*V(order).Length(2)-V(order).ab(3,1)*V(order).Length(3)+1/2*V(order).ab(4,1)*V(order).Length(4),...
        1/2*V(order).ab(1,2)*V(order).Length(1)+2*V(order).ab(2,2)*V(order).Length(2)+1/2*V(order).ab(3,2)*V(order).Length(3)-V(order).ab(4,2)*V(order).Length(4),...
        V(order).ab(1,1)*V(order).Length(1)-1/2*V(order).ab(2,1)*V(order).Length(2)-2*V(order).ab(3,1)*V(order).Length(3)-1/2*V(order).ab(4,1)*V(order).Length(4),...
        1/2*V(order).ab(1,2)*V(order).Length(1)+V(order).ab(2,2)*V(order).Length(2)-1/2*V(order).ab(3,2)*V(order).Length(3)-2*V(order).ab(4,2)*V(order).Length(4)];
    
    % Local constant
    V(order).Lc=[-V(order).c5c6(1,1)+V(order).c5c6(1,2);
        -V(order).c5c6(2,1)+V(order).c5c6(2,2);
        -V(order).c5c6(3,1)+V(order).c5c6(3,2);
        -V(order).c5c6(4,1)+V(order).c5c6(4,2);]...
        -[V(order).Length(1) V(order).Length(2) V(order).Length(3) V(order).Length(4)]*[-V(order).c5c6(1,1)+V(order).c5c6(1,2);
        -V(order).c5c6(2,1)+V(order).c5c6(2,2);
        -V(order).c5c6(3,1)+V(order).c5c6(3,2);
        -V(order).c5c6(4,1)+V(order).c5c6(4,2)]...
        *1/(V(order).ab(1,1)*V(order).Length(1)+V(order).ab(2,2)*V(order).Length(2)-V(order).ab(3,1)*V(order).Length(3)-V(order).ab(4,2)*V(order).Length(4))...
        *[V(order).ab(1,1);V(order).ab(2,2);V(order).ab(3,1);V(order).ab(4,2)];
end

%% Global system
K=zeros(8*Ntotal+1,4*Ntotal);
C=zeros(8*Ntotal+1,1);

% Traction continuity and boundary condition
for order=1:Ntotal
    % Face 1
    if V(order).cont(1)==1
        K(4*(order-1)+1,4*(order-1)+1)=V(order).L(1,1);
        K(4*(order-1)+1,4*(order-1)+2)=V(order).L(1,2);
        K(4*(order-1)+1,4*(order-1)+3)=V(order).L(1,3);
        K(4*(order-1)+1,4*(order-1)+4)=V(order).L(1,4);
        K(4*(order-1)+1,4*(order+NS-2)+1)=V(order+NS-1).L(3,1);
        K(4*(order-1)+1,4*(order+NS-2)+2)=V(order+NS-1).L(3,2);
        K(4*(order-1)+1,4*(order+NS-2)+3)=V(order+NS-1).L(3,3);
        K(4*(order-1)+1,4*(order+NS-2)+4)=V(order+NS-1).L(3,4);
        C(4*(order-1)+1)=-V(order).Lc(1)-V(order+NS-1).Lc(3);
    elseif V(order).cont(1)==2
        K(4*(order-1)+1,4*(order-1)+1)=V(order).L(1,1);
        K(4*(order-1)+1,4*(order-1)+2)=V(order).L(1,2);
        K(4*(order-1)+1,4*(order-1)+3)=V(order).L(1,3);
        K(4*(order-1)+1,4*(order-1)+4)=V(order).L(1,4);
        K(4*(order-1)+1,4*(order+(M-1)/2-1)+1)=V(order+(M-1)/2).L(1,1);
        K(4*(order-1)+1,4*(order+(M-1)/2-1)+2)=V(order+(M-1)/2).L(1,2);
        K(4*(order-1)+1,4*(order+(M-1)/2-1)+3)=V(order+(M-1)/2).L(1,3);
        K(4*(order-1)+1,4*(order+(M-1)/2-1)+4)=V(order+(M-1)/2).L(1,4);
        C(4*(order-1)+1)=-V(order).Lc(1)-V(order+(M-1)/2).Lc(1);
    elseif V(order).cont(1)==3
        if order-(NS-1)^2<=(M-1)/4
            K(4*(order-1)+1,4*(order-1)+1)=V(order).L(1,1);
            K(4*(order-1)+1,4*(order-1)+2)=V(order).L(1,2);
            K(4*(order-1)+1,4*(order-1)+3)=V(order).L(1,3);
            K(4*(order-1)+1,4*(order-1)+4)=V(order).L(1,4);
            K(4*(order-1)+1,4*(order-(NS-1)^2-1)*(NS-1)+1)=V((order-(NS-1)^2-1)*(NS-1)+1).L(4,1);
            K(4*(order-1)+1,4*(order-(NS-1)^2-1)*(NS-1)+2)=V((order-(NS-1)^2-1)*(NS-1)+1).L(4,2);
            K(4*(order-1)+1,4*(order-(NS-1)^2-1)*(NS-1)+3)=V((order-(NS-1)^2-1)*(NS-1)+1).L(4,3);
            K(4*(order-1)+1,4*(order-(NS-1)^2-1)*(NS-1)+4)=V((order-(NS-1)^2-1)*(NS-1)+1).L(4,4);
            C(4*(order-1)+1)=-V(order).Lc(1)-V((order-(NS-1)^2-1)*(NS-1)+1).Lc(4);
        elseif order-(NS-1)^2<=(M-1)/2
            K(4*(order-1)+1,4*(order-1)+1)=V(order).L(1,1);
            K(4*(order-1)+1,4*(order-1)+2)=V(order).L(1,2);
            K(4*(order-1)+1,4*(order-1)+3)=V(order).L(1,3);
            K(4*(order-1)+1,4*(order-1)+4)=V(order).L(1,4);
            K(4*(order-1)+1,4*(order-(M-1)/2-1)+1)=V(order-(M-1)/2).L(1,1);
            K(4*(order-1)+1,4*(order-(M-1)/2-1)+2)=V(order-(M-1)/2).L(1,2);
            K(4*(order-1)+1,4*(order-(M-1)/2-1)+3)=V(order-(M-1)/2).L(1,3);
            K(4*(order-1)+1,4*(order-(M-1)/2-1)+4)=V(order-(M-1)/2).L(1,4);
            C(4*(order-1)+1)=-V(order).Lc(1)-V(order-(M-1)/2).Lc(1);
        elseif order-(NS-1)^2<=3*(M-1)/4
            K(4*(order-1)+1,4*(order-1)+1)=V(order).L(1,1);
            K(4*(order-1)+1,4*(order-1)+2)=V(order).L(1,2);
            K(4*(order-1)+1,4*(order-1)+3)=V(order).L(1,3);
            K(4*(order-1)+1,4*(order-1)+4)=V(order).L(1,4);
            K(4*(order-1)+1,4*((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)-1)+1)=V((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)).L(2,1);
            K(4*(order-1)+1,4*((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)-1)+2)=V((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)).L(2,2);
            K(4*(order-1)+1,4*((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)-1)+3)=V((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)).L(2,3);
            K(4*(order-1)+1,4*((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)-1)+4)=V((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)).L(2,4);
            C(4*(order-1)+1)=-V(order).Lc(1)-V((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)).Lc(2);
        elseif order-(NS-1)^2<=(M-1)
            K(4*(order-1)+1,4*(order-1)+1)=V(order).L(1,1);
            K(4*(order-1)+1,4*(order-1)+2)=V(order).L(1,2);
            K(4*(order-1)+1,4*(order-1)+3)=V(order).L(1,3);
            K(4*(order-1)+1,4*(order-1)+4)=V(order).L(1,4);
            K(4*(order-1)+1,4*((NS-1)^2+M-1-order)+1)=V((NS-1)^2+M-1-order+1).L(3,1);
            K(4*(order-1)+1,4*((NS-1)^2+M-1-order)+2)=V((NS-1)^2+M-1-order+1).L(3,2);
            K(4*(order-1)+1,4*((NS-1)^2+M-1-order)+3)=V((NS-1)^2+M-1-order+1).L(3,3);
            K(4*(order-1)+1,4*((NS-1)^2+M-1-order)+4)=V((NS-1)^2+M-1-order+1).L(3,4);
            C(4*(order-1)+1)=-V(order).Lc(1)-V((NS-1)^2+M-1-order+1).Lc(3);
        end
    elseif V(order).cont(1)==4
        K(4*(order-1)+1,4*(order-1)+1)=V(order).L(1,1);
        K(4*(order-1)+1,4*(order-1)+2)=V(order).L(1,2);
        K(4*(order-1)+1,4*(order-1)+3)=V(order).L(1,3);
        K(4*(order-1)+1,4*(order-1)+4)=V(order).L(1,4);
        K(4*(order-1)+1,4*(order-M)+1)=V(order-M+1).L(3,1);
        K(4*(order-1)+1,4*(order-M)+2)=V(order-M+1).L(3,2);
        K(4*(order-1)+1,4*(order-M)+3)=V(order-M+1).L(3,3);
        K(4*(order-1)+1,4*(order-M)+4)=V(order-M+1).L(3,4);
        C(4*(order-1)+1)=-V(order).Lc(1)-V(order-M+1).Lc(3);
    end
    % Face 2
    if V(order).cont(2)==1
        K(4*(order-1)+2,4*(order-1)+1)=V(order).L(2,1);
        K(4*(order-1)+2,4*(order-1)+2)=V(order).L(2,2);
        K(4*(order-1)+2,4*(order-1)+3)=V(order).L(2,3);
        K(4*(order-1)+2,4*(order-1)+4)=V(order).L(2,4);
        K(4*(order-1)+2,4*order+1)=V(order+1).L(4,1);
        K(4*(order-1)+2,4*order+2)=V(order+1).L(4,2);
        K(4*(order-1)+2,4*order+3)=V(order+1).L(4,3);
        K(4*(order-1)+2,4*order+4)=V(order+1).L(4,4);
        C(4*(order-1)+2)=-V(order).Lc(2)-V(order+1).Lc(4);
    elseif V(order).cont(2)==2
        K(4*(order-1)+2,4*(order-1)+1)=V(order).L(2,1);
        K(4*(order-1)+2,4*(order-1)+2)=V(order).L(2,2);
        K(4*(order-1)+2,4*(order-1)+3)=V(order).L(2,3);
        K(4*(order-1)+2,4*(order-1)+4)=V(order).L(2,4);
        K(4*(order-1)+2,4*((NS-1)^2+3*(M-1)/4-order/(NS-1))+1)=V((NS-1)^2+3*(M-1)/4+1-order/(NS-1)).L(1,1);
        K(4*(order-1)+2,4*((NS-1)^2+3*(M-1)/4-order/(NS-1))+2)=V((NS-1)^2+3*(M-1)/4+1-order/(NS-1)).L(1,2);
        K(4*(order-1)+2,4*((NS-1)^2+3*(M-1)/4-order/(NS-1))+3)=V((NS-1)^2+3*(M-1)/4+1-order/(NS-1)).L(1,3);
        K(4*(order-1)+2,4*((NS-1)^2+3*(M-1)/4-order/(NS-1))+4)=V((NS-1)^2+3*(M-1)/4+1-order/(NS-1)).L(1,4);
        C(4*(order-1)+2)=-V(order).Lc(2)-V((NS-1)^2+3*(M-1)/4+1-order/(NS-1)).Lc(1);
    elseif V(order).cont(2)==4
        K(4*(order-1)+2,4*(order-1)+1)=V(order).L(2,1);
        K(4*(order-1)+2,4*(order-1)+2)=V(order).L(2,2);
        K(4*(order-1)+2,4*(order-1)+3)=V(order).L(2,3);
        K(4*(order-1)+2,4*(order-1)+4)=V(order).L(2,4);
        K(4*(order-1)+2,4*(order-2)+1)=V(order-1).L(4,1);
        K(4*(order-1)+2,4*(order-2)+2)=V(order-1).L(4,2);
        K(4*(order-1)+2,4*(order-2)+3)=V(order-1).L(4,3);
        K(4*(order-1)+2,4*(order-2)+4)=V(order-1).L(4,4);
        C(4*(order-1)+2)=-V(order).Lc(2)-V(order-1).Lc(4);
    elseif V(order).cont(2)==5
        K(4*(order-1)+2,4*(order-1)+1)=V(order).L(2,1);
        K(4*(order-1)+2,4*(order-1)+2)=V(order).L(2,2);
        K(4*(order-1)+2,4*(order-1)+3)=V(order).L(2,3);
        K(4*(order-1)+2,4*(order-1)+4)=V(order).L(2,4);
        K(4*(order-1)+2,4*(order+M-3)+1)=V(order+M-2).L(4,1);
        K(4*(order-1)+2,4*(order+M-3)+2)=V(order+M-2).L(4,2);
        K(4*(order-1)+2,4*(order+M-3)+3)=V(order+M-2).L(4,3);
        K(4*(order-1)+2,4*(order+M-3)+4)=V(order+M-2).L(4,4);
        C(4*(order-1)+2)=-V(order).Lc(2)-V(order+M-2).Lc(4);
    end
    % Face 3
    if V(order).cont(3)==1
        K(4*(order-1)+3,4*(order-1)+1)=V(order).L(3,1);
        K(4*(order-1)+3,4*(order-1)+2)=V(order).L(3,2);
        K(4*(order-1)+3,4*(order-1)+3)=V(order).L(3,3);
        K(4*(order-1)+3,4*(order-1)+4)=V(order).L(3,4);
        K(4*(order-1)+3,4*(order-NS)+1)=V(order-NS+1).L(1,1);
        K(4*(order-1)+3,4*(order-NS)+2)=V(order-NS+1).L(1,2);
        K(4*(order-1)+3,4*(order-NS)+3)=V(order-NS+1).L(1,3);
        K(4*(order-1)+3,4*(order-NS)+4)=V(order-NS+1).L(1,4);
        C(4*(order-1)+3)=-V(order).Lc(3)-V(order-NS+1).Lc(1);
    elseif V(order).cont(3)==2
        K(4*(order-1)+3,4*(order-1)+1)=V(order).L(3,1);
        K(4*(order-1)+3,4*(order-1)+2)=V(order).L(3,2);
        K(4*(order-1)+3,4*(order-1)+3)=V(order).L(3,3);
        K(4*(order-1)+3,4*(order-1)+4)=V(order).L(3,4);
        K(4*(order-1)+3,4*((NS-1)^2+M-order-1)+1)=V((NS-1)^2+M-order).L(1,1);
        K(4*(order-1)+3,4*((NS-1)^2+M-order-1)+2)=V((NS-1)^2+M-order).L(1,2);
        K(4*(order-1)+3,4*((NS-1)^2+M-order-1)+3)=V((NS-1)^2+M-order).L(1,3);
        K(4*(order-1)+3,4*((NS-1)^2+M-order-1)+4)=V((NS-1)^2+M-order).L(1,4);
        C(4*(order-1)+3)=-V(order).Lc(3)-V((NS-1)^2+M-order).Lc(1);
    elseif V(order).cont(3)==4
        K(4*(order-1)+3,4*(order-1)+1)=V(order).L(3,1);
        K(4*(order-1)+3,4*(order-1)+2)=V(order).L(3,2);
        K(4*(order-1)+3,4*(order-1)+3)=V(order).L(3,3);
        K(4*(order-1)+3,4*(order-1)+4)=V(order).L(3,4);
        K(4*(order-1)+3,4*(order+M-2)+1)=V(order+M-1).L(1,1);
        K(4*(order-1)+3,4*(order+M-2)+2)=V(order+M-1).L(1,2);
        K(4*(order-1)+3,4*(order+M-2)+3)=V(order+M-1).L(1,3);
        K(4*(order-1)+3,4*(order+M-2)+4)=V(order+M-1).L(1,4);
        C(4*(order-1)+3)=-V(order).Lc(3)-V(order+M-1).Lc(1);
    elseif V(order).cont(3)==0 % boundary condition
        K(4*(order-1)+3,4*(order-1)+1)=V(order).L(3,1);
        K(4*(order-1)+3,4*(order-1)+2)=V(order).L(3,2);
        K(4*(order-1)+3,4*(order-1)+3)=V(order).L(3,3);
        K(4*(order-1)+3,4*(order-1)+4)=V(order).L(3,4);
        C(4*(order-1)+3)=-V(order).Lc(3);
    end
    % Face 4
    if V(order).cont(4)==1
        K(4*(order-1)+4,4*(order-1)+1)=V(order).L(4,1);
        K(4*(order-1)+4,4*(order-1)+2)=V(order).L(4,2);
        K(4*(order-1)+4,4*(order-1)+3)=V(order).L(4,3);
        K(4*(order-1)+4,4*(order-1)+4)=V(order).L(4,4);
        K(4*(order-1)+4,4*(order-2)+1)=V(order-1).L(2,1);
        K(4*(order-1)+4,4*(order-2)+2)=V(order-1).L(2,2);
        K(4*(order-1)+4,4*(order-2)+3)=V(order-1).L(2,3);
        K(4*(order-1)+4,4*(order-2)+4)=V(order-1).L(2,4);
        C(4*(order-1)+4)=-V(order).Lc(4)-V(order-1).Lc(2);
    elseif V(order).cont(4)==2
        K(4*(order-1)+4,4*(order-1)+1)=V(order).L(4,1);
        K(4*(order-1)+4,4*(order-1)+2)=V(order).L(4,2);
        K(4*(order-1)+4,4*(order-1)+3)=V(order).L(4,3);
        K(4*(order-1)+4,4*(order-1)+4)=V(order).L(4,4);
        K(4*(order-1)+4,4*((NS-1)^2+floor(order/(NS-1)))+1)=V((NS-1)^2+floor(order/(NS-1))+1).L(1,1);
        K(4*(order-1)+4,4*((NS-1)^2+floor(order/(NS-1)))+2)=V((NS-1)^2+floor(order/(NS-1))+1).L(1,2);
        K(4*(order-1)+4,4*((NS-1)^2+floor(order/(NS-1)))+3)=V((NS-1)^2+floor(order/(NS-1))+1).L(1,3);
        K(4*(order-1)+4,4*((NS-1)^2+floor(order/(NS-1)))+4)=V((NS-1)^2+floor(order/(NS-1))+1).L(1,4);
        C(4*(order-1)+4)=-V(order).Lc(4)-V((NS-1)^2+floor(order/(NS-1))+1).Lc(1);
    elseif V(order).cont(4)==4
        K(4*(order-1)+4,4*(order-1)+1)=V(order).L(4,1);
        K(4*(order-1)+4,4*(order-1)+2)=V(order).L(4,2);
        K(4*(order-1)+4,4*(order-1)+3)=V(order).L(4,3);
        K(4*(order-1)+4,4*(order-1)+4)=V(order).L(4,4);
        K(4*(order-1)+4,4*order+1)=V(order+1).L(2,1);
        K(4*(order-1)+4,4*order+2)=V(order+1).L(2,2);
        K(4*(order-1)+4,4*order+3)=V(order+1).L(2,3);
        K(4*(order-1)+4,4*order+4)=V(order+1).L(2,4);
        C(4*(order-1)+4)=-V(order).Lc(4)-V(order+1).Lc(2);
    elseif V(order).cont(4)==5
        K(4*(order-1)+4,4*(order-1)+1)=V(order).L(4,1);
        K(4*(order-1)+4,4*(order-1)+2)=V(order).L(4,2);
        K(4*(order-1)+4,4*(order-1)+3)=V(order).L(4,3);
        K(4*(order-1)+4,4*(order-1)+4)=V(order).L(4,4);
        K(4*(order-1)+4,4*(order-M+1)+1)=V(order-M+2).L(2,1);
        K(4*(order-1)+4,4*(order-M+1)+2)=V(order-M+2).L(2,2);
        K(4*(order-1)+4,4*(order-M+1)+3)=V(order-M+2).L(2,3);
        K(4*(order-1)+4,4*(order-M+1)+4)=V(order-M+2).L(2,4);
        C(4*(order-1)+4)=-V(order).Lc(4)-V(order-M+2).Lc(2);
    end
    % Constraint
    if V(order).cont(3)==0
        K(4*Ntotal+1,4*(order-1)+3)=1;
        C(4*Ntotal+1)=0;
    end
end

% Displacement continuity
for order=1:Ntotal
    % Face 1
    if V(order).cont(1)==1
        K(4*Ntotal+1+4*(order-1)+1,4*(order-1)+1)=1;
        K(4*Ntotal+1+4*(order-1)+1,4*(order+NS-2)+3)=-1;
        C(4*Ntotal+1+4*(order-1)+1)=0;
    elseif V(order).cont(1)==2
        K(4*Ntotal+1+4*(order-1)+1,4*(order-1)+1)=1;
        K(4*Ntotal+1+4*(order-1)+1,4*(order+(M-1)/2-1)+1)=-1;
        C(4*Ntotal+1+4*(order-1)+1)=-0;
    elseif V(order).cont(1)==3
        if order-(NS-1)^2<=(M-1)/4
            K(4*Ntotal+1+4*(order-1)+1,4*(order-1)+1)=1;
            K(4*Ntotal+1+4*(order-1)+1,4*(order-(NS-1)^2-1)*(NS-1)+1)=-1;
            C(4*Ntotal+1+4*(order-1)+1)=0;
        elseif order-(NS-1)^2<=(M-1)/2
            K(4*Ntotal+1+4*(order-1)+1,4*(order-1)+1)=1;
            K(4*Ntotal+1+4*(order-1)+1,4*(order-(M-1)/2-1)+1)=-1;
            C(4*Ntotal+1+4*(order-1)+1)=0;
        elseif order-(NS-1)^2<=3*(M-1)/4
            K(4*Ntotal+1+4*(order-1)+1,4*(order-1)+1)=1;
            K(4*Ntotal+1+4*(order-1)+1,4*((NS-(order-(NS-1)^2-(M-1)/2))*(NS-1)-1)+2)=-1;
            C(4*Ntotal+1+4*(order-1)+1)=0;
        elseif order-(NS-1)^2<=(M-1)
            K(4*Ntotal+1+4*(order-1)+1,4*(order-1)+1)=1;
            K(4*Ntotal+1+4*(order-1)+1,4*((NS-1)^2+M-1-order)+3)=-1;
            C(4*Ntotal+1+4*(order-1)+1)=0;
        end
    elseif V(order).cont(1)==4
        K(4*Ntotal+1+4*(order-1)+1,4*(order-1)+1)=1;
        K(4*Ntotal+1+4*(order-1)+1,4*(order-M)+3)=-1;
        C(4*Ntotal+1+4*(order-1)+1)=0;
    end
    % Face 2
    if V(order).cont(2)==1
        K(4*Ntotal+1+4*(order-1)+2,4*(order-1)+2)=1;
        K(4*Ntotal+1+4*(order-1)+2,4*order+4)=-1;
        C(4*Ntotal+1+4*(order-1)+2)=0;
    elseif V(order).cont(2)==2
        K(4*Ntotal+1+4*(order-1)+2,4*(order-1)+2)=1;
        K(4*Ntotal+1+4*(order-1)+2,4*((NS-1)^2+3*(M-1)/4-order/(NS-1))+1)=-1;
        C(4*Ntotal+1+4*(order-1)+2)=0;
    elseif V(order).cont(2)==4
        K(4*Ntotal+1+4*(order-1)+2,4*(order-1)+2)=1;
        K(4*Ntotal+1+4*(order-1)+2,4*(order-2)+4)=-1;
        C(4*Ntotal+1+4*(order-1)+2)=0;
    elseif V(order).cont(2)==5
        K(4*Ntotal+1+4*(order-1)+2,4*(order-1)+2)=1;
        K(4*Ntotal+1+4*(order-1)+2,4*(order+M-3)+4)=-1;
        C(4*Ntotal+1+4*(order-1)+2)=0;
    end
    % Face 3
    if V(order).cont(3)==1
        K(4*Ntotal+1+4*(order-1)+3,4*(order-1)+3)=1;
        K(4*Ntotal+1+4*(order-1)+3,4*(order-NS)+1)=-1;
        C(4*Ntotal+1+4*(order-1)+3)=0;
    elseif V(order).cont(3)==2
        K(4*Ntotal+1+4*(order-1)+3,4*(order-1)+3)=1;
        K(4*Ntotal+1+4*(order-1)+3,4*((NS-1)^2+M-order-1)+1)=-1;
        C(4*Ntotal+1+4*(order-1)+3)=0;
    elseif V(order).cont(3)==4
        K(4*Ntotal+1+4*(order-1)+3,4*(order-1)+3)=1;
        K(4*Ntotal+1+4*(order-1)+3,4*(order+M-2)+1)=-1;
        C(4*Ntotal+1+4*(order-1)+3)=0;
    end
    % Face 4
    if V(order).cont(4)==1
        K(4*Ntotal+1+4*(order-1)+4,4*(order-1)+4)=1;
        K(4*Ntotal+1+4*(order-1)+4,4*(order-2)+2)=-1;
        C(4*Ntotal+1+4*(order-1)+4)=0;
    elseif V(order).cont(4)==2
        K(4*Ntotal+1+4*(order-1)+4,4*(order-1)+4)=1;
        K(4*Ntotal+1+4*(order-1)+4,4*((NS-1)^2+floor(order/(NS-1)))+1)=-1;
        C(4*Ntotal+1+4*(order-1)+4)=0;
    elseif V(order).cont(4)==4
        K(4*Ntotal+1+4*(order-1)+4,4*(order-1)+4)=1;
        K(4*Ntotal+1+4*(order-1)+4,4*order+2)=-1;
        C(4*Ntotal+1+4*(order-1)+4)=0;
    elseif V(order).cont(4)==5
        K(4*Ntotal+1+4*(order-1)+4,4*(order-1)+4)=1;
        K(4*Ntotal+1+4*(order-1)+4,4*(order-M-3)+2)=-1;
        C(4*Ntotal+1+4*(order-1)+4)=0;
    end 
end

% Solve the linear system to obtain the displacement
K=sparse(K);
U=K\C;

% Assign surface-averaged displacement
T=0;
for order=1:Ntotal
    V(order).u(1)=U(4*(order-1)+1);
    V(order).u(2)=U(4*(order-1)+2);
    V(order).u(3)=U(4*(order-1)+3);
    V(order).u(4)=U(4*(order-1)+4);
    
    V(order).W(5)=(V(order).A(1)*V(order).A(6)-V(order).A(3)*V(order).A(4))/3/(V(order).ab(1,1)*V(order).Length(1)+V(order).ab(2,2)*V(order).Length(2)-V(order).ab(3,1)*V(order).Length(3)-V(order).ab(4,2)*V(order).Length(4))...
        *[V(order).Length(1) V(order).Length(2) V(order).Length(3) V(order).Length(4)]...
        *[-V(order).c5c6(1,1)+V(order).c5c6(1,2);
        -V(order).c5c6(2,1)+V(order).c5c6(2,2);
        -V(order).c5c6(3,1)+V(order).c5c6(3,2);
        -V(order).c5c6(4,1)+V(order).c5c6(4,2)]...
        +1/3/(V(order).ab(1,1)*V(order).Length(1)+V(order).ab(2,2)*V(order).Length(2)-V(order).ab(3,1)*V(order).Length(3)-V(order).ab(4,2)*V(order).Length(4))...
        *[2*V(order).ab(1,1)*V(order).Length(1)+1/2*V(order).ab(2,1)*V(order).Length(2)-V(order).ab(3,1)*V(order).Length(3)+1/2*V(order).ab(4,1)*V(order).Length(4),...
        1/2*V(order).ab(1,2)*V(order).Length(1)+2*V(order).ab(2,2)*V(order).Length(2)+1/2*V(order).ab(3,2)*V(order).Length(3)-V(order).ab(4,2)*V(order).Length(4),...
        V(order).ab(1,1)*V(order).Length(1)-1/2*V(order).ab(2,1)*V(order).Length(2)-2*V(order).ab(3,1)*V(order).Length(3)-1/2*V(order).ab(4,1)*V(order).Length(4),...
        1/2*V(order).ab(1,2)*V(order).Length(1)+V(order).ab(2,2)*V(order).Length(2)-1/2*V(order).ab(3,2)*V(order).Length(3)-2*V(order).ab(4,2)*V(order).Length(4)]...
        *[V(order).u(1);V(order).u(2);V(order).u(3);V(order).u(4)];
    V(order).W(1)=1/2*V(order).u(2)-1/2*V(order).u(4);
    V(order).W(2)=-1/2*V(order).u(1)+1/2*V(order).u(3);
    V(order).W(3)=1/2*V(order).u(2)+1/2*V(order).u(4)-V(order).W(5);
    V(order).W(4)=1/2*V(order).u(1)+1/2*V(order).u(3)-V(order).W(5);
    
    V(order).Jc=[V(order).A(1) V(order).A(4);V(order).A(3) V(order).A(6)];
    dwdxdwdy=inv(V(order).Jc)*[V(order).W(1);V(order).W(2)];
    V(order).strainc=[1/2*(dwdxdwdy(1)-theta*V(order).Cent(2));1/2*(dwdxdwdy(2)+theta*V(order).Cent(1))];
    V(order).stressc=[V(order).G(1)*(dwdxdwdy(1)-theta*V(order).Cent(2));V(order).G(2)*(dwdxdwdy(2)+theta*V(order).Cent(1))];
    T=T+(-V(order).stressc(1)*V(order).Cent(2)+V(order).stressc(2)*V(order).Cent(1))*V(order).Area;
end

fprintf('Angle of twist per unit length: %6.4f rad/m \n',theta)
fprintf('  Numerical Torsional rigidity: %6.4f N.m^2/rad \n',T/theta)
fprintf('                  Total torque: %6.4f N.m \n',T)

%% Plot
% Plot stress
x=zeros(Ntotal,1);
y=zeros(Ntotal,1);
z=zeros(Ntotal,1);
for order=1:Ntotal
    x(order)=V(order).Cent(1);
    y(order)=V(order).Cent(2);
    z(order)=sqrt(V(order).stressc(1)^2+V(order).stressc(2)^2);
end

xv = linspace(min(x), max(x), 200);
yv = linspace(min(y), max(y), 200);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);

figure('Name','Shear Stress resultant')
surf(X, Y, Z/10^6);
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Shear Stress resultant (MPa)')
set(gca,'FontSize',20)
grid on
colormap jet
shading interp

% Plot out-of-plane displacement
x=zeros(Ntotal,1);
y=zeros(Ntotal,1);
z=zeros(Ntotal,1);
for order=1:Ntotal
    x(order)=V(order).Cent(1);
    y(order)=V(order).Cent(2);
    z(order)=V(order).W(5);
end

xv = linspace(min(x), max(x), 200);
yv = linspace(min(y), max(y), 200);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);

figure('Name','Z-Displacement')
surf(X, Y, Z);
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z-Displacement (m)')
set(gca,'FontSize',20)
grid on
colormap jet
shading interp