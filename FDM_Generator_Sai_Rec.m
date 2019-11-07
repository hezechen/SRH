%% FDM Generator-Sai (Rectangular cross section)

%{
Topic: Saint-Venant's semi-inverse torsion problem
Advisors: Marek Pindera and Jose Gomez
Student: Heze Chen and Seth Engel
 
Cross section shape: Rectangular
Cross section material: Isotropic/transverse isotropic with PMC being xy
%}

% Clear all
clc; close all; clearvars;

tic
%% Sub-volume propertites
h=50;                       % The dimension of volume in x-direction
l=2;                        % The dimension of volume in y-direction
Nalpha=550;                   % Amount of subvolumes in x-direction
Nbeta=22;                     % Amount of subvolumes in y-direction
halpha=h/Nalpha;             % The dimension of each subvolume in x-direction
lbeta=l/Nbeta;               % The dimension of each subvolume in y-direction
Ntotal=(Nalpha+2)*(Nbeta+2); % Total number of subvolumes
theta=1;                     % Angle of twist per unit length
G=1;                         % Shear modulus of material

% Make structure for each subvolume
V(1:Ntotal)=struct('alphabeta',[],'xy',[],'FDMdisp',[],'FDMfi',[],'FDMstress',[],'FDMtauresultant',[]);

%% Subvolume locations
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        V(order).alphabeta=[j,i];    % The index location of each subvolume
        V(order).xy=[0,0];           % Pre-allocate the location of each subvolume
        
        if j==1
            V(order).xy(1)=-h/2;
        elseif j==Nalpha+2
            V(order).xy(1)=h/2;
        else
            V(order).xy(1)=-h/2+halpha/2+(j-2)*halpha;
        end
        
        if i==1
            V(order).xy(2)=-l/2;
        elseif i==Nbeta+2
            V(order).xy(2)=l/2;
        else
            V(order).xy(2)=-l/2+lbeta/2+(i-2)*lbeta;
        end
        
    end
end

% Out-of-plane displacement field
K=zeros(Ntotal+1,Ntotal);
C=zeros(Ntotal+1,1);
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        if i==1 && j==1
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            E=V(order+1).xy(1)-V(order).xy(1);
            S=N;
            W=E;
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order+1)=2/W/(E+W)+2/E/(E+W);
            K(order,order+Nalpha+2)=2/S/(N+S)+2/N/(N+S);
            C(order)=2*V(order).xy(2)/W-2*V(order).xy(1)/S;
        elseif  i==Nbeta+2 && j==1
            E=V(order+1).xy(1)-V(order).xy(1);
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            N=S;
            W=E;
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order+1)=2/W/(E+W)+2/E/(E+W);
            K(order,order-Nalpha-2)=2/S/(N+S)+2/N/(N+S);
            C(order)=2*V(order).xy(2)/W+2*V(order).xy(1)/N;
        elseif  j==Nalpha+2 && i==Nbeta+2
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            W=V(order).xy(1)-V(order-1).xy(1);
            E=W;
            N=S;
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order-1)=2/W/(E+W)+2/E/(E+W);
            K(order,order-Nalpha-2)=2/S/(N+S)+2/N/(N+S);
            C(order)=-2*V(order).xy(2)/E+2*V(order).xy(1)/N;
        elseif j==Nalpha+2 && i==1
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            W=V(order).xy(1)-V(order-1).xy(1);
            S=N;
            E=W;
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order-1)=2/W/(E+W)+2/E/(E+W);
            K(order,order+Nalpha+2)=2/S/(N+S)+2/N/(N+S);
            C(order)=-2*V(order).xy(2)/E-2*V(order).xy(1)/S;
        elseif i==Nbeta+2
            E=V(order+1).xy(1)-V(order).xy(1);
            W=V(order).xy(1)-V(order-1).xy(1);
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            N=S;
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order-1)=2/W/(E+W);
            K(order,order+1)=2/E/(E+W);
            K(order,order-Nalpha-2)=2/S/(N+S)+2/N/(N+S);
            C(order)=2*V(order).xy(1)/N;
        elseif i==1
            E=V(order+1).xy(1)-V(order).xy(1);
            W=V(order).xy(1)-V(order-1).xy(1);
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            S=N;
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order-1)=2/W/(E+W);
            K(order,order+1)=2/E/(E+W);
            K(order,order+Nalpha+2)=2/S/(N+S)+2/N/(N+S);
            C(order)=-2*V(order).xy(1)/S;
        elseif j==1
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            E=V(order+1).xy(1)-V(order).xy(1);
            W=E;
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order+1)=2/W/(E+W)+2/E/(E+W);
            K(order,order+Nalpha+2)=2/N/(N+S);
            K(order,order-Nalpha-2)=2/S/(N+S);
            C(order)=2*V(order).xy(2)/W;
        elseif j==Nalpha+2
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            W=V(order).xy(1)-V(order-1).xy(1);
            E=W;
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order-1)=2/W/(E+W)+2/E/(E+W);
            K(order,order+Nalpha+2)=2/N/(N+S);
            K(order,order-Nalpha-2)=2/S/(N+S);
            C(order)=-2*V(order).xy(2)/E;
        else
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            E=V(order+1).xy(1)-V(order).xy(1);
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            W=V(order).xy(1)-V(order-1).xy(1);
            K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
            K(order,order-1)=2/W/(E+W);
            K(order,order+1)=2/E/(E+W);
            K(order,order-Nalpha-2)=2/S/(N+S);
            K(order,order+Nalpha+2)=2/N/(N+S);
            C(order)=0;
        end
    end
end

% % Out-of-plane displacement field (Alternative)
% K=zeros(Ntotal+1,Ntotal);
% C=zeros(Ntotal+1,1);
% for i=1:Nbeta+2
%     for j=1:Nalpha+2
%         order=(i-1)*(Nalpha+2)+j;
%         if i==1
%             N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
%             K(order,order)=-1/N;
%             K(order,order+Nalpha+2)=1/N;
%             C(order)=-V(order).xy(1);
%         elseif  i==Nbeta+2
%             S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
%             K(order,order)=1/S;
%             K(order,order-Nalpha-2)=-1/S;
%             C(order)=-V(order).xy(1);
%         elseif  j==Nalpha+2
%             W=V(order).xy(1)-V(order-1).xy(1);
%             K(order,order)=1/W;
%             K(order,order-1)=-1/W;
%             C(order)=V(order).xy(2);
%         elseif j==1
%             E=V(order+1).xy(1)-V(order).xy(1);
%             K(order,order)=-1/E;
%             K(order,order+1)=1/E;
%             C(order)=V(order).xy(2);
%         else
%             N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
%             E=V(order+1).xy(1)-V(order).xy(1);
%             S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
%             W=V(order).xy(1)-V(order-1).xy(1);
%             K(order,order)=-(2/E/(E+W)+2/W/(E+W)+2/N/(N+S)+2/S/(N+S));
%             K(order,order-1)=2/W/(E+W);
%             K(order,order+1)=2/E/(E+W);
%             K(order,order-Nalpha-2)=2/S/(N+S);
%             K(order,order+Nalpha+2)=2/N/(N+S);
%             C(order)=0;
%         end
%     end
% end

% % Fix the origin
% K(Ntotal+1,(Nalpha+2)*((Nbeta+2)/2-1)+(Nalpha+2)/2)=1;
% K(Ntotal+1,(Nalpha+2)*((Nbeta+2)/2-1)+(Nalpha+2)/2+1)=1;
% K(Ntotal+1,(Nalpha+2)*((Nbeta+2)/2)+(Nalpha+2)/2)=1;
% K(Ntotal+1,(Nalpha+2)*((Nbeta+2)/2)+(Nalpha+2)/2+1)=1;

% Constraint 1
% for n=1:Ntotal
%     K((Ntotal+1),n)=1;
% end
% C(Ntotal+1)=0;

% Constraint 2
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        req=0;
        if i==1
            req=req+1;
        end
        if i==Nbeta+2
            req=req+1;
        end
        if j==1
            req=req+1;
        end
        if j==Nalpha+2
            req=req+1;
        end
        if req==1
            K((Ntotal+1),order)=1;
        end 
    end
end
C(Ntotal+1)=0;

% Calculate the displacements for every point
K=sparse(K);
Sai=K\C;

for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        V(order).FDMdisp=Sai(order)*theta;
    end
end


for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        % Sigma-xy1 field
        if j==1 || j==Nalpha+2
            V(order).FDMstress(1)=0;
        else
            E=V(order+1).xy(1)-V(order).xy(1);
            W=V(order).xy(1)-V(order-1).xy(1);
            V(order).FDMstress(1)=G*theta*((V(order+1).FDMdisp-V(order-1).FDMdisp)/(E+W)-V(order).xy(2));
        end
        % Sigma-yy1 field
        if i==1 || i==Nbeta+2
            V(order).FDMstress(2)=0;
        else
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            V(order).FDMstress(2)=G*theta*((V(order+Nalpha+2).FDMdisp-V(order-Nalpha-2).FDMdisp)/(S+N)+V(order).xy(1));
        end
        % Tau resultant field
        V(order).FDMtauresultant=sqrt(V(order).FDMstress(1)^2+V(order).FDMstress(2)^2);
    end
end

% Fi field
% for i=1:Nbeta+2
%     for j=1:Nalpha+2
%         order=(i-1)*(Nalpha+2)+j;
%         V(order).FDMfi=0;
%         if j==1 || j==Nalpha+2 || i==1 || i==Nbeta+2
%         else
%             for jp=2:j
%                 orderp=(i-1)*(Nalpha+2)+jp;
%                 E=V(orderp+1).xy(1)-V(orderp).xy(1);
%                 W=V(orderp).xy(1)-V(orderp-1).xy(1);
%                 V(order).FDMfi=V(order).FDMfi-V(orderp).FDMstress(2)*(E+W)/2;
%             end
%             V(order).FDMfi=V(order).FDMfi+V(order).FDMstress(2)*E/2;
%         end
%     end
% end

% W, E, S, N
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        V(order).WESN=[0,0,0,0];
        if order>1 && mod(order-1,Nalpha+2)~=0
            V(order).WESN(1)=V(order).xy(1)-V(order-1).xy(1);
        end
        if order+1<=Ntotal && mod(order,Nalpha+2)~=0
            V(order).WESN(2)=V(order+1).xy(1)-V(order).xy(1);
        end
        if order-Nalpha-2>=1
            V(order).WESN(3)=V(order).xy(2)-V(order-Nalpha-2).xy(2);
        end
        if order+Nalpha+2<=Ntotal
            V(order).WESN(4)=V(order+Nalpha+2).xy(2)-V(order).xy(2);
        end
    end
end

% Fi field
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        V(order).FDMfi=0;
        for k1=2:j
            orderp=(i-1)*(Nalpha+2)+k1;
            V(order).FDMfi=V(order).FDMfi-V(orderp).FDMstress(2)*(V(orderp).WESN(1)+V(orderp).WESN(2))/2;
        end
        if k1>=2
            V(order).FDMfi=V(order).FDMfi-V((i-1)*(Nalpha+2)+1).FDMstress(2)*V((i-1)*(Nalpha+2)+1).WESN(2)/2;
            V(order).FDMfi=V(order).FDMfi+V(order).FDMstress(2)*V(order).WESN(2)/2;
        end
        
        %         for k2=2:i
        %             orderp=(k2-1)*(Nalpha+2)+j;
        %             V(order).FDMfi=V(order).FDMfi+V(orderp).FDMstress(1)*(V(orderp).WESN(3)+V(orderp).WESN(4))/2;
        %         end
        %         if k2>=2
        %             V(order).FDMfi=V(order).FDMfi+V(j).FDMstress(1)*V(j).WESN(4)/2;
        %             V(order).FDMfi=V(order).FDMfi-V(order).FDMstress(1)*V(order).WESN(4)/2;
        %         end
    end
end

%% Comparing Results
trueMaxdisp=0;
trueMaxtau=0;
trueMaxfi=0;
M1=0;
M2=0;

Maxdisp=0;
Maxtau=0;
Maxfi=0;

sumfx=0;
sumfy=0;
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        if V(order).FDMdisp>trueMaxdisp
            trueMaxdisp=V(order).FDMdisp;         % True maximum out-of-plane displacement
        end
        if V(order).FDMtauresultant>trueMaxtau
            trueMaxtau=V(order).FDMtauresultant;   % True maximum tau resultant
        end
        if V(order).FDMfi>trueMaxfi
            trueMaxfi=V(order).FDMfi;              % True maximum Prandtl stress
        end
        % Testing points
        if V(order).xy(2)==-l/2 || V(order).xy(2)==l/2 || V(order).xy(1)==-h/2 || V(order).xy(1)==h/2
        else
            if mod(i-1,Nbeta/2)==floor((Nbeta/2)/2)+1 || Nbeta==2
                if mod(j-1,Nalpha/40)==floor((Nalpha/40)/2)+1 || Nalpha==40
                    if V(order).FDMdisp>Maxdisp
                        Maxdisp=V(order).FDMdisp;        % Maximum out-of-plane displacement
                    end
                    if V(order).FDMtauresultant>Maxtau
                        Maxtau=V(order).FDMtauresultant; % Maximum tau resultant
                    end
                    if V(order).FDMfi>Maxfi
                        Maxfi=V(order).FDMfi;            % Maximum Prandtl's stress
                    end
                end
            end
            % Check if the force along x or y is 0
            sumfx=sumfx+V(order).FDMstress(1);
            sumfy=sumfy+V(order).FDMstress(2);
            % Moment
            M1=M1+(V(order).xy(1)*V(order).FDMstress(2)-V(order).xy(2)*V(order).FDMstress(1))*(halpha*lbeta);
            M2=M2+V(order).FDMfi*(halpha*lbeta);
        end
        
    end
end

toc

save('FDMAR25(550x22).mat')