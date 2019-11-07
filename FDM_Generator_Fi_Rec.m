%% FDM Generator-Fi (Rectangular cross section)
 
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
h=10;                       % The dimension of volume in x-direction
l=10;                       % The dimension of volume in y-direction
Nalpha=30;                  % Amount of subvolumes in x-direction
Nbeta=30;                   % Amount of subvolumes in y-direction
halpha=h/Nalpha;            % The dimension of each subvolume in x-direction
lbeta=l/Nbeta;              % The dimension of each subvolume in y-direction
Ntotal=(Nalpha+2)*(Nbeta+2);% Total number of subvolumes
theta=1;                    % Angle of twist per unit length
G=1;                        % Shear modulus of material
 
% Make structure for each subvolume
V(1:Ntotal)=struct('alphabeta',[],'xy',[],'FDMfi',[],'FDMstress',[],'FDMtauresultant',[],'FDMdisp',[]);
 
%% Subvolume locations
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        V(order).alphabeta=[j,i];   % The index location of each subvolume
        V(order).xy=[0,0];          % Pre-allocate the location of each subvolume
        
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
 
% Fi field
K=zeros(Ntotal,Ntotal);
C=zeros(Ntotal,1);
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        if j==1 || j==Nalpha+2 || i==1 || i==Nbeta+2
            K(order,order)=1;
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
            C(order)=-2*G*theta;
        end
    end
end
 
fi=K\C;
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        V(order).FDMfi=fi(order);
    end
end
 
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        
        % Sigma-xz field
        if j==1 || j==Nalpha+2
            V(order).FDMstress(1)=0;
        elseif i==1
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            V(order).FDMstress(1)=(V(order+Nalpha+2).FDMfi-V(order).FDMfi)/N;
        elseif i==Nbeta+2
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            V(order).FDMstress(1)=(V(order).FDMfi-V(order-Nalpha-2).FDMfi)/S;
        else
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            V(order).FDMstress(1)=(V(order+Nalpha+2).FDMfi-V(order-Nalpha-2).FDMfi)/(S+N);
        end
        
        % Sigma-yz field
        if i==1 || i==Nbeta+2
            V(order).FDMstress(2)=0;
        elseif j==1
            E=V(order+1).xy(1)-V(order).xy(1);
            V(order).FDMstress(2)=-(V(order+1).FDMfi-V(order).FDMfi)/E;
        elseif j==Nalpha+2
            W=V(order).xy(1)-V(order-1).xy(1);
            V(order).FDMstress(2)=-(V(order).FDMfi-V(order-1).FDMfi)/W;
        else
            E=V(order+1).xy(1)-V(order).xy(1);
            W=V(order).xy(1)-V(order-1).xy(1);
            V(order).FDMstress(2)=-(V(order+1).FDMfi-V(order-1).FDMfi)/(E+W);
        end
        
        % Tau resultant field
        V(order).FDMtauresultant=sqrt(V(order).FDMstress(1)^2+V(order).FDMstress(2)^2);
    end
end
 
% Out-of-plane displacement field
k=zeros(Ntotal+1+Nalpha*Nbeta,Ntotal);
c=zeros(Ntotal+1+Nalpha*Nbeta,1);
 
count=0;
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        if i==1
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            k(order,order)=-1/N;
            k(order,order+Nalpha+2)=1/N;
            c(order)=-V(order).xy(1);
        elseif  i==Nbeta+2
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            k(order,order)=1/S;
            k(order,order-Nalpha-2)=-1/S;
            c(order)=-V(order).xy(1);
        elseif j==1
            E=V(order+1).xy(1)-V(order).xy(1);
            k(order,order)=-1/E;
            k(order,order+1)=1/E;
            c(order)=V(order).xy(2);
        elseif j==Nalpha+2
            W=V(order).xy(1)-V(order-1).xy(1);
            k(order,order)=1/W;
            k(order,order-1)=-1/W;
            c(order)=V(order).xy(2);
        else
            count=count+1;
            N=V(order+Nalpha+2).xy(2)-V(order).xy(2);
            E=V(order+1).xy(1)-V(order).xy(1);
            S=V(order).xy(2)-V(order-Nalpha-2).xy(2);
            W=V(order).xy(1)-V(order-1).xy(1);
            k(order,order-1)=-1/(E+W);
            k(order,order+1)=1/(E+W);
            c(order)=V(order).FDMstress(1)/G/theta+V(order).xy(2);
            k(Ntotal+1+count,order-Nalpha-2)=-1/(N+S);
            k(Ntotal+1+count,order+Nalpha+2)=1/(N+S);
            c(Ntotal+1+count)=V(order).FDMstress(2)/G/theta-V(order).xy(1);
        end
    end
end
 
% Fixing condition
% k(Ntotal+1,(Nalpha+2)*((Nbeta+2)/2-1)+(Nalpha+2)/2)=1;
% k(Ntotal+1,(Nalpha+2)*((Nbeta+2)/2-1)+(Nalpha+2)/2+1)=1;
% k(Ntotal+1,(Nalpha+2)*((Nbeta+2)/2)+(Nalpha+2)/2)=1;
% k(Ntotal+1,(Nalpha+2)*((Nbeta+2)/2)+(Nalpha+2)/2+1)=1;
% c(Ntotal+1)=0;

% Constraint 1
% for n=1:Ntotal
%     K((Ntotal+1),n)=1;
% end
% C(Ntotal+1)=0;

% Constraint 2
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        if i==1 || i==Nbeta+2 || j==1 || j==Nalpha+2
             K((Ntotal+1),order)=1;
        end
    end
end
C(Ntotal+1)=0;

Sai=k\c;
 
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        V(order).FDMdisp=Sai(order)*theta;
    end
end
 
%% Comparing Results
Maxdisp=0;
Maxtau=0;
Maxfi=0;
M1=0;
M2=0;
for i=1:Nbeta+2
    for j=1:Nalpha+2
        order=(i-1)*(Nalpha+2)+j;
        if V(order).xy(2)==-l/2 || V(order).xy(2)==l/2 || V(order).xy(1)==-h/2 || V(order).xy(1)==h/2
        else
            % 100 testing points
            if mod(i-1,Nbeta/10)==floor((Nbeta/10)/2)+1 || Nbeta==10
                if mod(j-1,Nalpha/10)==floor((Nalpha/10)/2)+1 || Nalpha==10
                    % Maximum out-of-plane displacement
                    if V(order).FDMdisp>Maxdisp
                        Maxdisp=V(order).FDMdisp;
                    end
                    % Maximum tau resultant
                    if V(order).FDMtauresultant>Maxtau
                        Maxtau=V(order).FDMtauresultant;
                    end
                    % Maximum Prandtl's stress
                    if V(order).FDMfi>Maxfi
                        Maxfi=V(order).FDMfi;
                    end
                end
            end
            % Moment
            M1=M1+(V(order).xy(1)*V(order).FDMstress(2)-V(order).xy(2)*V(order).FDMstress(1))*(halpha*lbeta);
            M2=M2+2*V(order).FDMfi*(halpha*lbeta);
        end
    end
end
  
toc