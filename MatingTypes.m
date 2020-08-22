function [FrB,FrU,FrMt1,FrMt2] = MatingTypes(mu,rsex,r)
    % Returns equilibrium frequencies of B, U, UMt1 and BMt2 alleles
    % mu: mtDNA mutation rate. (0.05)
    % rsex: sex rate. (0.9)
    % r: recombination rate between mating type and uniparentla inheritance regulation loci (0)
    
    M=50; % the number of mitochondria per cell
    n=6; 
    x=1.5; % fitness nonlinearity, x>1
    s=1; % selection strength
    l=0.; % paternal leakage
    delta=0.001; % gamete death rate
    
    %
    P=zeros(M+1,n); % B U UMt1 BMt2 UMt2 BMt1
    P(1,1)=1;

    U=zeros(M+1,M+1);
    for i=0:M
        for j=0:M
            U(i+1,j+1)=binopdf(i-j,M-j,mu);
        end
    end
    w=zeros(M+1,n);
    for i=0:M
        w(i+1,:)=1 - s*(i / M)^x;
    end
    L1=zeros(2*M+1,2*M+1);
    for i=0:2*M
        for j=0:2*M
            L1(i+1,j+1)=hygepdf(i,4*M,2*j,2*M);
        end
    end
    L2=zeros(M+1,2*M+1);
    for i=0:M
        for j=0:2*M
            L2(i+1,j+1)=hygepdf(i,2*M,j,M);
        end
    end
    PSI=zeros(0,1);
    for i=0:2*M
        for j=0:round(l*M)+M
            PSI(i+1,j+1)=binopdf(i,2*M,j/(round(l*M)+M));
        end
    end
    PSIP=zeros(0,1);
    for i=0:l*M
        for j=0:M
            PSIP(i+1,j+1)=hygepdf(i,M,j,round(l*M));
        end
    end
    A=zeros(0,1);
    for i=0:M
        for j=0:M
            A(i+1,j+1)=hygepdf(i,2*M,2*j,M);
        end
    end
    
    for k=1:10000
        if k==100
            P(:,2)=0.05*P(:,1);
            P(:,1)=P(:,1)-0.01*P(:,1);
        end
        if k==500
            P(:,3)=0.1*P(:,2);
            P(:,2)=P(:,2)-0.1*P(:,2);
        end
        if k==1000
            P(:,4)=0.1*P(:,1);
            P(:,1)=P(:,1)-0.1*P(:,1);
        end
               
        P = U*P;
        P = (w.*P) / sum(w.*P,'All');

        fB0=sum(P(:,1)); fU0=sum(P(:,2));
        fUMt1=sum(P(:,3)); fBMt2=sum(P(:,4));
        fUMt2=sum(P(:,5)); fBMt1=sum(P(:,6));
        FB0=fB0; FU0=fU0;
        FUMt1=fUMt1; FBMt2=fBMt2;
        FUMt2=fUMt2; FBMt1=fBMt1;
        dt=0.01;
        for i=1:2000
            sqsum=1-FUMt1^2-FBMt2^2-FUMt2^2-FBMt1^2;
            FB0=FB0+( -delta*FB0-FB0+fB0*(delta+sqsum) )*dt;
            FU0=FU0+( -delta*FU0-FU0+fU0*(delta+sqsum) )*dt;
            FUMt1=FUMt1+( -delta*FUMt1-FUMt1*(1-FUMt1)+fUMt1*(delta+sqsum) )*dt;
            FBMt2=FBMt2+( -delta*FBMt2-FBMt2*(1-FBMt2)+fBMt2*(delta+sqsum) )*dt;
            FUMt2=FUMt2+( -delta*FUMt2-FUMt2*(1-FUMt2)+fUMt2*(delta+sqsum) )*dt;
            FBMt1=FBMt1+( -delta*FBMt1-FBMt1*(1-FBMt1)+fBMt1*(delta+sqsum) )*dt;
        end
        
        if fB0==0; fB0=1; assert(FB0==0); end
        if fU0==0; fU0=1; assert(FU0==0); end
        if fUMt1==0; fUMt1=1; assert(FUMt1==0); end
        if fBMt2==0; fBMt2=1; assert(FBMt2==0); end
        if fUMt2==0; fUMt2=1; assert(FUMt2==0); end
        if fBMt1==0; fBMt1=1; assert(FBMt1==0); end
   
        % zygotes
        B0B0=conv(P(:,1),P(:,1))*FB0*FB0/(fB0*fB0); %
        B0U0=conv(PSIP*P(:,1),P(:,2))*2*FB0*FU0/(fB0*fU0); %
        B0U1=conv(PSIP*P(:,1),P(:,3))*2*FB0*FUMt1/(fB0*fUMt1); %
        B0U2=conv(PSIP*P(:,1),P(:,5))*2*FB0*FUMt2/(fB0*fUMt2); %
        B0B1=conv(P(:,1),P(:,6))*2*FB0*FBMt1/(fB0*fBMt1); %
        B0B2=conv(P(:,1),P(:,4))*2*FB0*FBMt2/(fB0*fBMt2); %
        U0U0=conv(P(:,2),P(:,2))*FU0*FU0/(fU0*fU0); %
        U0U1=conv(P(:,2),P(:,3))*2*FU0*FUMt1/(fU0*fUMt1); %
        U0U2=conv(P(:,2),P(:,5))*2*FU0*FUMt2/(fU0*fUMt2); %
        B1U0=conv(PSIP*P(:,6),P(:,2))*2*FBMt1*FU0/(fBMt1*fU0); %
        B2U0=conv(PSIP*P(:,4),P(:,2))*2*FBMt2*FU0/(fBMt2*fU0); %
        B2U1=conv(PSIP*P(:,4),P(:,3))*2*FBMt2*FUMt1/(fBMt2*fUMt1); %
        B1U2=conv(PSIP*P(:,6),P(:,5))*2*FBMt1*FUMt2/(fBMt1*fUMt2); %
        B1B2=conv(P(:,6),P(:,4))*2*FBMt1*FBMt2/(fBMt1*fBMt2); %
        U1U2=conv(P(:,3),P(:,5))*2*FUMt1*FUMt2/(fUMt1*fUMt2); %     
       
        B0U0=PSI*B0U0; B1U0=PSI*B1U0;
        B0U1=PSI*B0U1; B2U0=PSI*B2U0;
        B0U2=PSI*B0U2; B2U1=PSI*B2U1;
        B1U2=PSI*B1U2;
        
        Z=[B0B0,B0U0,B0U1,B0U2,B0B1,B0B2,U0U0,U0U1,U0U2,B1U0,B2U0,B2U1,B1U2,B1B2,U1U2];
        S = sum(Z,'All');
        Z=Z./S;
        
        % recombination
        Z0 = Z;
        Z(:,3)=(1-r)*Z0(:,3)+r*Z0(:,10);
        Z(:,4)=(1-r)*Z0(:,4)+r*Z0(:,11);
        Z(:,10)=(1-r)*Z0(:,10)+r*Z0(:,3);
        Z(:,11)=(1-r)*Z0(:,11)+r*Z0(:,4);
        Z(:,12)=(1-r)*Z0(:,12)+r*Z0(:,13);
        Z(:,13)=(1-r)*Z0(:,13)+r*Z0(:,12);
             
        SP(:,1)=L2*L1*( Z(:,1)+0.5*(Z(:,2)+Z(:,3)+Z(:,4)+Z(:,5)+Z(:,6)) );
        SP(:,2)=L2*L1*( Z(:,7)+0.5*(Z(:,2)+Z(:,8)+Z(:,9)+Z(:,10)+Z(:,11)) );
        SP(:,3)=L2*L1*0.5*(Z(:,3)+Z(:,8)+Z(:,12)+Z(:,15));
        SP(:,4)=L2*L1*0.5*(Z(:,6)+Z(:,11)+Z(:,12)+Z(:,14));
        SP(:,5)=L2*L1*0.5*(Z(:,4)+Z(:,9)+Z(:,13)+Z(:,15));
        SP(:,6)=L2*L1*0.5*(Z(:,5)+Z(:,10)+Z(:,13)+Z(:,14));
        
        P=A*P*(1-rsex)+SP*rsex;
        
        FrB(k)=sum(P(:,1));
        FrU(k)=sum(P(:,2));
        FrMt1(k)=sum(P(:,3))+sum(P(:,6));
        FrMt2(k)=sum(P(:,4))+sum(P(:,5));
       
    end
end
