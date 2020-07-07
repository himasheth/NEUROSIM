%File network_example.m
%Figure 2.9 Numerical Simulation of network

function final_model

    %declare model parameters
    global SmVi;
    global gNa ;
    global F ;
    global RTF ;
    global Nae ; 
    global Vm;
    global kpump;
    global kmpump;
    global TmaxGLC;
    global KtGLC;
    global kHKPFK;
    global ktATP;
    global nH;
    global Kg;
    global kPGK;
    global N;
    global kPK;
    global kLDH1; %positive
    global kLDH2; %negative
    global TmaxLAC;
    global KtLAC;
    global kCK1; %positive
    global kCK2; %negative
    global C;
    global PScapVi;
    global Ko2;
    global HbOP;
    global nh;
    global Vcap;
    global O2a;
    global GLCa;
    global LACa;
    global F0;
    global A;
    global VvO;
    global Tv;
    global nOP;
    global alpha;
    global qAK;
    global nAero;
    global rc;
    global Fin;
    global vstim;
    global KO2i;
    global KiMito;
    global vMito0;
    global KmMito;
    global rCBF;
    global AMP;
    global NAD;
    global ADP;
    global Fout;
    global dAMP_dATP;
    global alphaf;
    global k1;
    global k2;
    global k3;
    global E0;
    global dHb0;
    
    
    
    %set initial condition
    S0=[15,1.2,0.0057,0.02,0.16,1,0.026,2.2,5,0.0262,7.01,4.56,0.35,0.0237,0.063];
    %declare right-hand-side function (defined below)
    ODEFUN=@netexddt;
    %perform simulation over time-interval [0,800]
    [t,S]=ode15s(ODEFUN, [0,900], S0);
    
    
    %generate plot
    % figure(1);
    % set(gca,'fontsize',14)
    % plot(t, S(:,1), 'k', t, S(:,6), 'k--',t, S(:,11), 'k:','Linewidth', 3)
    % legend('Nai', 'LACi', 'O2c')
    % xlabel('Time (sec)')
    % ylabel('Concentration (mM)')
    
    figure(2);
    set(gca,'fontsize',14)
    plot(t, S(:,2), 'k', t, S(:,3), 'k--',t, S(:,5), 'k:',t, S(:,7), 'g-.',t, S(:,14), 'y-.', t, S(:,13), 'g:',t, S(:,15),'r--',t, S(:,4), 'r', 'Linewidth', 3)
    legend('GLCi', 'GAP', 'PYR', 'NADH', 'Vv','LACc','dHb','PEP')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')
    
    % figure(3);
    % set(gca,'fontsize',14)
    % plot(t, S(:,10),'k-', 'Linewidth', 3)
    % legend('O2i')
    % xlabel('Time (sec)')
    % ylabel('Concentration (mM)')
    
    figure(4);
    set(gca,'fontsize',14)
    plot(t, S(:,2), 'k',t, S(:,6), 'k--',t, S(:,8), 'y-.','Linewidth', 3)
    legend('GLCI', 'LACI','ATP')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')
    
    % figure(5);
    % set(gca,'fontsize',14)
    % plot(t, S(:,12), 'g--',t, S(:,13), 'g:','Linewidth', 3)
    % legend('GLCc', 'LAC')
    % xlabel('Time (sec)')
    % ylabel('Concentration (mM)')
    
    figure(6);
    set(gca,'fontsize',14)
    plot(t, S(:,1), 'r' ,'Linewidth', 3)
    legend('Nai')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')
    
    Vv=S(:,14);
    dHb=S(:,15);
    
    BOLD = (VvO.*(k1*(dHb/dHb0)) + (k2*(1-dHb/dHb0).*(Vv/VvO))+ (k3*(1-(Vv/VvO))));
    
    figure(7);
    set(gca,'fontsize',14)
    plot(t, BOLD, 'r' ,'Linewidth', 3)
    legend('BOLD')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')
    
        
    
    
    
    end
    
    
    
    function dS=netexddt(t,S)
    
    %define right-hand-side of differential equations
    
    Nai =S(1);
    GLCi=S(2); 
    GAP = S(3);
    PEP=S(4);
    PYR=S(5);
    LACi=S(6); 
    NADH=S(7);
    ATP=S(8);
    PCr=S(9);
    O2i=S(10); 
    O2c=S(11);
    GLCc=S(12);
    LACc=S(13);
    Vv=S(14);
    dHb=S(15);
   
    
    % ================================ assign parameter values =======================================
    
    SmVi=9*10^4 ; 
    gNa=0.0039; %%og 0.0039
    F = 9.65*10^4 ;
    RTF = 26.73 ;
    Nae = 150;
    Vm = -70; 
    kpump = 0.29*10^-6 ; 
    kmpump = 0.5; 
    TmaxGLC = 0.0476;
    KtGLC = 9; %%og 9
    kHKPFK = 0.12;
    ktATP = 1;
    nH = 4;
    Kg = 0.05;
    kPGK = 42.6;
    N = 0.1; %%n=0.01 when using mitochondrial respiration %0.1 in figure 4
    kPK = 86.7;
    kLDH1 = 2000; %positive
    kLDH2 = 44.8; %negative
    TmaxLAC = 0.00628;
    KtLAC = 0.5; %%og 0.5
    kCK1 = 3666;
    kCK2 = 20;
    C = 10;
    PScapVi = 1.6;
    Ko2 = 0.0361;
    HbOP = 8.6;
    nh = 2.73;
    Vcap = 0.0055;
    O2a = 8.34;
    
    %variables of interest
    GLCa=3;       %original 4.8
    LACa=8;     %original 0.313
    
    F0=0.012;
    alpha=0.5;
    VvO=0.0237;
    Tv=35;
    nOP = 15;
    A=2.212;
    qAK=0.92;
    nAero=3;
    rc=0.01;
    KO2i = 0.001;
    KiMito = 183.3;
    vMito0=0.0192;
    KmMito = 0.05;
    NAD=N-NADH;
    ADP=(ATP/2)*(-qAK+sqrt(qAK^2+4*qAK*(A/ATP-1)));
    AMP=A-ATP-ADP;
    dAMP_dATP= (-1+((qAK/2)-1/2))*(sqrt(qAK^2+4*qAK*(A/ATP-1))) + ((qAK)^A)/((ATP)*(sqrt(qAK^2+4*qAK*((A/ATP-1)))));
    alphaf = 0.5;
    rCBF = (Fin)/F0;
    
    %BOLD (Blood Oyxgen Level Dependant) variables 
    E0 = 0.4;
    k1 = 7;
    k2 = 2;
    k3 = (2*E0)-0.2;
    dHb0 = S(15);
    
    
    
    if (0<=t) && (t<=360) 
             vstim = 0.23 ;
        else , vstim =0 ; 
    end
    
    
    
    if (5 <= t) && (t <= 600) 
             Fin = ((1+alphaf)*F0);
    elseif (t <= 0) && (t>=605) 
             Fin = vMito0;
    else
            Fin = F0;
    end
    
    
    
    %vMito(t) 
    %define velocity equations
    
    vLeakNa = (SmVi)*(gNa/F)*(RTF*log(Nae/Nai) - Vm);
    vpump = ((SmVi)*kpump)*(ATP*Nai*((1+(ATP/kmpump))^-1)); %NADH+NAD=N %GLCen 
    vGLCm=(TmaxGLC)*(((GLCc)/(GLCc+KtGLC))-((GLCi)/(GLCi+KtGLC))) ;
    vHKPFK= ((kHKPFK)*ATP)*((1+((ATP/ktATP)^nH))^-1)*((GLCi/(GLCi+Kg)));
    vPGK= (kPGK*GAP*((ADP)*(NAD/NADH))+N);
    vPK=(kPK)*(PEP*ADP);
    vLDH= ((kLDH1)*(PYR)*(NADH))-((kLDH2)*(LACi)*(NAD));
    vLACm= (TmaxLAC)*((((LACi)/(LACi+KtLAC))-((LACc)/(LACc+KtLAC))));
    vATPases=0.149;
    vMaxMito=0.025;
    vMito= vMaxMito*(PYR/(KmMito+PYR))*(1/(1+(ATP/(ADP*KiMito))))*(O2i/(KO2i+O2i));
    Cr=C-PCr;
    vCK=((kCK1)*(PCr*ADP))-((kCK2)*(Cr)*(ATP));
    vO2m=((PScapVi)*((Ko2)*((HbOP/O2c)-1)^(-1/nh))-O2i);
    vO2c= (((2)*Fin)/Vcap)*(O2a-O2c);
    vGLCc= (((2)*Fin)/Vcap)*(GLCa-GLCc);
    vLACc= (((2)*Fin)/Vcap)*(LACa-LACc);
    %Fout=(F0)*(((Vv/VvO)^(1/alpha))+((Tv)*(Vv/VvO)^(-1/2))*((1/VvO)*(Fin-(dVv_dt))));
    Fout = (((F0*((Vv/VvO)^(1/alpha)))+(Tv*((Vv/VvO)^(-1/2))*(1/VvO))*Fin))/(1+(Tv*((Vv/VvO)^(-1/2))*(1/VvO)));
    O2cc=((2)*(O2c))-O2a;
    
    
    
        
    dS=[   vLeakNa - (3)*(vpump) + vstim,... 
           vGLCm - vHKPFK,...
           (2*vHKPFK) - vPGK,...
           vPGK - vPK,...
           vPK-vLDH-vMito,...
           vLDH - vLACm,...
           vPGK-vLDH-vMito,...
           (((-2)*(vHKPFK))+vPGK+vPK-vATPases-vpump+((nOP)*(vMito))+vCK)*(1-(dAMP_dATP)),...
           -vCK,...
           vO2m - ((nAero)*(vMito)),...
           vO2c - ((1/rc)*(vO2m)),...
           vGLCc - ((1/rc)*(vGLCm)),...
           vLACc + ((1/rc)*(vLACm)),...
           Fin - Fout,...
           (Fin)*(O2a-O2cc)-((Fout)*(dHb/Vv))]';
    
    end
    