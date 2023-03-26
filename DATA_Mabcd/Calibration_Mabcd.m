clc
clear
load("DATA_Mabcd_Sample.mat")


for t=1:10

tStart1 = tic;
    S0=300;
    G0=100;
    M0=50;


    A = [0 0 0 0 0 0 0 1 -1 0 0];
    b = [0];
    nvars=11;
    LB = [0 0.00 0 0 0 0 0.00 -20.0 -10.0 0 0];
    UB = [1 2000 1 1 5 1 20.0 10.00 20.00 1 1];
    options = optimoptions('ga','PopulationSize',200,'MaxGenerations',3000);
    ObjectiveFunction = @(x)my_fun(x,P,R,Rb,PET,S0,G0,T,M0);
    [x,~,~,~,~,~]= ga(ObjectiveFunction,nvars,A,b,[],[],LB,UB,[],options);
    n=NSE_Mabcd(x,P,R,Rb,PET,S0,G0,T,M0);
    tEnd1 = toc(tStart1)/60;
    NN(t,1:4)=n(1:4);
    NN(t,5:15)=x;


    disp(['runing: ',num2str(t),'- ',num2str(tEnd1,'%.2f')])


end



