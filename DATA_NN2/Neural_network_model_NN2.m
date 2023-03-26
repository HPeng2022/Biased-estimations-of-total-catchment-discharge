
clc 
clear 
load('DATA_NN2_Sample.mat')


input=[SL,NDVI,LAT];
output=net(input');
plot(1:80,output,'-o')
ylabel('ω','fontsize',20,'fontname','Times New Roman');%
xlabel('Catchments','fontsize',20,'fontname','Times New Roman');%