
function nse=NSE_Mabcd(x,P,R,Rb,PET,S0,G0,TM,M0)


Pwarm=repmat(P(1:12),8,1);
PETwarm=repmat(PET(1:12),8,1);
Rwarm=repmat(R(1:12),8,1);
Rbwarm=repmat(Rb(1:12),8,1);
TMwarm=repmat(TM(1:12),8,1);


[n0,~]=size(P);
n=n0 ;
n1=round((n0-12)*0.5)+96;

P1=P(13:n);
PET1=PET(13:n);
Ro=R(13:n);
Rb1=Rb(13:n);
T=TM(13:n);


P1=[Pwarm;P1];
PET1=[PETwarm;PET1];
Ro=[Rwarm;Ro];
Rb1=[Rbwarm;Rb1];
T=[TMwarm;T];



Rdo1=Ro-Rb1;
[m,~]=size(P1);
P1=P1';
PET1=PET1';
SM=zeros(m,1);
SF=zeros(m,1);
G10=G0;
G20=G0;

for j=1:m

    if T(j)<x(8)
        FT=P1(j)*x(10);
        MT=0;

    elseif T(j)>=x(8)&&T(j)<=x(9)
        FT=P1(j)*x(10)*(x(9)-T(j))/(x(9)-x(8));
        MT=(FT+M0)*x(11)*(T(j)-x(8))/(x(9)-x(8)); 

    else
        FT=0;
        MT=M0*x(11);
    end

        M(1,j)=M0+FT-MT;
        Pt(1,j)=P1(1,j)-FT+MT;

    W(1,j)=Pt(1,j)+S0;
    b(1,j)=PET1(1,j)+x(2);
    Y(1,j)=(W(1,j)+b(1,j))*0.5/x(1)-((((W(1,j)+b(1,j))*0.5/x(1))^2)-(W(1,j)*b(1,j)/x(1)))^0.5;

    S(1,j)=Y(1,j)*exp(-PET1(1,j)/b(1,j));

    E(1,j)=Y(1,j)*(1-exp(-PET1(1,j)/b(1,j)));
    Rd(1,j)=(1-x(3))*(W(1,j)-Y(1,j));
    Re(1,j)=x(3)*(W(1,j)-Y(1,j));
k1=x(4)*x(5);
k2=(1-x(4))*x(5)/x(4);
    G1(1,j)=(G10+x(4)*Re(1,j))/(1+k1);
    G2(1,j)=(G20+(1-x(4))*Re(1,j))/(1+k2);


G(1,j)=G1(1,j)+G2(1,j);
    GD(1,j)=k1*G1(1,j)+k2*G2(1,j);
    IGF(1,j)=x(7)*G(1,j);
    Rbs(1,j)=GD(1,j)*x(6)+IGF(1,j);
    G10=G1(1,j);
    G20=G2(1,j);
    S0=S(1,j);
    M0=M(1,j);
end

Rs=Rd+Rbs;
HCI=sum(Ro(97:m))/(sum(Rd(97:m))+sum(GD(97:m)));
Rs=Rs';
Rs1=Rs(97:n1);
Ro1=Ro(97:n1);
Rs2=Rs(n1+1:end);
Ro2=Ro(n1+1:end);



n1=1-sum( (Rs1-Ro1).^2 )/sum( (Ro1-mean(Ro1)).^2 ); 
n2=1-sum( (Rs2-Ro2).^2 )/sum( (Ro2-mean(Ro2)).^2 ); 
n3=1-sum( (Rs(97:end)-Ro(97:end)).^2 )/sum( (Ro(97:end)-mean(Ro(97:end))).^2 ); 
nse=[n1,n2,n3,HCI];



end

