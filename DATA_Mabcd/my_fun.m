
function z=my_fun(x,P,R,Rb,PET,S0,G0,TM,M0)


Pwarm=repmat(P(1:12),8,1);
PETwarm=repmat(PET(1:12),8,1);
Rwarm=repmat(R(1:12),8,1);
Rbwarm=repmat(Rb(1:12),8,1);
TMwarm=repmat(TM(1:12),8,1);

[n0,~]=size(P);
n=round((n0-12)*0.5);


P1=P(13:n);
E01=PET(13:n);
Ro=R(13:n);
Rb1=Rb(13:n);
T=TM(13:n);


P1=[Pwarm;P1];
E01=[PETwarm;E01];
Ro=[Rwarm;Ro];
Rb1=[Rbwarm;Rb1];
T=[TMwarm;T];


Rdo1=Ro-Rb1;
[m,~]=size(P1);
P1=P1';
E01=E01';

M=zeros(1,m);
Pt=zeros(1,m);
W=zeros(1,m);
b=zeros(1,m);
Y=zeros(1,m);
S=zeros(1,m);
E=zeros(1,m);
Rd=zeros(1,m);
Re=zeros(1,m);
G=zeros(1,m);
G1=zeros(1,m);
G2=zeros(1,m);
G10=G0;
G20=G0;
Rbs0=zeros(1,m);
IGF=zeros(1,m);
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
    b(1,j)=E01(1,j)+x(2);
    Y(1,j)=(W(1,j)+b(1,j))*0.5/x(1)-((((W(1,j)+b(1,j))*0.5/x(1))^2)-(W(1,j)*b(1,j)/x(1)))^0.5;

    S(1,j)=Y(1,j)*exp(-E01(1,j)/b(1,j));

    E(1,j)=Y(1,j)*(1-exp(-E01(1,j)/b(1,j)));
    
    Rd(1,j)=(1-x(3))*(W(1,j)-Y(1,j));
    Re(1,j)=x(3)*(W(1,j)-Y(1,j));
k1=x(4)*x(5);
k2=(1-x(4))*x(5)/x(4);
    G1(1,j)=(G10+x(4)*Re(1,j))/(1+k1);
    G2(1,j)=(G20+(1-x(4))*Re(1,j))/(1+k2);


G(1,j)=G1(1,j)+G2(1,j);
    Rbs0(1,j)=k1*G1(1,j)+k2*G2(1,j);
    IGF(1,j)=x(7)*G(1,j);
    Rbs(1,j)=Rbs0(1,j)*x(6)+IGF(1,j);
    G10=G1(1,j);
    G20=G2(1,j);
    S0=S(1,j);
    M0=M(1,j);
end



f11=sum((log(Rdo1(97:m)+0.001)-log(Rd(97:m)'+0.001)).^2)/sum((log(Rdo1(97:m)+0.001)-log(mean(Rdo1(97:m))+0.001)).^2);
f12=sum((log(Rb1(97:m)+0.001)-log(Rbs(97:m)'+0.001)).^2)/sum((log(Rb1(97:m)+0.001)-log(mean(Rb1(97:m))+0.001)).^2);

f21=sum((Rdo1(97:m)-Rd(97:m)').^2)/sum((Rdo1(97:m)-mean(Rdo1(97:m))).^2);
f22=sum((Rb1(97:m)-Rbs(97:m)').^2)/sum((Rb1(97:m)-mean(Rb1(97:m))).^2);

f31=1-sum( (Rdo1(97:m)-mean(Rdo1(97:m))).*(Rd(97:m)'-mean(Rd(97:m))) )/( sum( (Rdo1(97:m)-mean(Rdo1(97:m))).^2 )*  sum( (Rd(97:m)-mean(Rd(97:m))).^2 )   ).^0.5;
f32=1-sum( (Rb1(97:m)-mean(Rb1(97:m))).*(Rbs(97:m)'-mean(Rbs(97:m))) )/( sum( (Rb1(97:m)-mean(Rb1(97:m))).^2 )*  sum( (Rbs(97:m)-mean(Rbs(97:m))).^2 )   ).^0.5;

f41=abs(   log( sum(Rd(97:m))/ sum(Rdo1(97:m)))  );
f42=abs(   log( sum(Rbs(97:m))/ sum(Rb1(97:m)))  );

F1=(f11+f21+f31+f41);
F2=(f12+f22+f32+f42);
z=F1+F2;

end
