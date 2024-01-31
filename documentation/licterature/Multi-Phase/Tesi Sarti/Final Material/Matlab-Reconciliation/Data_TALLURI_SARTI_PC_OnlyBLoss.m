sheet=4; %% number of Excel worksheets analyzed%%
%%disp([num2str('CoolProp Fluids:'), ' ', 
count=0;
for f=1:sheet
    Matrix=xlsread('Data_TALLURI.xlsx',f);
    Number=Matrix(:,66);
    Y=3000;
    for y=1:length(Number)
        if Number(y)~=Y && Number(y)>0
            Y=Number(y);
            count=count+1;
            Point_vector(count,1)=Y;
        end
    end
    
end
N=count;

Flow=zeros(N,1); %% Performance indicator
T_exp_su=zeros(N,1);
T_exp_ex=zeros(N,1);
SH=zeros(N,1);
P_exp_su=zeros(N,1);
P_exp_ex=zeros(N,1);
H_exp_su=zeros(N,1);
H_exp_ex=zeros(N,1);
H_exp_ex_is=zeros(N,1);
Beta=zeros(N,1);
TORQ=zeros(N,1);
Pot=zeros(N,1);
Pot_th=zeros(N,1);
Pot_id=zeros(N,1);
Rot_speed=zeros(N,1);
Eta_th=zeros(N,1);
Eta_mech=zeros(N,1);
Point=zeros(N,1);


ctrl_pump=zeros(N,1); %% Setting
T_boiler=zeros(N,1);
fz_fan=zeros(N,1);

std_M=zeros(N,1); %% standard deviation
std_Pin=zeros(N,1);
std_Pout=zeros(N,1);
std_Tin=zeros(N,1);
std_Tout=zeros(N,1);
std_Tamb=zeros(N,1);
std_Torq=zeros(N,1);

Delta_Flow=zeros(N,1);
Delta_T_exp_su=zeros(N,1);
Delta_T_exp_ex=zeros(N,1);
Delta_P_exp_su=zeros(N,1);
Delta_P_exp_ex=zeros(N,1);
Delta_H_exp_su=zeros(N,1);
Delta_H_exp_ex=zeros(N,1);
Delta_TORQ=zeros(N,1);


Flow_mean=zeros(N,1); %% Mean performance indicator
T_exp_su_mean=zeros(N,1);
T_exp_su_mean_K=zeros(N,1);
T_exp_ex_mean=zeros(N,1);
T_exp_ex_mean_K=zeros(N,1);
SH_mean=zeros(N,1);
P_exp_su_mean=zeros(N,1);
P_exp_ex_mean=zeros(N,1);
H_exp_su_mean=zeros(N,1);
H_exp_ex_mean=zeros(N,1);
H_exp_ex_is_mean=zeros(N,1);
Beta_mean=zeros(N,1);
TORQ_mean=zeros(N,1);
Pot_mean=zeros(N,1);
Pot_th_mean=zeros(N,1);
Pot_id_mean=zeros(N,1);
Eta_th_mean=zeros(N,1);
Eta_mech_mean=zeros(N,1);


Flow_weight=zeros(N,1);   %% valutation of the weight of the reconciliation method
T_exp_su_weight=zeros(N,1);
T_exp_ex_weight=zeros(N,1);
P_exp_su_weight=zeros(N,1);
P_exp_ex_weight=zeros(N,1);
Torq_weight=zeros(N,1);
H_exp_su_weight=zeros(N,1);
H_exp_ex_weight=zeros(N,1);

Bearing_Loss=zeros(N,1);

constr_viol=zeros(N,1);
constr_viol_mean=zeros(N,1);
weights_mean=zeros(N,1);
N_meas_variables=8;

V=0;
for s=1:sheet %% Reading each worksheet%%
fprintf('Worksheet number %d\n',s);
Dati=xlsread('Data_TALLURI.xlsx',s);
test_number=Dati(:,66);
Index=3000;
N_index=0;
for test=1:length(test_number) %% search data in steady state condition
    if test_number(test)~=Index && test_number(test)>0 %% during the tests, points with the same test_number value(different from zero) correspond to points with the same boundaries
        Index=test_number(test); %% if test_number=0 means we weren't in steady state conditions 
        N_index=N_index+1;
        Start_vector(N_index)=test; %% each value of the vector is the Start index for a condition
    end
    if test_number(test)==Index && test_number(test+1)==0
        Stop_vector(N_index)=test; %% each value of the vector is the Stop index for a condition
    end
end
for n=1:N_index
    Start=Start_vector(n);
    Stop=Stop_vector(n);
M=Dati(Start:Stop,49)/1000; %% kg/s
Pin=Dati(Start:Stop,35)*100000; %%Pa because CoolProp use Pascal
Pout=Dati(Start:Stop,36)*100000;
Tin=Dati(Start:Stop,6)+273.16; %% K because CoolProp use Kelvin
Tout=Dati(Start:Stop,8)+273.16;
Torq=Dati(Start:Stop,51);
rpm=Dati(Start:Stop,70);
rad_s=Dati(Start:Stop,70)*2*pi/60;
Hin=zeros(length(M),1);
Hout=zeros(length(M),1);
Hout_s=zeros(length(M),1);
Hin_s=zeros(length(M),1);
Tamb=(15+273.16)*ones(length(M),1);

r_sensor=0.03; %% radius of pipe(m)--->ask
A_sensor=pi*(r_sensor^(2));

for h=1:length(M)
Hin(h)=py.CoolProp.CoolProp.PropsSI('H','P',Pin(h),'T',Tin(h),'R1233zd(E)'); %%J/kg
Hout(h)=py.CoolProp.CoolProp.PropsSI('H','P',Pout(h),'T',Tout(h),'R1233zd(E)');
end

H_tot=0.024; %% total expander's height (m)
D_ex=0.217;  %% Stator External Diamter(m)
D_ex_rot=0.216; %%Rotor External Diameter(m)
D_int_rot=0.055; %% Rotor internal Diameter(m)
d_2=D_ex;
d_3=D_int_rot;
r_2=d_2/2;
r_3=d_3/2;
N_noz=4; %number of nozzle 
TW=0.001; %throat width(m)
n_diskA=30; %% number of disk A(m)
n_diskB=30; %% number of disk B(m)
h_diskA=0.001; %% height of disk A(m)
h_diskB=0.0008; %% height of disk B(m)
n_ch=n_diskA+n_diskB; %% number of channels
t=H_tot-n_diskA*h_diskA-n_diskB*h_diskB; %total height for the flow between disks 
A_in_rot=2*pi*r_3*H_tot;
A_out_noz=TW*H_tot;
Ar=pi*D_ex*H_tot+pi*((D_ex/2)^(2))*2; %%Heat exchange area
eps=1-(N_noz*A_out_noz/A_in_rot);

T_max=400; %%K
M_max=0.9; %%kg/s
Torq_max=100; %%Nm
Pin_max=3000000; %% Pa
Pout_max=600000;
FS_T=0.5;
FS_M=0.0002;
FS_Torq=0.0025;
FS_P=0.003;
Delta_M=FS_M*M_max; %%standard deviation of each instruments%%
Delta_Torq=FS_Torq*Torq_max;
Delta_T=FS_T;
Delta_Pin_sens=Pin_max*FS_P;
Delta_Pout_sens=Pout_max*FS_P;

Max_Hin=0;
Max_Hout=0;
Delta=zeros(length(M),2);
for j=1:length(M)
    for k=1:length(M)
        Delta(k,1)=abs(Hin(j)-Hin(k));
        Delta(k,2)=abs(Hout(j)-Hout(k));
    end
    if max(Delta(:,1))>Max_Hin
    Delta_Hin=max(Delta(:,1))/2;
    end
    if max(Delta(:,2))>Max_Hout
    Delta_Hout=max(Delta(:,2))/2;
    end
    Max_Hin=max(Delta(:,1));
    Max_Hout=max(Delta(:,2));
end
%%Mean value of measurements%%
M_mean=mean(M);
Pin_mean=mean(Pin);
Pout_mean=mean(Pout);
Tout_mean=mean(Tout);
Tin_mean=mean(Tin);
Torq_mean=mean(Torq); %% 0.33= medium offset
rpm_mean=mean(rpm);
rad_s_mean=rpm_mean*2*pi/60;
Hin_mean=mean(Hin);
Hout_mean=mean(Hout);
Hout_mean_s=mean(Hout_s);
Tamb_mean=mean(Tamb);

%%standard deviation
std_M(n+V)=std(M);
if std_M(n+V)>Delta_M
    Delta_M=std_M(n+V);
end

std_Pin(n+V)=std(Pin);
if std_Pin(n+V)>Delta_Pin_sens
    Delta_Pin=std_Pin(n+V);
else
    Delta_Pin=Delta_Pin_sens;
end
 
std_Pout(n+V)=std(Pout);
if std_Pout(n+V)>Delta_Pout_sens
    Delta_Pout=std_Pout(n+V);
else
    Delta_Pout=Delta_Pout_sens;
end

std_Tin(n+V)=std(Tin);
if std_Tin(n+V)>Delta_T
    Delta_Tin=std_Tin(n+V);
else
    Delta_Tin=Delta_T;
end

std_Tout(n+V)=std(Tout);
if std_Tout(n+V)>Delta_T
    Delta_Tout=std_Tout(n+V);
else
    Delta_Tout=Delta_T;
end

std_Tamb(n+V)=std(Tamb);
if std_Tamb(n+V)>Delta_T
    Delta_Tamb=std_Tamb(n+V);
else
    Delta_Tamb=Delta_T;
end

std_Torq(n+V)=std(Torq);
if std_Torq(n+V)>Delta_Torq
    Delta_Torq=std_Torq(n+V);
end
std_Hin=std(Hin);
if std_Hin<Delta_Hin
    Delta_Hin=std_Hin;
end
std_Hout=std(Hout);
if std_Hout<Delta_Hout
    Delta_Hout=std_Hout;
end

%%Bearing losses%%
Corr=xlsread('Bearing losses.xlsx');
r=0.0175;
F=250;
f=0.15;
rpm_bloss=Corr(:,1);
cost=Corr(:,2);
Bloss=0;
for i=1:length(cost)
 if rpm_bloss(i)==rpm_mean
        Bloss=cost(i)+r*f*F*rad_s_mean;
 end
end

%%elimination of outliers
out_m=0;
for o_m=1:length(M)
    if M(o_m)>M_mean+Delta_M || M(o_m)<M_mean-Delta_M 
        M(o_m)=0;
        out_m=out_m+1;
    end
end
N_M_in=length(M)-out_m;
M_in=zeros(N_M_in,1);
in_M=0;
for z_m=1:length(M)
    if M(z_m)~=0
        in_M=in_M+1;
        M_in(in_M)=M(z_m);
    end
end
std_Min=std(M_in);

out_torq=0;
for o_t=1:length(Torq)
    if Torq(o_t)>Torq_mean+Delta_Torq || Torq(o_t)<Torq_mean-Delta_Torq
        Torq(o_m)=0;
        out_torq=out_torq+1;
    end
end
N_Torq_in=length(Torq)-out_torq;
Torq_in=zeros(N_Torq_in,1);
in_Torq=0;
for z_t=1:length(Torq)
    if Torq(z_t)~=0
        in_Torq=in_Torq+1;
        Torq_in(in_Torq)=Torq(z_t);
    end
end
std_Torq_in=std(Torq_in);

out_Pin=0;
for o_pin=1:length(Pin)
    if Pin(o_pin)>Pin_mean+Delta_Pin || Pin(o_pin)<Pin_mean-Delta_Pin
        Pin(o_pin)=0;
        out_Pin=out_Pin+1;
    end
end
N_Pin_in=length(Pin)-out_Pin;
Pin_in=zeros(N_Pin_in,1);
in_Pin=0;
for z_Pin=1:length(Pin)
    if Pin(z_Pin)~=0
        in_Pin=in_Pin+1;
        Pin_in(in_Pin)=Pin(z_Pin);
    end
end
std_Pin_in=std(Pin_in);

out_Pout=0;
for o_pout=1:length(Pout)
    if Pout(o_pout)>Pout_mean+Delta_Pout || Pout(o_pout)<Pout_mean-Delta_Pout
        Pout(o_pout)=0;
        out_Pout=out_Pout+1;
    end
end
N_Pout_in=length(Pout)-out_Pout;
Pout_in=zeros(N_Pout_in,1);
in_Pout=0;
for z_Pout=1:length(Pout)
    if Pout(z_Pout)~=0
        in_Pout=in_Pout+1;
        Pout_in(in_Pout)=Pout(z_Pout);
    end
end
std_Pout_in=std(Pout_in);

out_Tin=0;
for o_tin=1:length(Tin)
    if Tin(o_tin)>Tin_mean+Delta_Tin || Tin(o_tin)<Tin_mean-Delta_Tin
        Tin(o_tin)=0;
        out_Tin=out_Tin+1;
    end
end
N_Tin_in=length(Tin)-out_Tin;
Tin_in=zeros(N_Tin_in,1);
in_Tin=0;
for z_Tin=1:length(Tin)
    if Tin(z_Tin)~=0
        in_Tin=in_Tin+1;
        Tin_in(in_Tin)=Tin(z_Tin);
    end
end
std_Tin_in=std(Tin_in);

out_Tout=0;
for o_tout=1:length(Tout)
    if Tout(o_tout)>Tout_mean+Delta_Tout || Tout(o_tout)<Tout_mean-Delta_Tout
        Tout(o_tout)=0;
        out_Tout=out_Tout+1;
    end
end
N_Tout_in=length(Tout)-out_Tout;
Tout_in=zeros(N_Tout_in,1);
in_Tout=0;
for z_Tout=1:length(Tout)
    if Tout(z_Tout)~=0
        in_Tout=in_Tout+1;
        Tout_in(in_Tout)=Tout(z_Tout);
    end
end
std_Tout_in=std(Tout_in);

out_Tamb=0;
for o_tamb=1:length(Tamb)
    if Tamb(o_tamb)>Tamb_mean+Delta_Tamb || Tamb(o_tamb)<Tamb_mean-Delta_Tamb
        Tamb(o_tamb)=0;
        out_Tamb=out_Tamb+1;
    end
end
N_Tamb_in=length(Tamb)-out_Tamb;
Tamb_in=zeros(N_Tamb_in,1);
in_Tamb=0;
for z_Tamb=1:length(Tamb)
    if Tamb(z_Tamb)~=0
        in_Tamb=in_Tamb+1;
        Tamb_in(in_Tamb)=Tamb(z_Tamb);
    end
end
std_Tamb_in=std(Tamb_in);

out_Hin=0;
for o_hin=1:length(Hin)
    if Hin(o_hin)>Hin_mean+Delta_Hin || Hin(o_hin)<Hin_mean-Delta_Hin
        Hin(o_hin)=0;
        out_Hin=out_Hin+1;
    end
end
N_Hin_in=length(Hin)-out_Hin;
Hin_in=zeros(N_Hin_in,1);
in_Hin=0;
for z_Hin=1:length(Hin)
    if Hin(z_Hin)~=0
        in_Hin=in_Hin+1;
        Hin_in(in_Hin)=Hin(z_Hin);
    end
end
std_Hin_in=std(Hin_in);

out_Hout=0;
for o_hout=1:length(Hout)
    if Hout(o_hout)>Hout_mean+Delta_Hout || Hout(o_hout)<Hout_mean-Delta_Hout
        Hout(o_hout)=0;
        out_Hout=out_Hout+1;
    end
end
N_Hout_in=length(Hout)-out_Hout;
Hout_in=zeros(N_Hout_in,1);
in_Hout=0;
for z_Hout=1:length(Hout)
    if Hout(z_Hout)~=0
        in_Hout=in_Hout+1;
        Hout_in(in_Hout)=Hout(z_Hout);
    end
end
std_Hout_in=std(Hout_in);



%%Reconciliation method-Find the minimum%%
x0=[Torq_mean;M_mean;Pin_mean;Pout_mean;Tin_mean;Tout_mean;Tamb_mean;Bloss;Ar;Hin_mean;Hout_mean;rpm_mean]; %%First guess values%%
lb=[Torq_mean-Delta_Torq;M_mean-Delta_M;Pin_mean-Delta_Pin;Pout_mean-Delta_Pout;Tin_mean-Delta_Tin;Tout_mean-Delta_Tout;Tamb_mean-Delta_Tamb;Bloss;Ar;Hin_mean-Delta_Hin;Hout_mean-Delta_Hout;rpm_mean];
ub=[Torq_mean+Delta_Torq;M_mean+Delta_M;Pin_mean+Delta_Pin;Pout_mean+Delta_Pout;Tin_mean+Delta_Tin;Tout_mean+Delta_Tout;Tamb_mean+Delta_Tamb;Bloss;Ar;Hin_mean+Delta_Hin;Hout_mean+Delta_Hout;rpm_mean];
obj=@(x)SumRes_Bloss(Torq_in,M_in,Pin_in,Pout_in,Tin_in,Tout_in,Tamb_in,Hin_in,Hout_in,Delta_Torq,Delta_M,Delta_Pin,Delta_Pout,Delta_Tin,Delta_Tout,Delta_Tamb,Delta_Hin,Delta_Hout,x,length(Torq_in),length(M_in),length(Pin_in),length(Pout_in),length(Tin_in),length(Tout_in),length(Tamb_in),length(Hin_in),length(Hout_in));
nonlcon=@bound_Bloss;
A=[];
b=[];
Aeq=[];
beq=[];
Tol_bound =0.1; %%tolerance for constraits%%
Max_Iter=1000;
options = optimset('Algorithm','interior-point','TolCon',Tol_bound,'MaxIter',Max_Iter);
[x,fval,exitflag,output]=fmincon(obj,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%%[x,fval,exitflag,output]=fmincon(obj,x0,A,b,Aeq,beq,lb,ub); no constraints
s_in=py.CoolProp.CoolProp.PropsSI('S','P',x(3),'T',x(5),'R1233zd(E)');
s_in_mean=py.CoolProp.CoolProp.PropsSI('S','P',Pin_mean,'T',Tin_mean,'R1233zd(E)');
Hout_is=py.CoolProp.CoolProp.PropsSI('H','P',x(4),'S',s_in,'R1233zd(E)');
Hout_is_mean=py.CoolProp.CoolProp.PropsSI('H','P',Pout_mean,'S',s_in_mean,'R1233zd(E)');
T_sat=py.CoolProp.CoolProp.PropsSI('T','P',x(3),'Q',1,'R1233zd(E)');
rho_out=py.CoolProp.CoolProp.PropsSI('D','P',x(4),'T',x(6),'R1233zd(E)');
v_out=x(2)/(rho_out*A_sensor);
Hout_s_is=Hout_is-((1/2)*v_out^(2)); %static isoentropic enthalpy
DH=x(10)-x(11);
rho_in=py.CoolProp.CoolProp.PropsSI('D','P',x(3),'T',x(5),'R1233zd(E)');
v_in=x(2)/(rho_in*A_sensor);

Bearing_Loss(n+V)=Bloss;

Point(n+V,1)=Point_vector(n+V,1);  %% Perfomance indicator

Flow(n+V)=x(2); 
Flow_mean(n+V)=M_mean;
Flow_weight(n+V)=abs(Flow(n+V)-Flow_mean(n+V))/Delta_M;
Delta_Flow(n+V)=Delta_M;

Rot_speed(n+V)=x(12);

Pot(n+V)=(x(1)+0.33)*x(12)*2*pi/60;
Pot_mean(n+V)=(Torq_mean+0.33)*x(12)*2*pi/60;

T_exp_su(n+V)=x(5)-273.16;
T_exp_su_mean(n+V)=Tin_mean-273.16;
T_exp_su_weight(n+V)=abs(T_exp_su(n+V)-T_exp_su_mean(n+V))/Delta_Tin;
Delta_T_exp_su(n+V)=Delta_Tin;

T_exp_ex(n+V)=x(6)-273.16;
T_exp_ex_mean(n+V)=Tout_mean-273.16;
T_exp_ex_weight(n+V)=abs(T_exp_ex(n+V)-T_exp_ex_mean(n+V))/Delta_Tout;
Delta_T_exp_ex(n+V)=Delta_Tout;

SH(n+V)=x(5)-T_sat;
SH_mean(n+V)=Tin_mean-T_sat;

P_exp_su(n+V)=x(3);
P_exp_su_mean(n+V)=Pin_mean;
P_exp_su_weight(n+V)=abs(P_exp_su(n+V)-P_exp_su_mean(n+V))/Delta_Pin;
Delta_P_exp_su(n+V)=Delta_Pin;

P_exp_ex(n+V)=x(4);
P_exp_ex_mean(n+V)=Pout_mean;
P_exp_ex_weight(n+V)=abs(P_exp_ex(n+V)-P_exp_ex_mean(n+V))/Delta_Pout;
Delta_P_exp_ex(n+V)=Delta_Pout;

H_exp_su(n+V)=x(10);
H_exp_su_mean(n+V)=Hin_mean;
H_exp_su_weight(n+V)=abs(H_exp_su(n+V)-H_exp_su_mean(n+V))/Delta_Hin;
Delta_H_exp_su(n+V)=Delta_Hin;

H_exp_ex(n+V)=x(11);
H_exp_ex_mean(n+V)=Hout_mean;
H_exp_ex_weight(n+V)=abs(H_exp_ex(n+V)-H_exp_ex_mean(n+V))/Delta_Hout;
Delta_H_exp_ex(n+V)=Delta_Hout;

H_exp_ex_is(n+V)=Hout_is;
H_exp_ex_is_mean(n+V)=Hout_is_mean;

Beta(n+V)=x(3)/x(4);
Beta_mean(n+V)=Pin_mean/Pout_mean;

TORQ(n+V)=x(1);
TORQ_mean(n+V)=Torq_mean;
Torq_weight(n+V)=abs(TORQ(n+V)-TORQ_mean(n+V))/Delta_Torq;
Delta_TORQ(n+V)=Delta_Torq;

Pot_th(n+V)=x(2)*(x(10)-x(11));
Pot_th_mean(n+V)=M_mean*(Hin_mean-Hout_mean);

Pot_id(n+V)=x(2)*(x(10)-Hout_is);
Pot_id_mean(n+V)=M_mean*(Hin_mean-Hout_is_mean);


Eta_th(n+V)=Pot_th(n+V)/Pot_id(n+V);
Eta_th_mean(n+V)=Pot_th_mean(n+V)/Pot_id_mean(n+V);

Eta_mech(n+V)=Pot(n+V)/Pot_id(n+V);
Eta_mech_mean(n+V)=Pot_mean(n+V)/Pot_id_mean(n+V);

Eta_tot_tot=Pot_th(n+V)/(x(10)-Hout_is);
Eta_tot_st=Pot_th(n+V)/(x(10)-Hout_s_is);

ctrl_pump(n+V)=mean(Dati(Start:Stop,66));
T_boiler(n+V)=mean(Dati(Start:Stop,67));
fz_fan(n+V)=mean(Dati(Start:Stop,68));

constr_viol(n+V)=output.constrviolation;
constr_viol_mean(n+V)=Pot_th_mean(n+V)-Pot_mean(n+V)-Bearing_Loss(n+V);
weights_mean(n+V)=(Flow_weight(n+V)+T_exp_su_weight(n+V)+T_exp_ex_weight(n+V)+P_exp_su_weight(n+V)+P_exp_ex_weight(n+V)+H_exp_su_weight(n+V)+H_exp_ex_weight(n+V)+Torq_weight(n+V))/N_meas_variables;


fprintf('fval=%8.4f\n',fval);
output
fprintf('Variables Values\nP=%5.3f Watt\nHin=%6.3f\nHout=%6.3f\nHout_i_s=%6.3f\nDelta_H=%6.3f\nm=%8.6f kg/s\nSH=%f\nPin=%6.3f Pa\nPout=%6.3fPa\nTin=%6.2f K\nTout=%6.2f K\nrpm=%6.3f\nBearinglosses=%6.3f\nBeta=%f\nEta_th=%f\nEta_th_mean=%f\nEta_mech=%f\n',Pot(n+V),x(10),x(11),Hout_is,DH,x(2),SH(n+V),x(3),x(4),x(5),x(6),x(12),x(8),Beta(n+V),Eta_th(n+V),Eta_th_mean(n+V),Eta_mech(n+V))
end
V=V+N_index;
end

Table=table(Point,Flow,std_M,T_exp_su,std_Tin,T_exp_ex,std_Tout,SH,P_exp_su,std_Pin,P_exp_ex,std_Pout,H_exp_su,H_exp_ex,H_exp_ex_is,Beta,TORQ,std_Torq,Rot_speed,Pot,Pot_th,Pot_id,Eta_th,Eta_mech,ctrl_pump,T_boiler,fz_fan);
Table_mean=table(Point,Flow_mean,std_M,T_exp_su_mean,std_Tin,T_exp_ex_mean,std_Tout,SH_mean,P_exp_su_mean,std_Pin,P_exp_ex_mean,std_Pout,H_exp_su_mean,H_exp_ex_mean,H_exp_ex_is_mean,Beta_mean,TORQ_mean,std_Torq,Rot_speed,Pot_mean,Pot_th_mean,Pot_id_mean,Eta_th_mean,Eta_mech_mean,ctrl_pump,T_boiler,fz_fan);
Table_weights=table(Point,Flow,Flow_mean,Delta_Flow,Flow_weight,T_exp_su,T_exp_su_mean,Delta_T_exp_su,T_exp_su_weight,T_exp_ex,T_exp_ex_mean,Delta_T_exp_ex,T_exp_ex_weight,P_exp_su,P_exp_su_mean,Delta_P_exp_su,P_exp_su_weight,P_exp_ex,P_exp_ex_mean,Delta_P_exp_ex,P_exp_ex_weight,H_exp_su,H_exp_su_mean,Delta_H_exp_su,H_exp_su_weight,H_exp_ex,H_exp_ex_mean,Delta_H_exp_ex,H_exp_ex_weight,TORQ,TORQ_mean,Delta_TORQ,Torq_weight);
Table_valutation=table(Point,weights_mean,constr_viol,constr_viol_mean,Eta_th,Eta_mech);


n_row=0;
for g=1:length(Point)
    if constr_viol(g)<0.1 && weights_mean(g)<0.3 && Bearing_Loss(g)~=0
        n_row=n_row+1;
        Table_good(n_row,:)=Table(g,:);
    end   
end

Table(50,:) = []; %% just because of tests problems 
Table(50,:) = [];
Table_mean(50,:) = [];
Table_mean(50,:) = [];
Table_weights(50,:) = [];
Table_weights(50,:) = [];
Table_valutation(50,:) = [];
Table_valutation(50,:) = [];
 
writetable(Table,'Results_TALLURI_OnlyBloss.xlsx','Sheet',1);
writetable(Table_mean,'Results_TALLURI_OnlyBloss.xlsx','Sheet',2);
writetable(Table_weights,'Results_TALLURI_OnlyBloss.xlsx','Sheet',3);
writetable(Table_valutation,'Results_TALLURI_OnlyBloss.xlsx','Sheet',4);
writetable(Table_good,'Results_TALLURI_OnlyBloss.xlsx','Sheet',5);


figure;
scatter(Beta,Flow);
hold on
scatter(Beta_mean,Flow_mean);
title('Beta-Flow');

figure;
scatter(Beta,SH);
hold on
scatter(Beta_mean,SH_mean);
title('Beta-SH');

figure;
scatter(Beta,Rot_speed);
hold on
scatter(Beta_mean,Rot_speed);
title('Beta-rpm');

figure;
scatter(Beta,P_exp_su);
hold on
scatter(Beta_mean,P_exp_su_mean);
title('Beta-P_i_n');

figure;
scatter(Beta,Pot_th);
hold on
scatter(Beta_mean,Pot_th_mean);
title('Beta-Power_t_h');

figure;
scatter(Beta,Pot_id);
hold on
scatter(Beta_mean,Pot_id_mean);
title('Beta-Power_i_d');

figure;
scatter(Beta,Eta_th);
hold on
scatter(Beta_mean,Eta_th_mean);
title('Beta-Eta_t_h');

figure;
scatter(Beta,Eta_mech);
hold on
scatter(Beta_mean,Eta_mech_mean);
title('Beta-Eta_m_e_c_h');

figure;
scatter(Flow,Pot_th);
hold on
scatter(Flow_mean,Pot_th_mean);
title('Power_t_h-Flow');

figure;
scatter(Table_valutation.Point,Table_valutation.constr_viol);
hold on
scatter(Table_valutation.Point,Table_valutation.constr_viol_mean);
title('Constraints violation')

figure;
scatter(Table_valutation.Point,Table_valutation.weights_mean);
title('Weights')

figure;
scatter(Table_valutation.constr_viol,Table_valutation.weights_mean)
title('Constr_viol-weights')

