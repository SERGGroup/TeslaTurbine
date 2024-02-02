function V=SumRes_Bloss(Torq,M,Pin,Pout,Tin,Tout,Tamb,Hin,Hout,Delta_Torq,Delta_M,Delta_Pin,Delta_Pout,Delta_Tin,Delta_Tout,Delta_Tamb,Delta_Hin,Delta_Hout,x,N_Torq,N_M,N_Pin,N_Pout,N_Tin,N_Tout,N_Tamb,N_Hin,N_Hout)
Torq_tot=0;
M_tot=0;
Pin_tot=0;
Pout_tot=0;
Tin_tot=0;
Tout_tot=0;
Hin_tot=0;
Hout_tot=0;
Tamb_tot=0;
for to=1:N_Torq
    Torq_tot=Torq_tot+(((Torq(to)-x(1))^(2))/((Delta_Torq)^(2)));
end
for m=1:N_M
     M_tot=M_tot+(((M(m)-x(2))^(2))/((Delta_M)^(2)));
end
for pin=1:N_Pin
    Pin_tot=Pin_tot+(((Pin(pin)-x(3))^(2))/((Delta_Pin)^(2)));
end
for pout=1:N_Pout
    Pout_tot=Pout_tot+(((Pout(pout)-x(4))^(2))/((Delta_Pout)^(2)));
end
for tin=1:N_Tin
Tin_tot=Tin_tot+(((Tin(tin)-x(5))^(2))/((Delta_Tin)^(2)));
end
for tout=1:N_Tout
Tout_tot=Tout_tot+(((Tout(tout)-x(6))^(2))/((Delta_Tout)^(2)));
end
for tamb=1:N_Tamb
    Tamb_tot=((Tamb(tamb)-x(7))^(2))/((Delta_Tamb)^(2));
end
for hin=1:N_Hin
    Hin_tot=Hin_tot+(((Hin(hin)-x(10))^(2))/((Delta_Hin)^(2)));
end
for hout=1:N_Hout
    Hout_tot=Hout_tot+(((Hout(hout)-x(11))^(2))/((Delta_Hout)^(2)));
end

V=Torq_tot+M_tot+Pin_tot+Pout_tot+Tin_tot+Tout_tot+Tamb_tot+Hin_tot+Hout_tot;