%tesla comparision data
clear
clc

%averaging step (Two lines are
%                               one for efficiency-matched set of data
%                               one for Power     -matched set of data


%%%%%%%%   PW       EFF           Ma     Flow    Pin  Pout   Ratio  quality RPM
%case 1-1
DATA(:,:,1)=[ 1        2         3      4       5       6       7   8   9
              1        2         3      4       5       6       7   8   8];
%case 1-2          
DATA(:,:,2)=[ 1        2         3      4       5       6       7   8  9
              1        2         3      4       5       6       7   8  9];  
%Case 1-3          
DATA(:,:,3)=[ 1        2         3      4       5       6       7   8  9
              1        2         3      4       5       6       7   8  9];  
          
%
%
%
%Case 1-n
n=3; %update to number of cases
DATAbar=zeros(n);

for i=1:n     
DATAbar=mean(DATA(:,:,i));
end
%% experimental data
%each line is for one set
%26 line:cases expected
EXP=        [ 1          2            3      4       5       6 
              1          2            3      4       5       6 
              1          2            3      4       5       6 ];
%%
figure
plot(1:n,DATAbar(:,1),'r*',1:n,EXP(:,1),'b-s')
Flow=[];
xlabel('cases')  % should be updated to different parameters
ylabel('power (W)')
grid on