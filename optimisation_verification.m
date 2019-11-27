clc;
clear all;
close all;
%% Initialising
distance = 1:1:48; %distance from emitter
L = 10;

a0=0.1;                                                                     
a1=9*a0;

%% Declaring a and b values of noiseless signal in the channel, then calculating non-centrality parameters
for i=1:length(distance)
    coeff(i)=  sqrt((0.01*0.01)/((distance(i)^2)*(48-distance(i))^2));
    
    a(i)=a0*coeff(i); %a value  times channel coefficient 
    sa2(i)=L*(a(i))^2; 
    sa(i) = sqrt(sa2(i));
    
    b(i)=a1*coeff(i); %sb value 
    sb2(i)= L*(b(i))^2;
    sb(i) = sqrt(sb2(i));
    
    n(i)= ((1/(sb(i)-sa(i)))*(L-0.5)*log(sb(i)/sa(i))+(sb(i)+sa(i))/2)^2;  

end

theoryBER= 0.5-0.5*marcumq(sb,sqrt(n),L)+0.5*(marcumq(sa,sqrt(n),L));

close all
figure
semilogy(distance,theoryBER,'bO-','LineWidth',2);
grid on
legend('AWGN-Theory');
xlabel('distance, m');
ylabel('Bit Error Rate');
title('BER for OOK modulation in AWGN channel');