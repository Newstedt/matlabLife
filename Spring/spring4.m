K42=6.92*10^(-4);
K4=2607587399;
K5=1668855935;

U=[120,140,160,180,200,220,240,260,280,300];
U2=[120,140,180,200,220,240,260,280,300];
I5=[0.840,1.018,1.112,1.201,1.278,1.344,1.403,1.479,1.533,1.590];
I4=[1.074,1.294,1.418,1.544,1.625,1.718,1.771,1.878,1.940,2.020];
I42=[1.074,1.294,1.544,1.625,1.718,1.771,1.878,1.940,2.020];
epsi1=K4.*U./(I4.^2);
epsi2=K5.*U./(I5.^2);
% epsi2=K2.*U./(I2.^2)

Iavgsquare=((I5+I4)./2).^2;

Isquare4=I4.^2;
Isquare5=I5.^2;

Isquare42=I42.^2;

UK4=U.*K4;
UK5=U.*K5;
UK42=U2.*K4;

UKavg=((UK4+UK5)./2);

Ilong = [Isquare4 Isquare5];
UKlong = [UK4 UK5];

Ilong2 = [Isquare42 Isquare5];

UKlong2 = [UK42 UK5];
