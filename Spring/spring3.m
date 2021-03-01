L=0.1;
z=linspace(L/2,L/2*3);
my0=4*pi*10^(-7);
m=2;
mp=2;
Bdip=my0.*m./(2*pi.*z.^(3));

plot(z,Bdip,'r')
hold on

a=0.01;


Bp=(my0*m/(2*pi*a^2*L))*(((z+(L/2))./sqrt((z+(L/2)).^2+a^2))-...
                        ((z-(L/2))./sqrt((z-(L/2)).^2+a^2)));
Bp2=(my0*mp/(2*pi*a^2*L))*(((z+(L/2))./sqrt((z+(L/2)).^2+a^2))-...
                        ((z-(L/2))./sqrt((z-(L/2)).^2+a^2)));                    

% plot(z,Bp,'y')
% hold on
plot(z,Bp2,'b')


A=L.^(2)*((R.*sqrt(R.^2+L.^2))-R.*sqrt(R.^2+L.^2)-4.*R.*l_2)/...
    ((R.*sqrt(R.^2+L.^2)+R.^2+L.^2).*(R.*sqrt(R.^2+L.^2)+R.^2+L.^2));
