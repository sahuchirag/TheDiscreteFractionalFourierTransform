t=-3:0.03:3;
k=8;
x=sqrt(2*pi)*t;
hk=hermiteH(k,x);
phik=power(2,1/4)/sqrt(power(2,k)*factorial(k)) *hk.* exp(-pi*(t.^2));
plot(t,phik);
title(['Hermite-Gaussian Function for k = ', num2str(k)]);
grid on;