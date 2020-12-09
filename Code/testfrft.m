% DFrFT for Rectangular Pulse

x=-1.0:0.002:1.0; 
y = rectangularPulse(x);
for a=0:0.2:2
    yfunc=disFrFT(y,a,2);
    plot(x,real(yfunc),x,imag(yfunc),x,abs(yfunc)),legend('Real','Imaginary','Magnitude');
    title(['a = ',num2str(a)]);
    pause(1.0);
end
