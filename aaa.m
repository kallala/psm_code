N=64;
x=zeros(N,1);
e=zeros(N,N);
#for k=0:127
#  for j=0:N-1
#    e(k+1,j+1)=exp(-2*i*pi*k*j/N);
#  end
#end


a = -2;
b = 2;

dx = (b-a)/N;
t = a + dx*(0:N-1);
f = exp(-t.*t/0.05);
dfx=-2/0.05*t.*f;
fftx = fft(f);
k = 2*pi/(b-a)*[0:N/2-1, 0, -N/2+1:-1];
dffft = 1i*k.*fftx;
df = ifft(dffft);
plot(t,df)



#e*x-y
#plot(real(y),'r')