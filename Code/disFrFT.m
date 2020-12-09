function y = disFrFT(f,a,p)
% Function for DFrFT
% Inputs are the samples of signal (f), order (a), and approximation (p)

N = length(f); 
even = ~rem(N,2);
shft = rem((0:N-1) + fix(N/2),N)+1;
f = f(:);
if (nargin == 2)
    p = N/2;
end
E = dFRFT(N,p);
% DFrFT is calculated and returned to main program
y(shft,1) = E*(exp(-1i*pi/2*a*([0:N-2 N-1+even])).' .*(E'*f(shft)));

function E = dFRFT(N,p)
% Function determines the eigenvectors of the Fourier transform matrix
% It returns NxN eigenvectors of the  Fourier transform matrix
global E_final p_final

if (length(E_final) ~= N || p_final ~= p)
    E = make_E(N,p);
    E_final = E; 
    p_final = p;
else
    E = E_final; 
end

function E = make_E(N,p)
% This function returns sorted eigenvectors and eigenvalues of corresponding vectors

% Construct matrix H, use approximate order p
d2 = [1 -2 1]; 
d_p = 1; 
s = 0; 
st = zeros(1,N);
for k = 1:p/2
    d_p = conv(d2,d_p);
    st([N-k+1:N,1:k+1]) = d_p;
    st(1) = 0;
    temp = [1,1:k-1;1,1:k-1]; 
    temp = temp(:)'./[1:2*k];
    s = s + (-1)^(k-1)*prod(temp)*2*st;        
end
% H = circulant(C_) + diagonal(D_)
col = (0:N-1)'; 
row = (N:-1:1);
idx = col(:,ones(N,1)) + row(ones(N,1),:);
st = [s(N:-1:2).';s(:)];
H = st(idx) + diag(real(fft(s)));

% Construct transformation matrix V
r = floor(N/2);
even = ~rem(N,2);
V1 = (eye(N-1) + flipud(eye(N-1))) / sqrt(2);
V1(N-r:end,N-r:end) = -V1(N-r:end,N-r:end);
if (even) 
    V1(r,r) = 1; 
end
V = eye(N); 
V(2:N,2:N) = V1;

% Compute eigenvectors
VHV = V*H*V';
E = zeros(N);
Ev = VHV(1:r+1,1:r+1);
Od = VHV(r+2:N,r+2:N);
[ve,ee] = eig(Ev);               
[vo,eo] = eig(Od); 
[d,inde] = sort(diag(ee));      
[d,indo] = sort(diag(eo));
ve = ve(:,inde');               
vo = vo(:,indo');
E(1:r+1,1:r+1) = fliplr(ve);     
E(r+2:N,r+2:N) = fliplr(vo);
E = V*E;

% interlacing eigenvectors
ind = [1:r+1;r+2:2*r+2]; 
ind = ind(:);
if (even) 
    ind([N,N+2]) = []; 
else
    ind(N+1) = [];
end
E = E(:,ind');
