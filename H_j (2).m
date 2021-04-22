function H = H_j(N2,w)
H = zeros(length(w),2*N2+1);
for i=1:length(w)
    H(i,:) = besselh(-N2:N2,w(i));
end
end