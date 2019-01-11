function k4 = k4_est(data)

N = length(data);

N4 = N*(N-1)*(N-2)*(N-3);

for i = 1:4
    S(i) = sum( power(data, i) );
end

k4 = (1/N4)*( (N^3 + N^2)*S(4) - 4*(N^2 + N)*S(3)*S(1) - ...
     3*(N^2 - N)*(S(2)^2) + 12*N*S(2)*(S(1)^2) - 6*(S(1)^4));