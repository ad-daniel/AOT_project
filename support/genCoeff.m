function coeff = genCoeff(n,p)
% generates AR coefficients ensuring a stable VAR(p) process. Given
%    X_t = A X(t-1) + B X(t-2) + eps_t 
% Reduce VAR(p) to VAR(1) process
%    M = [ A , B ; I , 0 ]
% Ensure that eig(M) are in unit circle for stability


%% reduce VAR(p) to VAR(1) while ensuring stability
% generate components of matrix M
M_row = [];
for pp = 1:p
   % generate random coeffs
   phi{pp} = rand(n,n)-0.5*ones(n,n);   
   M_row = cat(2,M_row,phi{pp});
end

% generate matrix M
if(p == 1)
   M = M_row;
else
   I_row = [eye( (p-1)*n ), zeros((p-1)*n,n)];
   M = [M_row; I_row];
end

% find biggest eigenvalue
if (n == 1)
   a = 0.95*rand / max(abs(eig(M))); % max eigenvalue will be 0.95*rand
   % NB: we add rand otherwise it would always generate the same
else
   a = 0.95 / max(abs(eig(M))); % max eigenvalue will be 0.95
end

% scale back elements of M accordingly to ensure stability, example:
% M = [A B; eye(n) zeros(n)];
% a = 0.9 / max(abs(eig(M)));  NB: max eigenvalue of M will be 0.9
% M = [a*A a^2*B; eye(n) zeros(n)];
% generate components of matrix M
M_row = [];
for pp = 1:p
   M_row = cat(2,M_row,phi{pp}.*(a^pp)); % scale AR coefficient 
end
% generate matrix M
if(p == 1)
   M = M_row;
else
   I_row = [eye( (p-1)*n ), zeros((p-1)*n,n)];
   M = [M_row; I_row];
end

% Recheck if stable for safety after scaling
eigVal = abs(eig(M));

LI = eigVal >= 1 - 1000*eps; % check if any eig outside unit circle
if( nnz(LI) )
   error('Should be stable, but created nonstable!')
end

% All is good, return coeffs
coeff = zeros(n,n,p);
for pp = 1:p
   coeff(:,:,pp) = phi{pp}.*(a^pp);
end