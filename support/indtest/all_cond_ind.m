function res = all_cond_ind(X,indtest,pars,m)
% function res = all_cond_ind(X,indtest,pars,m)
%
% Checks all conditional independences with #cond_set <= m
%
% INPUT:
%   X           data matrix (n x p);
%   indtest     favourite indtest (see folder indtest)
%   m           maximal size of conditioning set (optional)
%
% OUTPUT:
%   res         not sure yet
%
% EXAMPLE:
%

if nargin < 4
    m = size(X,2) - 2;
end

p = size(X,2);
allPairs = nchoosek(1:p,2);

for i = 1:size(allPairs,1)
    currentPair = allPairs(i,:);
    remainingVars = setdiff(1:p,currentPair);
    for j = 1:min(m, p - 2)
        allCondsets = nchoosek(remainingVars, j);
        for k = 1:size(allCondsets, 1)
            currentCondset = allCondsets(k,:);
            a = feval(indtest,X(:,currentPair(1)),X(:,currentPair(2)),X(:,currentCondset),pars);
            fprintf('p-val of ind.test %i against %i given %s: %d', currentPair(1), currentPair(2), vector2str(currentCondset), a);
            if a > 0.05
                fprintf(' ****\n');
            else
                fprintf('\n');
            end
        end
    end
end

res = 0;
