function [dir] = indtest_decision(pval_fw, pval_bw, sig1, sig2)

settings;

if (pval_fw > sig1  && pval_bw < sig2)
    dir = DIR_FW;
elseif (pval_bw > sig1  && pval_fw < sig2) 
    dir = DIR_BW;
else
   dir = DIR_UNKNOWN;
end

end