function decFactors = computeDecimationFactors(coeff)

if (coeff < 2)
    error('Cannot decimate by %g', coeff)
end

decFactors = [];
testFactors = 13:-1:2;

while (coeff > 13)
    rems = mod(coeff, testFactors);
    [~, idx] = min(rems);
    decFactors = [decFactors, testFactors(idx)]; %#ok<AGROW>
    coeff = coeff / testFactors(idx);
end

coeff = floor(coeff);

if (coeff >= 2)
    decFactors = [decFactors, coeff];
end