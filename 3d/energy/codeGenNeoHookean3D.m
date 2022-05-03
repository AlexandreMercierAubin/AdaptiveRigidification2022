clear;
F = sym( 'F', [3,3], 'real' );
mu = sym( 'mu' );
lambda = sym( 'lambda' );
volume = sym( 'volume' );

J = det(F);

C = F'*F;
Ic = trace(C);

%https://graphics.pixar.com/library/StableElasticity/paper.pdf eq 14
alpha = 1 + mu/lambda - mu/(4*lambda);
Jn = J-alpha;
psi = 0.5 * mu * (Ic-3) + 0.5 * lambda *Jn*Jn - 0.5 * mu * log(Ic + 1);

dpsidF = sym( zeros( 9, 1 ) );
for i = 1:numel(F)
    dpsidF(i) = -diff( psi, F(i) ) * volume;
end

d2psidF2 = sym( zeros( 9, 9 ) );
for i = 1:numel(F)
    d2psidF2(:,i) = diff( dpsidF, F(i) );
end

% asking matlab to make a file for all parts simultaneously allows for a
% better use (elimination) of common sub expressions.  Here, the force is 
% put following Hessian so it is easy to pull off the end of the code gen.
out = [ reshape(d2psidF2,[],1); dpsidF ; psi ]; 
ccode( out, 'File', 'NeoHookean3DPsiGradHess.c' );
