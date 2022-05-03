clear;
F = sym( 'F', [3,3], 'real' );
mu = sym( 'mu' );
lambda = sym( 'lambda' );
volume = sym( 'volume' );

E = 0.5 *( F' * F - eye(3) );
froE2 = sum( E.*E, 'all');
psi = mu * froE2 + lambda * 0.5 * trace(E)^2;
psi = simplify(psi);

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
ccode( out, 'File', 'STVK3DGradHess.c' );
