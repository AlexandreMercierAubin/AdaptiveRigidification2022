clear;
F = sym( 'F', [3,3], 'real' );
mu = sym( 'mu' );
lambda = sym( 'lambda' );
volume = sym( 'volume' );

%cheap S of svd decomposition
S = eig(F'*F);
sumS = sum((S-ones(3,1)).^2);
prodS = sum(S)-3;
psi = mu *sumS + 0.5*lambda*prodS*prodS;
psi = simplify(psi);

dpsidF = sym( zeros( 9, 1 ) );
for i = 1:numel(F)
    dpsidF(i) = simplify(-diff( psi, F(i) ) * volume);
end

d2psidF2 = sym( zeros( 9, 9 ) );
for i = 1:numel(F)
    d2psidF2(:,i) = simplify(diff( dpsidF, F(i) ));
end

% asking matlab to make a file for all parts simultaneously allows for a
% better use (elimination) of common sub expressions.  Here, the force is 
% put following Hessian so it is easy to pull off the end of the code gen.
out = [reshape(d2psidF2,[],1); dpsidF; psi ]; 
ccode( out, 'File', 'corotational.c' ); 
% save('corot.c', 'out');
