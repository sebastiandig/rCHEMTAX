# %  x    Matrix to be factored as a*b (best if a larger than b)
# %  sdx  Matrix of standard deviations for x
# %  b0   Stabilising value for b, non zero locations, & initial value
# %  sdb  Matrix of standard deviations for b
# %
# %  * is matrix multiply, 
# %  .* & ./ are matrix element by element multiply & divide operations

# % iteration without sd
# a=a .* (x*b') ./ (a* (b*b') +1e-100);
# b=b.*(a'*x+b0)./(a'*(a*b)+b+1e-100);

a <- a * (x %*% t(b)) / (a * (b %*% t(b)) + 1e-100)
b <- b * (t(a) %*% x + b0) / (t(a) %*% (a %*% b) + b + 1e-100)

# % iteration with sd
# w2x=1./sdx.^2;
# w2b=1./sdb.^2;
# a=a.*((x.*w2x)*b')./(((a*b).*w2x)*b');
# b=b.*((a'*(x.*w2x))+b0.*w2b)./(a'*(w2x.*(a*b))+b.*w2b);    

w2x <- 1 / (sdx^2)
w2b <- 1 / (sdb^2)
a   <- a * ((x * w2x) %*% t(b)) / (((a %*% b) * w2x) %*% t(b))
b   <- b *((t(a) %*% (x * w2x)) + b0 * w2b) / (t(a) %*% (w2x* (a %*% b)) + b * w2b);

# % precalc constant walues
# xw=x.*w2x;
# b0w=b0.*w2b;
# a=a.*(xw*b')./(((a*b).*w2x)*b');
# b=b.*((a'*xw)+b0w)./(a'*((a*b).*w2x)+b.*w2b);

xw <- x * w2x
b0w <- b0 * w2b
a <- a * (xw %*% t(b)) / (((a %*% b) * w2x) %*% t(b))
b <- b * ((t(a) %*% xw) + b0w) / (t(a) %*% ((a %*% b) * w2x) + b * w2b)

# % iteration for sd case as loops
# for i=1:ns
# a(i,:)=a(i,:).*(xw(i,:)*b')./(((a(i,:)*b).*w2x(i,:))*b');
# end

for (i in 1:ns) {
  a[i,] <- a[i,] * (xw[i,] %*% t(b)) / (((a[i,] %*% b) * w2x[i,]) %*% t(b));
}

# for j=1:np
# b(:,j)=b(:,j).*((a'*xw(:,j))+b0w(:,j))./  ...
# (a'*((a*b(:,j)).*w2x(:,j))+b(:,j).*w2b(:,j));
# end

for (j in 1:np) {
  b[,j] <- b[,j] * ((t(a) %*% xw[,j]) + b0w[,j]) / 
           (t(a) %*% ((a %*% b[,j]) * w2x[,j]) + b[,j] * w2b[,j]);
}