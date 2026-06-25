function sq=sqrm(M)

%compute the square root matrix of a symmetric positive definte matrix@

[ve,va] = eig(M);
dd=va;
sq=ve*sqrt(dd)*ve';
