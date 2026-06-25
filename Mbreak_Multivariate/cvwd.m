function retp=cvwd(q,trm)

cv90=(7.316+3.128*q-8.624*trm)*exp(-0.029*q-0.412*trm);
cv95=(9.039+3.318*q-9.969*trm)*exp(-0.030*q-0.327*trm);
cv975=(10.703+3.465*q-11.119*trm)*exp(-0.031*q-0.250*trm);
cv99=(13.189+3.346*q-11.870*trm)*exp(-0.026*q-0.173*trm);
retp=[cv90,cv95,cv975,cv99];
