function spearman_permutation,X,Y,N_TRIALS,seed=seed

   dist=dblarr(N_TRIALS)
   X1=X
   Y1=Y
   R=(R_CORRELATE(X,Y))[0]
   N=N_ELEMENTS(X)
;
   for i=0L,N_TRIALS-1L do begin
;
     r1=randomu(seed,N,/double)
     r2=randomu(seed,N,/double)
     X1=X[sort(r1)]
     Y1=Y[sort(r2)]
     dist[i]=(R_CORRELATE(X1,Y1))[0]

   endfor

return,dist

end
