FUNCTION laplace, feld
return,-4*feld+shift(feld,1,0)+shift(feld,-1,0)+shift(feld,0,1)+shift(feld,0,-1)
END