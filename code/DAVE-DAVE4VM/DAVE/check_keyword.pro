function check_keyword,key,knot=knot,value=value
;
;    this routine checks for the presence of keyword "key"
;    IT DOES NOT MODIFY "key"
;
  present=not((n_elements(key) eq 0))
  if (n_elements(knot) eq 0) then knot=0
  if present then value=key
  if (knot ne 0) then present=not(present)
 
return, present

end
