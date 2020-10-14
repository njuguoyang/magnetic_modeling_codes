function check_tags,tag_list,data,knot=knot,verbose=verbose
;
   if check_keyword(knot,/knot) then knot=0
   if check_keyword(verbose,/knot) then verbose=0

;     structure?
   ok=(size(data,/type) eq 8)
   if not(ok) then goto,skip
;
;     check for presence of tag
   NT=N_ELEMENTS(tag_list)
   ID=LONARR(NT)
   WORK=TAG_NAMES(data)
   for i=0,NT-1 do ID[i]=where(tag_list[i] eq work) 
   BAD=where(ID eq -1)
   ok= (BAD[0] eq -1) 
   if (verbose and not(ok)) then begin
      text='missing tags in magnetogram structure: '+TAG_LIST[BAD]
      message,text
      goto,skip
   endif
;
;     check sizes
   SZ0=size(data.(ID[0]))
   DIFF=dblarr(NT)
   for i=1,NT-1 do DIFF[i]=total(((size(data.(ID[i]))-SZ0)[0:2])^2,/double)
   bad=where(diff ne 0.d0)
   ok= (BAD[0] eq -1) 
   if (verbose and not(ok)) then begin
      text='size of '+WORK[0]+' does not match '+WORK[ID[bad]] 
      message, text
      goto,skip
   endif
;
;
skip:
   if (knot) then ok=not(ok)
   return,ok

end
