pro test_diff

    DIR='SDO/'
    get_fits_list,DIR+'files.txt',list
    N=N_ELEMENTS(list)
    for i=0,N-1L do begin
        work=mrdfits(DIR+list[i],0,hd)
        odiffxy5,work,workx,worky
        mwrfits,float(workx),DIR+'diff_'+list[i],/create
        mwrfits,float(worky),DIR+'diff_'+list[i]
    endfor
stop

end
