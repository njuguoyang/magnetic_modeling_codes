.com poincare.pro
.com box.pro
.com corner.pro
.com xitp.pro
.com field_interp.pro
.com matrix_interp.pro
.com null_newton_raphson.pro
.com cube_index.pro
.com null.pro
restore,'testfield.sav'
cube_index,bx,by,bz
find_null,bx,by,bz,'./null_para.sav','./cube_index.sav'
