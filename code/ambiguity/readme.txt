
IDL> restore,'mapbxyz20050527_1017sub.sav',/verbose
IDL> current_flip,mapbx,mapby,mapbz
The present scale factor is 2. Do you want to change it (Y/N)?
Input Y/N:  n (ENTER)
;================================================  
;If 'y', you can input a larger factor to enlarge
;================================================
...
Select a subregion (Region 1) to analyze. 
;================================================
;hold left button to drag, middle button to change
;the size of the box, right to select the region
;================================================
...
x1,y1,nx1,ny1:         109         142          75          66
The present scale factor is 5. Do you want to change it (Y/N)?
Input Y/N:  n (ENTER)
Do you want to zoom in (Y/N)?
Input Y/N:  n (ENTER)
% DEFROI: Left button to mark point
% DEFROI: Middle button to erase previous point
% DEFROI: Right button to close region
;================================================
;Select a region to flip the tranverse field in 
;window 1
;================================================
Do you want to analyze another region (Y/N)?
Input Y/N:  y (ENTER)
Do you want to zoom in (Y/N)?
Input Y/N:  y (ENTER)
;================================================
;hold left button to drag, middle button to change
;the size of the box, right to select the region 
;in window 2. One can select a small region along 
;the large current region.
;================================================
The present scale factor is 5. Do you want to change it (Y/N)?
Input Y/N:  n (ENTER)
% DEFROI: Left button to mark point
% DEFROI: Middle button to erase previous point
% DEFROI: Right button to close region
;================================================
;Select a region to flip the tranverse field in 
;window 3. But in this mode, the field vectors 
;change according the direction defined by the user.
;So, in the following click two points in window
;3 to define a direction.
;================================================
Click a START point for defining a reference direction in this region.
xsta,ysta:          76         143
Click an END point for defining a reference direction in this region.
xend,yend:         132          87
Do you want to analyze another region (Y/N)?
Input Y/N:  n (ENTER)
;================================================
;Or, input 'y' to repeat untill all the large currents
;are removed. If 'n', the ambiguity removed data
;are stored in 'ambiguity_removed_data.sav'.
;================================================
