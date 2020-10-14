The IDL procedures in this directory are aimed to track magnetic field lines.

fieldline3d.pro is the main procedure. Please read the head of this procedure to start.

There are two example 3D magnetic fields in "bxyz_titov_demoulin.sav" and "example.sav".
The simplest way to start tracking field lines in IDL is:
IDL>  restore,'example.sav',/ver
IDL>  fieldline3d,bx,by,bz 

Then, you can configure the properties of magnetic field lines using the full feature of iplot in IDL,
for example, changing the line color, style, and thickness etc.
