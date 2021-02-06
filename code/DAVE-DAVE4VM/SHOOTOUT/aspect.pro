;+
; NAME:
;  ASPECT
;
; PURPOSE:
;
;  This function calculates and returns the normalized position
;  coordinates necessary to put a plot with a specified aspect ratio
;  into the currently active graphics window. It works on the display
;  output window as well as in a PostScript output window.
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;  Graphics
;
; CALLING SEQUENCE:
;
;  position = ASPECT(aspectRatio)
;
; INPUTS:
;
;  aspectRatio: A floating point value that is the desired aspect
;     ratio (ratio of heigth to width) of the plot in the current
;     graphics output window. If this parameter is missing, an aspect
;     ratio of 1.0 (a square plot) is assumed.
;
; KEYWORD PARAMETERS:
;
;  MARGIN:  The margin around the edges of the plot. The value must be
;     a floating point value between 0.0 and 0.5. It is expressed in
;     normalized coordinate units. The default margin is 0.15.
;
;  WINDOWASPECT: The aspect ratio of the target window. If not provided,
;     the value is obtained from the current graphics window.
;
; OUTPUTS:
;
;  position: A four-element floating array of normalized coordinates.
;     The order of the elements is [x0, y0, x1, y1], similar to the
;     !P.POSITION system variable or the POSITION keyword on any IDL
;     graphic command.
;
; EXAMPLE:
;
;  To create a plot with an aspect ratio of 1:2 and a margin of
;  0.10 around the edge of the output window, do this:
;
;     plotPosition = ASPECT(0.5, Margin=0.10)
;     PLOT, Findgen(11), POSITION=plotPosition
;
;  Notice this can be done in a single IDL command, like this:
;
;     PLOT, Findgen(11), POSITION=ASPECT(0.5, Margin=0.10)
;
; MODIFICATION HISTORY:
;
;  Written by: David Fanning, November 1996.
;       Added better error checking, 18 Feb 1997, DWF.
;       Added WindowAspect keyword. 10 Feb 2000. DWF
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright 2000 Fanning Software Consulting
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################

FUNCTION ASPECT, aspectRatio, MARGIN=margin, WindowAspect=wAspectRatio

ON_ERROR, 1

   ; Check for aspect ratio parameter and possibilities.

IF N_PARAMS() EQ 0 THEN aspectRatio = 1.0

IF aspectRatio EQ 0 THEN BEGIN
   MESSAGE, 'Aspect Ratio of 0. Changing to 1...', /Informational
   aspectRatio = 1.0
ENDIF

s = SIZE(aspectRatio)
IF s(s(0)+1) NE 4 THEN $
   MESSAGE, 'Aspect Ratio is not a FLOAT. Take care...', /Informational

   ; Check for margins.

IF N_ELEMENTS(margin) EQ 0 THEN margin = 0.15

   ; Error checking.

IF margin LT 0 OR margin GE 0.5 THEN $
   MESSAGE, 'The MARGIN keyword value must be between 0.0 and 0.5.'

   ; Calculate the aspect ratio of the current window.

IF N_Elements(wAspectRatio) EQ 0 THEN wAspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE

   ; Calculate normalized positions in window.

IF (aspectRatio LE wAspectRatio) THEN BEGIN
   xstart = margin
   ystart = 0.5 - (0.5 - margin) * (aspectRatio / wAspectRatio)
   xend = 1.0 - margin
   yend = 0.5 + (0.5 - margin) * (aspectRatio / wAspectRatio)
ENDIF ELSE BEGIN
   xstart = 0.5 - (0.5 - margin) * (wAspectRatio / aspectRatio)
   ystart = margin
   xend = 0.5 + (0.5 - margin) * (wAspectRatio / aspectRatio)
   yend = 1.0 - margin
ENDELSE

position = [xstart, ystart, xend, yend]

RETURN, position
END
