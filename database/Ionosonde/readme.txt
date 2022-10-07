NOTE:  The height values are not included here - I'll work on packaging them up later.

bd.pro takes Boulder_scaled.csv and creates an IDL sav file called bd.sav with each of the fields:
time=data.FIELD01  ;time from 16:00 (increments of 2 minutes)
foe=data.FIELD03   
fof1=data.FIELD04
fxf1=data.FIELD05
fof2=data.FIELD06
fxf2=data.FIELD07
ofoe=data.FIELD08   ;these are the oblique echoes
ofof1=data.FIELD09
ofxf1=data.FIELD10
ofof2=data.FIELD11
ofxf2=data.FIELD12
of3o=data.FIELD13
of3x=data.FIELD14
NOTE:  The outputs from bd.sav are range values and should correspond to the frequencies seen in the ionograms - they are not yet scaled
NOTE:  The oblique values were redone because I wasn't happy with the method I was using to identify the cusp. The values from boulder_redo are considered better (and those are what I used for publication).  The main difference is the cusp is often hard to see in the ionograms - the first method I used was to identify the highest frequency pixel that I could identify as an echo - in the redone values I extrapolated where the cusp should be when the echoes became too faint.  So in general, the redone values should have equal or higher frequency values.  You can always compare them by looking at the ionogram (I can send the ionograms to you separately when you are ready for them).

scale_boulder.pro takes the range values in bd.sav and boulder_redo.sav and scales them.
NOTE:  scale_boulder.pro also does the same thing for the Lusk data - I'm including the save files so the program works properly, but I'm leaving out lk.pro and Lusk_scaled.csv to avoid confusion (I can give them to you later if you want).
NOTE:  I broke the observations into two sections, each of which have different scaling coefficients (BD1_final.sav and BD2_final.sav).  I can provide the method for computing these coefficients later if you need to know that.

Outputs are boulder.sav and lusk.sav.  If you do a restore and help, this is what the structures look like:

IDL> restore,'boulder.sav'
IDL> help
FOE             FLOAT     = Array[106]
FOF1            FLOAT     = Array[106]
FOF2            FLOAT     = Array[106]
FOF3            FLOAT     = Array[106]
FXF1            FLOAT     = Array[106]
FXF2            FLOAT     = Array[106]
FXF3            FLOAT     = Array[106]
OFOF1           FLOAT     = Array[106]
OFOF2           FLOAT     = Array[106]
OFXF1           FLOAT     = Array[106]
OFXF2           FLOAT     = Array[106]
TIMEVALS        FLOAT     = Array[106]



