pro scale_boulder

restore,'bd.sav'
;0-70 is group 1
;71-105 is group 2

restore,'BD1_final.sav' ;returns A - coefficients

;ALGORITHM
;for i=0,70 do begin ;compute new y for value of x
; bx = EXP(a[1] * x[i])
;  yy = a[0] * bx + a[2]
;newy[i]=yy
;endfor

;VARIABLES
;data.field
;time,startindex,foe,fof1,fxf1,fof2,fxf2,ofoe,ofof1,ofxf1,ofof2,ofxf2,f3o,f3x

;try fof2
foe=make_array(106,value=0.0)
fof2=make_array(106,value=0.0)
fxf2=make_array(106,value=0.0)
fof1=make_array(106,value=0.0)
fxf1=make_array(106,value=0.0)
ofof1=make_array(106,value=0.0)
ofxf1=make_array(106,value=0.0)
ofof2=make_array(106,value=0.0)
ofxf2=make_array(106,value=0.0)
fof3=make_array(106,value=0.0)
fxf3=make_array(106,value=0.0)
f3o=make_array(106,value=0.0)
f3x=make_array(106,value=0.0)
time=data.field01
ifoe=data.field03
ifof1=data.field04
ifxf1=data.field05
ifof2=data.field06
ifxf2=data.field07
iofof1=data.field09
iofxf1=data.field10
iofof2=data.field11
iofxf2=data.field12
ifof3=data.field13
ifxf3=data.field14
timevals=findgen(106)
timevals=timevals*2.0



newy=make_array(106,value=0.0)

for i=0,70 do begin ;do all
;foe
  if (ifoe[i] gt 0) then begin

 bx = EXP(a[1] * ifoe[i])
  yy = a[0] * bx + a[2]
newy[i]=yy
foe[i]=newy[i]
  endif
;fof1
  if (ifof1[i] gt 0) then begin
 bx = EXP(a[1] * ifof1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fof1[i]=yy
endif
;fxf1
  if (ifxf1[i] gt 0) then begin
 bx = EXP(a[1] * ifxf1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fxf1[i]=yy
  endif
;fof2  
if (ifof2[i] gt 0) then begin
 bx = EXP(a[1] * ifof2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fof2[i]=yy
  endif
;fxf2
  if (ifxf2[i] gt 0) then begin
 bx = EXP(a[1] * ifxf2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fxf2[i]=yy
  endif
;ofof1
  if (iofof1[i] gt 0) then begin
 bx = EXP(a[1] * iofof1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofof1[i]=yy
  endif
;ofxf1
  if (iofxf1[i] gt 0) then begin
 bx = EXP(a[1] * iofxf1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofxf1[i]=yy
  endif
;ofof2
  if (iofof2[i] gt 0) then begin
 bx = EXP(a[1] * iofof2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofof2[i]=yy
  endif
;ofxf2
  if (iofxf2[i] gt 0) then begin
 bx = EXP(a[1] * iofxf2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofxf2[i]=yy
  endif
;of3
  if (ifof3[i] gt 0) then begin
 bx = EXP(a[1] * ifof3[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fof3[i]=yy
endif
;xf3
  if (ifxf3[i] gt 0) then begin
 bx = EXP(a[1] * ifxf3[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fxf3[i]=yy
  endif


endfor

restore,'BD2_final.sav'
for i=71,105 do begin
;foe
  if (ifoe[i] gt 0) then begin

 bx = EXP(a[1] * ifoe[i])
  yy = a[0] * bx + a[2]
newy[i]=yy
foe[i]=newy[i]
  endif
;fof1
  if (ifof1[i] gt 0) then begin
 bx = EXP(a[1] * ifof1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fof1[i]=yy
endif
;fxf1
  if (ifxf1[i] gt 0) then begin
 bx = EXP(a[1] * ifxf1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fxf1[i]=yy
  endif
;fof2  
if (ifof2[i] gt 0) then begin
 bx = EXP(a[1] * ifof2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fof2[i]=yy
  endif
;fxf2
  if (ifxf2[i] gt 0) then begin
 bx = EXP(a[1] * ifxf2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fxf2[i]=yy
  endif
;ofof1
  if (iofof1[i] gt 0) then begin
 bx = EXP(a[1] * iofof1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofof1[i]=yy
  endif
;ofxf1
  if (iofxf1[i] gt 0) then begin
 bx = EXP(a[1] * iofxf1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofxf1[i]=yy
  endif
;ofof2
  if (iofof2[i] gt 0) then begin
 bx = EXP(a[1] * iofof2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofof2[i]=yy
  endif
;ofxf2
  if (iofxf2[i] gt 0) then begin
 bx = EXP(a[1] * iofxf2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofxf2[i]=yy
  endif
;of3
  if (ifof3[i] gt 0) then begin
 bx = EXP(a[1] * ifof3[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fof3[i]=yy
endif
;xf3
  if (ifxf3[i] gt 0) then begin
 bx = EXP(a[1] * ifxf3[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fxf3[i]=yy
  endif


endfor


save,file='boulder.sav',timevals,foe,fof1,fxf1,fof2,fxf2,ofof1,ofxf1,ofof2,ofxf2,fof3,fxf3


;;plot,fof1,psym=1,yrange=[2,7]

;print,ifof2
;print,fof2

;lets oplot lusk
restore,'lk.sav'
time=data.field01
ifoe=data.field03
ifof1=data.field04
ifxf1=data.field05
ifof2=data.field06
ifxf2=data.field07
iofof1=data.field09
iofxf1=data.field10
iofof2=data.field11
iofxf2=data.field12
ifof3=data.field13
ifxf3=data.field14

foe[*]=0.0
fof1[*]=0.0
fxf1[*]=0.0
fof2[*]=0.0
fxf2[*]=0.0
ofof1[*]=0.0
ofxf1[*]=0.0
ofof2[*]=0.0
ofxf2[*]=0.0

restore,'LK_final.sav'

for i=0,105 do begin

;foe
  if (ifoe[i] gt 0) then begin

 bx = EXP(a[1] * ifoe[i])
  yy = a[0] * bx + a[2]
newy[i]=yy
foe[i]=newy[i]
  endif
;fof1
  if (ifof1[i] gt 0) then begin
 bx = EXP(a[1] * ifof1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fof1[i]=yy
endif
;fxf1
  if (ifxf1[i] gt 0) then begin
 bx = EXP(a[1] * ifxf1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fxf1[i]=yy
  endif
;fof2  
if (ifof2[i] gt 0) then begin
 bx = EXP(a[1] * ifof2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fof2[i]=yy
  endif
;fxf2
  if (ifxf2[i] gt 0) then begin
 bx = EXP(a[1] * ifxf2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
fxf2[i]=yy
  endif
;ofof1
  if (iofof1[i] gt 0) then begin
 bx = EXP(a[1] * iofof1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofof1[i]=yy
  endif
;ofxf1
  if (iofxf1[i] gt 0) then begin
 bx = EXP(a[1] * iofxf1[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofxf1[i]=yy
  endif
;ofof2
  if (iofof2[i] gt 0) then begin
 bx = EXP(a[1] * iofof2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofof2[i]=yy
  endif
;ofxf2
  if (iofxf2[i] gt 0) then begin
 bx = EXP(a[1] * iofxf2[i])
  yy = a[0] * bx + a[2]
;newy[i]=yy
ofxf2[i]=yy
  endif

endfor

plot,fof1,psym=1,yrange=[2,6],title='Lusk foF1 and fxF1'
;oplot,fof1,psym=2
oplot,fxf1,psym=2

save,file='lusk.sav',timevals,foe,fof1,fxf1,fof2,fxf2,ofof1,ofxf1,ofof2,ofxf2

end
