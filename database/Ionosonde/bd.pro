Pro bd

data=read_csv('Boulder_scaled.csv')

;UT-time,pixel-0,FoE,FoF1,FxF1,FoF2,FxF2,OFoE,OFoF1,OFxF1,OFoF2,OFxF2,F3,F3,
time=data.FIELD01
foe=data.FIELD03
fof1=data.FIELD04
fxf1=data.FIELD05
fof2=data.FIELD06
fxf2=data.FIELD07
ofoe=data.FIELD08
ofof1=data.FIELD09
ofxf1=data.FIELD10
ofof2=data.FIELD11
ofxf2=data.FIELD12
of3o=data.FIELD13
of3x=data.FIELD14

save,file='bd.sav',data
end
