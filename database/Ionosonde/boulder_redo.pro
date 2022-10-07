pro boulder_redo


data=read_csv('Boulder_scaled.csv')

;UT-time,pixel-0,FoE,FoF1,FxF1,FoF2,FxF2,OFoE,OFoF1,OFxF1,OFoF2,OFxF2,F3,F3,

time=data.FIELD01
ofof1=data.FIELD02
ofxf1=data.FIELD03
ofof2=data.FIELD04
ofxf2=data.FIELD05

save,file='boulder_redo.sav',data

end
