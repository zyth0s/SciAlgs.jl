
using Plots
using Distributions

function xps_simulation()

   x = [0:1:800]

   x1 = [0:1:800]
   a = 0.3
   b = 0.35
   c = 0.85
   d = 1
   cd1s = a
   cd2s2p = (b*7)+(c*2)
   cd3s3p = (b*7)+(c*8)+(d*2)
   cd3d = (b*9)+(d*18)
   cd4s4p = (b*7)+(c*18)+(d*10)
   cd4d = (b*9)+(d*36)
   cd5s = (b*1)+(c*18)+(d*28)
   cdZ = 48

   Ecd1s=(((cdZ-cd1s)/1)^2)*2
   Ecd2s2p=(((cdZ-cd2s2p)/2)^2)*8
   Ecd3s3p=(((cdZ-cd3s3p)/3)^2)*8
   Ecd3d=(((cdZ-cd3d)/3)^2)*10
   Ecd4s4p=(((cdZ-cd4s4p)/4)^2)*8
   Ecd4d=(((cdZ-cd4d)/4)^2)*10
   Ecd5s=(((cdZ-cd5s)/5)^2)*2
   CdEtot=Ecd1s+Ecd2s2p+Ecd3s3p+Ecd3d+Ecd4s4p+Ecd4d+Ecd5s
   ccd1s=a
   ccd2s2p=(b*7)+(c*2)
   ccd3s3p=(b*7)+(c*8)+(d*2)
   ccd3d=(b*8)+(d*18)
   ccd4s4p=(b*7)+(c*17)+(d*10)
   ccd4d=(b*9)+(d*35)
   ccd5s=(b*1)+(c*18)+(d*27)

   Eccd1s=(((cdZ-ccd1s)/1)^2)*2
   Eccd2s2p=(((cdZ-ccd2s2p)/2)^2)*8
   Eccd3s3p=(((cdZ-ccd3s3p)/3)^2)*8
   Eccd3d=(((cdZ-ccd3d)/3)^2)*9
   Eccd4s4p=(((cdZ-ccd4s4p)/4)^2)*8
   Eccd4d=(((cdZ-ccd4d)/4)^2)*10
   Eccd5s=(((cdZ-ccd5s)/5)^2)*2

   cCdEtot=Eccd1s+Eccd2s2p+Eccd3s3p+Eccd3d+Eccd4s4p+Eccd4d+Eccd5s

   BE=(CdEtot-cCdEtot)*13.6

   ##############SHIELDING CONSTANTS############

   se1s=a
   se2s2p=(b*7)+(c*2)
   se3s3p=(b*7)+(c*8)+(d*2)
   se3d=(b*9)+(d*18)
   se4s4p=(b*5)+(c*18)+(d*10)

   ############################################

   seZ=34

   ########ENERGY OF ELECTRONS PER SHELL#######

   Ese1s=(((seZ-se1s)/1)^2)*2
   Ese2s2p=(((seZ-se2s2p)/2)^2)*8
   Ese3s3p=(((seZ-se3s3p)/3)^2)*8
   Ese3d=(((seZ-se3d)/3)^2)*10
   Ese4s4p=(((seZ-se4s4p)/4)^2)*6

   ############TOTAL ENERGY###################

   seEtot=Ese1s+Ese2s2p+Ese3s3p+Ese3d+Ese4s4p

   ############DEFECT ION CALCULATIONS#########

   sse1s=a
   sse2s2p=(b*7)+(c*2)
   sse3s3p=(b*7)+(c*8)+(d*2)
   sse3d=(b*8)+(d*18)
   sse4s4p=(b*5)+(c*17)+(d*10)

   Esse1s=(((seZ-sse1s)/1)^2)*2
   Esse2s2p=(((seZ-sse2s2p)/2)^2)*8
   Esse3s3p=(((seZ-sse3s3p)/3)^2)*8
   Esse3d=(((seZ-sse3d)/3)^2)*9
   Esse4s4p=(((seZ-sse4s4p)/4)^2)*6

   #########TOTAL ENERGY OF DEFECT ION##########

   sSeEtot=Esse1s+Esse2s2p+Esse3s3p+Esse3d+Esse4s4p

   ############BINDING ENERGY#################

   BE1=(seEtot-sSeEtot)*13.6

   x, BE, BE1
end

x, beCd, beSe = xps_simulation()

dist_Cd = Normal(beCd,10)
dist_Se = Normal(beSe,10)
xps_Cd = pdf.(dist_Cd,x)
xps_Se = pdf.(dist_Se,x)

plot( x,xps_Cd, label = "Cd")
plot!(x,xps_Se, label = "Se", xlabel = "Binding Energy (eV)", ylabel = "Intensity")
