#BEGIN PLOT /QCDAWAREJETS/C
LogY=0
NormalizeToIntegral=1
YLabel=$\mathrm{d}\sigma / \mathrm{d}\mathrm{C_2}$ [$\mu\text{b}$]
Title=Three-Point Energy Correlation Function double ratio $C_2^{\alpha}$
LegendXPos=0.74
XLabel=$C_2^\alpha$
#RatioPlot=1
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/CC
YLabel=$\mathrm{d}\sigma / \mathrm{d}\mathrm{C_1}$ [$\mu\text{b}$]
Title=Two-point Energy Correlation Function double ratio $C_1^{\alpha}$
#RatioPlot=0
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/ang
LogY=0
NormalizeToIntegral=1
YLabel=$\mathrm{d}\sigma / \mathrm{d}\lambda$ [$\mu\text{b}$]
XLabel=$\lambda_\alpha^1$
Title=IRC safe angularity $\lambda_\alpha^1$
LogX=1
#RatioPlot=1
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/eta
NormalizeToIntegral=1
LogY=1
LegendXPos=0.74
Title=Pseudorapidity of quark and gluon jets
YLabel=$\mathrm{d}\sigma / \mathrm{d}\eta$
XLabel=$\eta$
#RatioPlot=1
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/mass
NormalizeToIntegral=0
LogY=1
LegendXPos=0.74
Title=Masses of Jets
YLabel=$\mathrm{d}\sigma / \mathrm{d}m$
XLabel=$m$ [$\text{GeV}\!/\!c$]
#RatioPlot=1
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/pT
NormalizeToIntegral=0
LogY=1
LegendXPos=0.74
Title=Transverse momentum of Jets
YLabel=$\mathrm{d}\sigma / \mathrm{d}p_\mathrm{T}$ [$\mu\text{b}/(\text{GeV}\!/\!c)$]
XLabel=$p_\mathrm{T}$ [$\text{GeV}\!/\!c$]
XMin=20
#RatioPlot=1
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/C02_1
CustomLegend=$\alpha = 0.2$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/C05_1
CustomLegend=$\alpha = 0.5$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/C1_1
CustomLegend=$\alpha = 1$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/C2_1
CustomLegend=$\alpha = 2$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/CC02_1
CustomLegend=$\alpha = 0.2$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/CC05_1
CustomLegend=$\alpha = 0.5$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/CC1_1
CustomLegend=$\alpha = 1$
LogX=1
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/CC2_1
CustomLegend=$\alpha = 2$
LogX=1
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/ang02_1
LegendXPos=0.1
CustomLegend=$\alpha = 0.2$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/ang05_1
LegendXPos=0.1
CustomLegend=$\alpha = 0.5$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/ang1_1
LegendXPos=0.74
CustomLegend=$\alpha = 1$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/ang2_1
LegendXPos=0.74
CustomLegend=$\alpha = 2$
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/Delta_R
Title=Distance between FinalState and Parton jets
YLabel=$\mathrm{d}\sigma / \mathrm{d}\Delta\mathrm{R}$
XLabel=$\Delta\mathrm{R}$
XMax=2
#RatioPlot=1
#END PLOT

#BEGIN PLOT /QCDAWAREJETS/Delta_pT
Title=pT difference [FinalState - Parton]
YLabel=$\mathrm{d}\sigma / \mathrm{d}\Delta p_\mathrm{T}$ [$\mu\text{b}/(\text{GeV}\!/\!c)$]
XLabel=$\Delta p_\mathrm{T}$ [$\text{GeV}\!/\!c$]
#RatioPlot=1
#END PLOT




