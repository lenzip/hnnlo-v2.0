#!/bin/env python

import sys

from ROOT import TFile, TH1F

import math
import json

mass=sys.argv[1]#+"PtJetMin50"
what="incljets"


filece=TFile(mass+"_ce/output.root")
fileup=TFile(mass+"_up/output.root")
filedo=TFile(mass+"_do/output.root")

histo_ce=filece.Get(what)
histo_up=fileup.Get(what)
histo_do=filedo.Get(what)


#histo_ce.Scale(1./histo_ce.Integral())
#histo_up.Scale(1./histo_up.Integral())
#histo_do.Scale(1./histo_do.Integral())
#histo_ce.Scale(1./histo_ce.GetBinContent(1))
#histo_up.Scale(1./histo_up.GetBinContent(1))
#histo_do.Scale(1./histo_do.GetBinContent(1))

f=[]
delta=[]
theta=[]
k=[]
for i in range(4):
  f.append((histo_ce.GetBinContent(i+1) - histo_ce.GetBinContent(i+2))/histo_ce.GetBinContent(1))
  deltaUp = abs(histo_up.GetBinContent(i+1) - histo_ce.GetBinContent(i+1))
  deltaDo = abs(histo_do.GetBinContent(i+1) - histo_ce.GetBinContent(i+1))
  delta.append(max(deltaUp,deltaDo))
  theta.append(delta[-1]/histo_ce.GetBinContent(i+1))
  k.append(math.exp(theta[-1]))
#print f
#print theta
QCDscale_0 = math.pow(k[0],(1./f[0])) 
QCDscale1in_0 = math.pow(k[1], -(f[1]+f[2]+f[3])/f[0])
QCDscale1in_1 = math.pow(k[1],  (f[1]+f[2]+f[3])/f[1])
QCDscale2in_1 = math.pow(k[2], -(f[2]+f[3])/f[1])
QCDscale2in_2 = math.pow(k[2],  (f[2]+f[3])/f[2])
QCDscale3in_2 = math.pow(k[3], -f[3]/f[2])
QCDscale3in_3 = k[3]

#print   ("%.2f\t-\t-\t-" % QCDscale_0)
#print   ("%.2f\t%.2f\t-\t-" % (QCDscale1in_0, QCDscale1in_1))
#print   ("-\t%.2f\t%.2f\t-" % (QCDscale2in_1, QCDscale2in_2))
#print   ("-\t-\t%.2f\t%.2f" % (QCDscale3in_2, QCDscale3in_3))

dictionary={int(mass):{
                      "QCDscale":{"0jet":QCDscale_0, "1jet":1.0, "2jet":1.0, "VBF":1.0},
                      "QCDscale1in": {"0jet":QCDscale1in_0, "1jet":QCDscale1in_1, "2jet":1.0, "VBF":1.0},
                      "QCDscale2in": {"0jet":1.0, "1jet":QCDscale2in_1, "2jet":QCDscale2in_2, "VBF":1.0},
                      "QCDscale3in": {"0jet":1.0, "1jet":1.0, "2jet":QCDscale3in_2, "VBF":QCDscale3in_3},
                      }
            }

#print dictionary
json_string = json.dumps(dictionary, indent=4)
print json_string

ratiodo=histo_do.Clone()
ratiodo.Divide(histo_ce)

ratioup=histo_up.Clone()
ratioup.Divide(histo_ce)

#print "*"*50
#print mass
#print "scale Up:", ratioup.GetBinContent(1),ratioup.GetBinContent(2),ratioup.GetBinContent(3),ratioup.GetBinContent(4)
#print "scale Do:", ratiodo.GetBinContent(1),ratiodo.GetBinContent(2),ratiodo.GetBinContent(3),ratiodo.GetBinContent(4)
#print "*"*50


histo_ce.SetLineColor(1)
histo_up.SetLineColor(2)
histo_do.SetLineColor(4)

#histo_ce.Draw()
#histo_up.Draw("sames")
#histo_do.Draw("sames")

#a=raw_input("press any key..")
