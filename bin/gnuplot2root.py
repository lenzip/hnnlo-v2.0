#!/usr/bin/env python

from ROOT import *
from array import *
import sys



class MyPlot:
  def __init__(self, line):
    self.title = line.split()[1]
    print "building plot", self.title  
    self.x  = []
    self.y  = []
    self.ey = []

  
  def fillStartX(self,x):
    self.x.append(x)


  def fill(self,x,y,ey):
    self.x.append(x)
    self.y.append(y)
    self.ey.append(ey)
   
  def toTH1F(self):
    xa = array('d', self.x)
    plot = TH1F(self.title, self.title, len(self.x) - 1, xa)

    for i in range(1, len(self.x)):
      plot.SetBinContent(i, self.y[i-1])
      plot.SetBinError(i, self.ey[i-1])
    
    return plot

def lineSplit(line):
  line = line.replace("D", "E").lstrip(" ")
  ls = line.split(" ")
  if len(ls)>=4:
    return (float(ls[0]), float(ls[1]), float(ls[2]), float(ls[3]))
  else:
    return None





filein=open(sys.argv[1])

plots = []
for line in filein.readlines():
  if line.startswith("#"):
    plots.append(MyPlot(line))
    justStarted=True
  else:
    content=lineSplit(line)
    if justStarted:
      plots[-1].fillStartX(content[0])
      justStarted = False
    if content!=None:
      plots[-1].fill(content[1], content[2], content[3])  

outfile=TFile("output.root", "RECREATE")
outfile.cd()
for plot in plots:
  plot.toTH1F().Write()
