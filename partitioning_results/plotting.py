import readline
from xml.dom.expatbuilder import parseString
import matplotlib.pyplot as plt
import numpy as np

path = "/home/timo/Coding/mykahypar/partitioning_results/data/"
km1path = path + "km1.txt"
lowerpath = path + "lower_bounds.txt"
targetlowerpath = path + "target_lower_bounds.txt"
upperpath = path + "upper_bounds.txt"
targetupperpath = path + "target_upper_bounds.txt"
infopath = path + "info.txt"
resultspath = path + "other_results.txt"
standarddivspath = path + "standard_divs.txt"
rebalancestepspath = path + "rebalance_steps.txt"

def main():
  plt.figure()
  plt.subplot(211)
  showImbalance()
  plt.subplot(212)
  showKm1()
  plt.show()

def readLines(path):
  with open(path) as file:
    return file.read().splitlines()

def showActualAndTarget(actualpath, actuallabel, targetpath, targetlabel, showfinal):
  actuals = readLines(actualpath)
  num = len(actuals)
  if (num <= 0) : return
  targets = readLines(targetpath)
  if (len(targets) > 0): 
    plt.plot(range(num), floatify(targets), 'r', label=targetlabel)
    if showfinal:
      plt.hlines(y=float(targets[num - 1]), xmin=0, xmax=num - 1 , linewidth=1.5, color='g', label='final target weight')
  plt.plot(range(num), floatify(actuals), 'b', label=actuallabel)

#Expecting the following format for imbalance.txt: one value per line
#Expecting the following format for target_imbalances.txt: one value per line
def showImbalance():
  showActualAndTarget(lowerpath, 'smallest block weight', targetlowerpath, 'target smallest block', False)
  showActualAndTarget(upperpath, 'heaviest block weight', targetupperpath, 'target heaviest block', True)
  
  standard_divs = readLines(standarddivspath)
  plt.plot(floatify(standard_divs), color='y', label='standard deviation')

  rebalancesteps = floatify(readLines(rebalancestepspath))
  if (len(rebalancesteps) != 0):
    plt.vlines(rebalancesteps, 0, floatify(readLines(targetupperpath))[0], 'g', 'dashed', 'rebalance steps')
  plt.legend()
  infoLines = readLines(infopath)
  infoLines.pop()
  infoLines.pop()
  infoLines.pop()
  infoLines.pop()
  plt.title(', '.join(infoLines))


#Expecting the following format: just one value per line, no comma
def showKm1():
  km1lines = floatify(readLines(km1path))
  num = len(km1lines)
  if num == 0: return
  initkm1 = km1lines[0]
  finalkm1 = km1lines[len(km1lines) - 1]
  plt.hlines(y=initkm1, xmin=0, xmax=num - 1 , linewidth=0.75, color='g', label='initial km1: ' + str(initkm1))
  plt.hlines(y=finalkm1, xmin=0, xmax=num - 1 , linewidth=0.75, color='r', label='final km1: ' + str(finalkm1))

  infoLines = readLines(infopath)
  graphFile = getLastWordFromLine(infoLines[0])
  k = getLastWordFromLine(infoLines[1])
  e = getLastWordFromLine(infoLines[2])
  kahyparkm1 = getKm1FromFile(graphFile, k, e)
  if kahyparkm1 != 0:
    plt.hlines(y=kahyparkm1, xmin=0, xmax=num - 1 , linewidth=1, color='y', label='kahypar km1: ' + str(kahyparkm1))

  rebalancesteps = floatify(readLines(rebalancestepspath))
  if (len(rebalancesteps) != 0):
    plt.vlines(rebalancesteps, 0, max(km1lines), 'g', 'dashed', 'rebalance steps')

  plt.title('km1 goal - minimize')
  plt.plot(km1lines, color='b', label='km1')
  plt.legend()

def getLastWordFromLine(line):
  words = line.split(" ")
  return words[len(words) - 1]


def getKm1FromFile(graphFileName, k, e):
  lines = readLines(resultspath)
  for line in lines:
    lineElements = line.split(" ")
    if len(lineElements) != 4:
      print("some Error")
    if lineElements[0] == graphFileName and lineElements[1] == k and lineElements[2] == e:
      return float(lineElements[3])
  return 0

def floatify(arr):
  return [None if v == '' else float(v) for v in arr]

if __name__ == "__main__":
    main()
