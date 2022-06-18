import readline
import matplotlib.pyplot as plt
import numpy as np

path = "/home/timo/Coding/mykahypar/partitioning_results/data/"
km1path = path + "km1.txt"
imbpath = path + "imbalances.txt"
infopath = path + "info.txt"
targetpath = path + "target_imbalances.txt"
resultspath = path + "other_results.txt"

def main():
  fig, (ax1, ax2) = plt.subplots(2)
  showImbalance(fig, ax1)
  showKm1(fig, ax2)
  plt.show()

def readLines(path):
  with open(path) as file:
    return file.read().splitlines()

#Expecting the following format: just one value per line, no comma
def showKm1(fig, ax):
  lines = readLines(km1path)
  num = len(lines)
  if num == 0: return
  initkm1 = float(lines[0])
  finalkm1 = float(lines[len(lines) - 1])
  ax.hlines(y=initkm1, xmin=0, xmax=num - 1 , linewidth=1.5, color='g', label='initial km1: ' + str(initkm1))
  ax.hlines(y=finalkm1, xmin=0, xmax=num - 1 , linewidth=1.5, color='r', label='final km1: ' + str(finalkm1))

  infoLines = readLines(infopath)
  infoArr = infoLines[0].split(" ")
  graphFile = infoArr[len(infoArr) - 1]
  kahyparkm1 = getKm1FromFile(graphFile)
  if kahyparkm1 != 0:
    ax.hlines(y=kahyparkm1, xmin=0, xmax=num - 1 , linewidth=1.5, color='y', label='kahypar km1: ' + str(kahyparkm1))
    
  ax.set_title('km1 goal - minimize')
  ax.plot(floatify(lines), color='b', label='km1')
  ax.legend()


def getKm1FromFile(graphFileName):
  lines = readLines(resultspath)
  for line in lines:
    lineElements = line.split(" ")
    if len(lineElements) != 2:
      print("some Error")
    if lineElements[0] == graphFileName:
      return float(lineElements[1])
  return 0

#Expecting the following format for imbalance.txt: one value per line
#Expecting the following format for target_imbalances.txt: one value per line
def showImbalance(fig, ax):
  imbalances = readLines(imbpath)
  targetImbalances = readLines(targetpath)
  num = len(imbalances)
  if (num <= 0) : return
  ax.plot(floatify(imbalances), 'b', label='actual imbalances')
  if (len(targetImbalances) > 0):
    ax.plot(floatify(targetImbalances), 'r', label='target imbalances')
    ax.hlines(y=float(targetImbalances[num - 1]), xmin=0, xmax=num - 1 , linewidth=1.5, color='g', label='final target weight')
  ax.legend()
  infoLines = readLines(infopath)
  ax.set_title('       '.join(infoLines))

def floatify(arr):
  return [None if v == '' else float(v) for v in arr]

if __name__ == "__main__":
    main()
