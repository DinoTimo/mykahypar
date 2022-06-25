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
  plt.figure()
  plt.subplot(211)
  showImbalance()
  plt.subplot(212)
  showKm1()
  plt.show()

def readLines(path):
  with open(path) as file:
    return file.read().splitlines()

#Expecting the following format for imbalance.txt: one value per line
#Expecting the following format for target_imbalances.txt: one value per line
def showImbalance():
  imbalances = readLines(imbpath)
  num = len(imbalances)
  if (num <= 0) : return
  plt.plot(range(num), floatify(imbalances), 'b', label='actual imbalances')
  targetImbalances = readLines(targetpath)
  if (len(targetImbalances) > 0): 
    plt.plot(range(num), floatify(targetImbalances), 'r', label='target imbalances')
    plt.hlines(y=float(targetImbalances[num - 1]), xmin=0, xmax=num - 1 , linewidth=1.5, color='g', label='final target weight')
  plt.legend()
  infoLines = readLines(infopath)
  plt.title('       '.join(infoLines))


#Expecting the following format: just one value per line, no comma
def showKm1():
  lines = readLines(km1path)
  num = len(lines)
  if num == 0: return
  initkm1 = float(lines[0])
  finalkm1 = float(lines[len(lines) - 1])
  plt.hlines(y=initkm1, xmin=0, xmax=num - 1 , linewidth=0.75, color='g', label='initial km1: ' + str(initkm1))
  plt.hlines(y=finalkm1, xmin=0, xmax=num - 1 , linewidth=0.75, color='r', label='final km1: ' + str(finalkm1))

  infoLines = readLines(infopath)
  infoArr = infoLines[0].split(" ")
  graphFile = infoArr[len(infoArr) - 1]
  kahyparkm1 = getKm1FromFile(graphFile)
  if kahyparkm1 != 0:
    plt.hlines(y=kahyparkm1, xmin=0, xmax=num - 1 , linewidth=1, color='y', label='kahypar km1: ' + str(kahyparkm1))
    
  plt.title('km1 goal - minimize')
  plt.plot(floatify(lines), color='b', label='km1')
  plt.legend()


def getKm1FromFile(graphFileName):
  lines = readLines(resultspath)
  for line in lines:
    lineElements = line.split(" ")
    if len(lineElements) != 2:
      print("some Error")
    if lineElements[0] == graphFileName:
      return float(lineElements[1])
  return 0

def floatify(arr):
  return [None if v == '' else float(v) for v in arr]

if __name__ == "__main__":
    main()
