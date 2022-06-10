import readline
import matplotlib.pyplot as plt
import numpy as np

# Relevant general data: e, k, flow model, convergence speed?, max iteration
def main():
  path = "/home/timo/Coding/mykahypar/partitioning_results/data/"
  km1path = path + "km1.txt"
  imbpath = path + "imbalances.txt"
  targetpath = path + "target_balances.txt"
  fig, (ax1, ax2) = plt.subplots(2)
  showImbalance(fig, ax1, imbpath, targetpath)
  showKm1(fig, ax2, km1path)
  plt.show()


def readLines(path):
  with open(path) as file:
    return file.read().splitlines()

#Expecting the following format: just one value per line, no comma
def showKm1(fig, ax, km1path):
  lines = readLines(km1path)
  ax.set_title('(k - 1) goal - minimize')
  ax.plot(floatify(lines))

#Expecting the following format for imbalance.txt: target weight value in the first line, then the general data in the second line, afterwards one value per line
#Expecting the following format for target_balances.txt: one value per line
def showImbalance(fig, ax, imbpath, targetpath):
  imbalances = readLines(imbpath)
  targetImbalances = readLines(targetpath)
  targetWeight = imbalances.pop(0)
  generalInfo = imbalances.pop(0)
  ax.plot(floatify(targetImbalances), 'b', label='target imbalances')
  ax.plot(floatify(imbalances), 'g', label='actual imbalances')
  num = len(imbalances)
  ax.hlines(y=targetWeight, xmin=0, xmax=num - 1 , linewidth=1.5, color='r', label='epsilon')
  ax.legend()
  ax.set_title(generalInfo)

def floatify(arr):
  return [None if v == '' else float(v) for v in arr]

if __name__ == "__main__":
    main()
