import readline
import matplotlib.pyplot as plt
import numpy as np

# Relevant general data: e, k, flow model, convergence speed?, max iteration
def main():
  path = "/home/timo/Coding/mykahypar/partitioning_results/data/"
  km1path = path + "km1.txt"
  imbpath = path + "imbalances.txt"
  infopath = path + "info.txt"
  targetpath = path + "target_imbalances.txt"
  fig, (ax1, ax2) = plt.subplots(2)
  showImbalance(fig, ax1, imbpath, targetpath, infopath)
  showKm1(fig, ax2, km1path)
  plt.show()


def readLines(path):
  with open(path) as file:
    return file.read().splitlines()

#Expecting the following format: just one value per line, no comma
def showKm1(fig, ax, km1path):
  lines = readLines(km1path)
  if len(lines) == 0: return
  ax.set_title('(k - 1) goal - minimize')
  ax.plot(floatify(lines))

#Expecting the following format for imbalance.txt: one value per line
#Expecting the following format for target_imbalances.txt: one value per line
def showImbalance(fig, ax, imbpath, targetpath, infopath):
  imbalances = readLines(imbpath)
  targetImbalances = readLines(targetpath)
  num = len(imbalances)
  if (num <= 0) : return
  ax.plot(floatify(imbalances), 'g', label='actual imbalances')
  if (len(targetImbalances) > 0):
    ax.plot(floatify(targetImbalances), 'b', label='target imbalances')
    ax.hlines(y=float(targetImbalances[num - 1]), xmin=0, xmax=num - 1 , linewidth=1.5, color='r', label='final target weight')
  ax.legend()
  infoLines = readLines(infopath)
  ax.set_title(' - '.join(infoLines))



def floatify(arr):
  return [None if v == '' else float(v) for v in arr]

if __name__ == "__main__":
    main()
