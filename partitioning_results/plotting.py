import matplotlib.pyplot as plt
import pandas as pd

path = "/home/timo/Coding/mykahypar/partitioning_results/data/"
infopath = path + "info.txt"
resultspath = path + "other_results.txt"

def main():
  df = pd.read_csv('/home/timo/Coding/test.csv', delimiter=',')
  num_nodes = df['num_nodes'].to_list()
  target_imbalances = []
  try:
    target_imbalances = df['target_imbalance'].to_list()
  except:
    target_imbalances = []
  imbalances = df['imbalance'].to_list()
  km1s = df['km1'].to_list()
  rebalance_steps = df['rebalance_step'].to_list()
  rebalance_steps = list(filter(lambda step: step >= 0, rebalance_steps))
  plt.figure()
  
  plt.subplot(211)
  showImbalance(num_nodes, target_imbalances, imbalances, rebalance_steps)

  plt.subplot(212)
  showKm1(num_nodes, km1s, rebalance_steps)
  
  plt.show()

def readLines(path):
  try:
    with open(path) as file:
      return file.read().splitlines()
  except:
    return []

def showImbalance(num_nodes, target_imbalances, imbalances, rebalance_steps):
  highest_point = 0
  lowest_point = 0

  if len(target_imbalances) > 0:
    lowest_point = target_imbalances[len(target_imbalances) - 1]
    highest_point = target_imbalances[0]
    plt.plot(num_nodes, target_imbalances, 'r', label="target heaviest block")
    plt.hlines(y=lowest_point, xmin=num_nodes[0], xmax=num_nodes[len(num_nodes) - 1] - 1, linewidth=1.5, color='y', label='final target weight')
  else:
    lowest_point = imbalances[len(imbalances) - 1]
    highest_point = imbalances[0]

  plt.plot(num_nodes, imbalances, 'b', label="heaviest block")

  if (len(rebalance_steps) != 0):
    plt.vlines(rebalance_steps, lowest_point, highest_point, 'g', 'dashed', 'rebalance steps')
  plt.legend()
  infoLines = readLines(infopath)
  if (len(infoLines) > 3):
    infoLines.pop()
    infoLines.pop()
    infoLines.pop()
  plt.title(', '.join(infoLines))


#Expecting the following format: just one value per line, no comma
def showKm1(num_nodes, km1s, rebalance_steps):
  km1lines = km1s
  initkm1 = km1lines[0]
  finalkm1 = km1lines[len(km1lines) - 1]
  plt.hlines(y=initkm1, xmin=num_nodes[0], xmax=num_nodes[len(num_nodes) - 1], linewidth=0.75, color='g', label='initial km1: ' + str(initkm1))
  plt.hlines(y=finalkm1, xmin=num_nodes[0], xmax=num_nodes[len(num_nodes) - 1], linewidth=0.75, color='r', label='final km1: ' + str(finalkm1))


  infoLines = readLines(infopath)
  graphFile = getLastWordFromLine(infoLines[0])
  k = getLastWordFromLine(infoLines[1])
  e = getLastWordFromLine(infoLines[2])
  kahyparkm1 = getKm1FromFile(graphFile, k, e)
  if kahyparkm1 != 0:
    plt.hlines(y=kahyparkm1, xmin=num_nodes[0], xmax=num_nodes[len(num_nodes) - 1], linewidth=1, color='y', label='kahypar km1: ' + str(kahyparkm1))

  if (len(rebalance_steps) != 0):
    plt.vlines(rebalance_steps, max(km1s), min(km1s), 'g', 'dashed', 'rebalance steps')
 
  plt.title('km1 goal - minimize')
  plt.plot(num_nodes, km1lines, color='b', label='km1')
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

if __name__ == "__main__":
    main()
