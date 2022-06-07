import matplotlib.pyplot as plt
import os

def main():
  path = "/home/timo/Coding/mykahypar/partitioning_results/data"
  km1path = path + "km1.txt"
  imbpath = path + "imbalances.txt"
  showParallel(km1path, imbpath, "Km1 - minimize", "Imbalance - less than epsilon")

def showParallel(path, path2, title, title2) :
  newpid = os.fork()
  if (newpid == 0):
    openAndShow(path2, title2)
  else:
    openAndShow(path, title)

 

def openAndShow(path, title):
 with open(path) as file:
    lines = file.read().splitlines()
    _, ax = plt.subplots()
    plt.plot(lines)
    plt.title(title)
    for n, label in enumerate(ax.yaxis.get_ticklabels()):
      if n % 2 != 0 or n % 3 != 0:
          label.set_visible(False)
    plt.show()

if __name__ == "__main__":
    main()
