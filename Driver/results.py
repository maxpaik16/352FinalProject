
# file to plot stability results

# necessary imports
import csv
import matplotlib.pyplot as plt

# read in results from csv
stabilities = []
with open('training.csv') as datafile:
    csvReader = csv.reader(datafile)
    for row in csvReader:
        stabilities.append(int(row[-1]))

# find percentage of orbits that were unstable
percentage_unstable = sum(stabilities) / len(stabilities)
print(percentage_unstable)

# variables for bar chart of instabilities of our simulations
unstables = sum(stabilities)
stables = len(stabilities) - unstables
options = ("Stable Simulations", "Unstable Simulations")
xaxis = [0, 1]

# create plot
plt.bar(xaxis, [stables, unstables], align='center')
plt.xticks(xaxis, options)
plt.ylabel("Count")

plt.title("Stabilitiy of Simulated Planetary Motions")
plt.savefig("stabilities.png")
