import matplotlib.pyplot as plt

# Node locations
node_locations = [
    [0, 0, 0],
    [1, 0, 0],
    [2, 0, 0],
    [3, 0, 0],
    [4, 0, 0],
    [5, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [2, 1, 0],
    [3, 1, 0],
    [4, 1, 0],
    [5, 1, 0],
    [0, 2, 0],
    [1, 2, 0],
    [2, 2, 0],
    [3, 2, 0],
    [4, 2, 0],
    [5, 2, 0],
    [0, 3, 0],
    [1, 3, 0],
    [2, 3, 0],
    [3, 3, 0],
    [4, 3, 0],
    [5, 3, 0],
    [0, 4, 0],
    [1, 4, 0],
    [2, 4, 0],
    [3, 4, 0],
    [4, 4, 0],
    [5, 4, 0],
    [0, 5, 0],
    [1, 5, 0],
    [2, 5, 0],
    [3, 5, 0],
    [4, 5, 0],
    [5, 5, 0]
]

# Plotting node locations
fig, ax = plt.subplots()
for i, loc in enumerate(node_locations):
    ax.scatter(loc[0], loc[1], label=f"Node {i}", color='blue')
    ax.text(loc[0], loc[1], f"Node {i}", ha='center', va='bottom')

ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_title('Node Locations')
#ax.legend()
plt.show()