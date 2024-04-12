import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Initialize an empty dictionary to store the data for each node
node_data = {}

# Open the text file and read the data
with open('/Users/weizhiwang/Project-Chrono-Linker/ThirdParty/chrono_projects/hexa_test/build/node_positions.txt', 'r') as file:
    for line in file:
        # Split the line into words
        words = line.split()
        # Check if the first word is a number
        if words[0].rstrip(':').isdigit():
            # Extract the node number and coordinates
            node = int(words[0].rstrip(':'))
            x = float(words[1].rstrip(','))
            y = float(words[2].rstrip(','))
            z = float(words[3])
            # Store the data in the dictionary
            node_data[node] = (x, y, z)

# Create a new figure
fig = plt.figure()

# Create a 3D subplot
ax = fig.add_subplot(111, projection='3d')

# Plot each node
for node, coordinates in node_data.items():
    ax.scatter(*coordinates, label=f'Node {node}')

# Add labels and a title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Node Positions')

# Add a legend
ax.legend()

# Show the plot
plt.show()