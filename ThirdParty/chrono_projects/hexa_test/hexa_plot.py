import matplotlib.pyplot as plt

# Read the text file
with open('/Users/weizhiwang/Project-Chrono-Linker/ThirdParty/chrono_projects/hexa_test/build/node_positions.txt', 'r') as file:
    lines = file.readlines()

# Extract the coordinates
coordinates = []
nodes = []
for line in lines:
    if line.startswith('Node'):
        parts = line.split(':')
        node = int(parts[0].split()[1])
        coords = [float(x) for x in parts[1].split(',')]
        coordinates.append(coords)
        nodes.append(node)

# Separate the x, y, and z coordinates
x = [coord[0] for coord in coordinates]
y = [coord[1] for coord in coordinates]
z = [coord[2] for coord in coordinates]

# Plot the coordinates
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, z, y)
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_zlabel('Y')

# Adjust the aspect ratio of the axes
#ax.set_box_aspect([1, 1, 0.05])  # Modify the second value to adjust the scale of the y-axis

# Add text annotations for each node
for i in range(len(nodes)):
    ax.text(x[i], z[i], y[i], f'Node {nodes[i]}')

plt.show()