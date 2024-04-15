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

# Initialize a dictionary to store the data for each node
node_data = {49: []}

# Open the text file and read the data
with open('/Users/weizhiwang/Project-Chrono-Linker/ThirdParty/chrono_projects/hexa_test/build/tip_node.txt', 'r') as file:
    for line in file:
        # Split the line into words
        words = line.split()
        # Extract the time, node number, and z position
        time = float(words[1])
        node = int(words[4].rstrip(':'))
        y = float(words[-2].rstrip(',')) # Change this to -2 for vertical y position
        # Store the data in the dictionary
        node_data[node].append((time, y))

# Print the first few data points for each node to check if they are correct
#for node, data in node_data.items():
 #   print(f"First few data points for node {node}: {data[:5]}")

# Create a separate plot for each node
for node, data in node_data.items():
    # Unzip the list of tuples into two lists
    times, y_positions = zip(*data)
    # Create a new figure
    plt.figure()
    # Plot the z position as a function of time
    plt.plot(times, y_positions)
    # Add labels and a title
    plt.xlabel('Time')
    plt.ylabel('Y Position')
    plt.title(f'Node {node}')
    plt.xlim(0, max(times))
    #plt.ylim(-1.5, 1.5)
    # Show the plot
    plt.show()