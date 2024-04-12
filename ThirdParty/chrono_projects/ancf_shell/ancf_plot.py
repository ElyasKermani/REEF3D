import matplotlib.pyplot as plt

# Initialize a dictionary to store the data for each node
node_data = {9: [], 14: [], 21: [], 25: []}

# Open the text file and read the data
with open('/Users/weizhiwang/Project-Chrono-Linker/ThirdParty/chrono_projects/hexa_test/build/node_positions.txt', 'r') as file:
    for line in file:
        # Split the line into words
        words = line.split()
        # Extract the time, node number, and z position
        time = float(words[1])
        node = int(words[4].rstrip(':'))
        z = float(words[-1])
        # Store the data in the dictionary
        node_data[node].append((time, z))

# Print the first few data points for each node to check if they are correct
#for node, data in node_data.items():
 #   print(f"First few data points for node {node}: {data[:5]}")

# Create a separate plot for each node
for node, data in node_data.items():
    # Unzip the list of tuples into two lists
    times, z_positions = zip(*data)
    # Create a new figure
    plt.figure()
    # Plot the z position as a function of time
    plt.plot(times, z_positions)
    # Add labels and a title
    plt.xlabel('Time')
    plt.ylabel('Z Position')
    plt.title(f'Node {node}')
    plt.xlim(0, max(times))
    plt.ylim(-1.5, 1.5)
    # Show the plot
    plt.show()