import matplotlib.pyplot as plt
import sys
import numpy as np
import matplotlib.ticker as ticker

# Initialize empty lists for x, y2, and y3
x = []
y2 = []
y3 = []

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 5:
    print("Usage: python script.py input_file output_file second_ax_name third_ax_name")
    sys.exit(1)

# Extract input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
second_ax_name = sys.argv[3]
third_ax_name = sys.argv[4]

with open(input_file, 'r') as file:
    for line in file:
        data = line.strip().split('\t')
        x.append(float(data[0]))  # Convert x value to float
        y2.append(float(data[1]))  # Convert y2 value to float
        y3.append(float(data[2]))  # Convert y3 value to float

# Create a figure and axes for both plots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

bar_width = 0.2

# Plot using dots
ax1.plot(x, y2, '-', label=second_ax_name, color='navy', markersize=10)
ax1.plot(x, y3, '-', label=third_ax_name, color='firebrick', markersize=10)

ax1.legend(fontsize=18, loc='lower right')
ax1.tick_params(axis='both', which='major', labelsize=18)
ax1.set_ylabel('Recall Rate (%)', fontsize=20)
ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

# Specify the range for the zoomed-in plot
zoom_start = int(21 - x[0])
zoom_end = int(35 - x[0])

print("Zoom Start:", zoom_start, " Zoom End:", zoom_end)

# Plot the zoomed-in portion as a bar chart in the bottom subplot
ax2.bar(x[zoom_start:zoom_end], y2[zoom_start:zoom_end], label=second_ax_name, color='navy', width=bar_width)
ax2.bar([i + bar_width for i in x[zoom_start:zoom_end]], y3[zoom_start:zoom_end], label=third_ax_name, color='firebrick', width=bar_width)  # Shift x values for Dataset 3
ax2.legend(fontsize=18, loc='lower right')
ax2.tick_params(axis='both', which='major', labelsize=18)

# Set y-axis limits to 50-100 for the zoomed-in plot
ax2.set_ylim(20, 100)

# Format the y-axis ticks
# formatter = ticker.FuncFormatter(lambda x, pos: '0' if x == 0 else '{:.0f}k'.format(x*1e-3))
# ax2.yaxis.set_major_formatter(formatter)
ax2.set_ylabel('Recall Rate (%)', fontsize=20)
ax2.text(-0.1, 1.1, '(B)', transform=ax2.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

plt.xlabel('Minimum Alignment Size (R)', fontsize=20)

# Show the plot
plt.tight_layout()
plt.show()

plt.savefig(output_file)
