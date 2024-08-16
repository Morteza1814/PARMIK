import matplotlib.pyplot as plt
import sys
import numpy as np
import matplotlib.ticker as ticker

# Initialize empty lists for x, y1, y3, and y4
x = []
y1 = [] # t0+polish
y2 = [] # t0
y3 = [] # t5000+polish
y4 = [] # t5000
y5 = [] # blast

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 8:
    print("Usage: python script.py input_file output_file second_ax_name third_ax_name fourth_ax_name fifth_ax_name sixth_ax_name")
    sys.exit(1)

# Extract input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
second_ax_name = sys.argv[3]
third_ax_name = sys.argv[4]
fourth_ax_name = sys.argv[5]
fifth_ax_name = sys.argv[6]
sixth_ax_name = sys.argv[7]

# Read the data from the input file
with open(input_file, 'r') as file:
    for line in file:
        data = line.strip().split('\t')
        x.append(float(data[0]))  # Convert x value to float
        y1.append(float(data[1]))  # Convert y1 value to float
        y2.append(float(data[2]))  # Convert y2 value to float
        y3.append(float(data[3]))  # Convert y3 value to float
        y4.append(float(data[4]))  # Convert y4 value to float
        y5.append(float(data[5]))  # Convert y5 value to float

# Create a figure and axes for both plots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

bar_width = 0.15

# Plot using dots
ax1.plot(x, y1, 's-', label=second_ax_name, color='navy', markersize=10)
ax1.plot(x, y2, '^-', label=third_ax_name, color='darkseagreen', markersize=10)
ax1.plot(x, y5, 'v-', label=sixth_ax_name, color='firebrick', markersize=10)

ax1.legend(fontsize=18, loc='lower right')
ax1.tick_params(axis='both', which='major', labelsize=18)
ax1.set_ylabel('Recall Rate (%)', fontsize=20)
ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

# Specify the range for the zoomed-in plot
zoom_start = int(26 - x[0])
zoom_end = int(50 - x[0])

print("Zoom Start:", zoom_start, " Zoom End:", zoom_end)

# Plot the zoomed-in portion as a bar chart in the bottom subplot
ax2.bar(x[zoom_start:zoom_end], y1[zoom_start:zoom_end], label=second_ax_name, color='navy', width=bar_width)
ax2.bar([i + bar_width for i in x[zoom_start:zoom_end]], y3[zoom_start:zoom_end], label=fourth_ax_name, color='skyblue', width=bar_width)  # Shift x values for Dataset 3
ax2.bar([i + 2 * bar_width for i in x[zoom_start:zoom_end]], y2[zoom_start:zoom_end], label=third_ax_name, color='green', width=bar_width)  # Shift x values for Dataset 4
ax2.bar([i + 3 * bar_width for i in x[zoom_start:zoom_end]], y4[zoom_start:zoom_end], label=fifth_ax_name, color='darkseagreen', width=bar_width)  # Shift x values for Dataset 4
ax2.bar([i + 4 * bar_width for i in x[zoom_start:zoom_end]], y5[zoom_start:zoom_end], label=sixth_ax_name, color='firebrick', width=bar_width)  # Shift x values for Dataset 4

ax2.legend(fontsize=18, loc='lower right')
ax2.tick_params(axis='both', which='major', labelsize=18)

# Set y-axis limits to 50-100 for the zoomed-in plot
ax2.set_ylim(80, 100)

# Format the y-axis ticks
# formatter = ticker.FuncFormatter(lambda x, pos: '0' if x == 0 else '{:.0f}k'.format(x*1e-3))
# ax2.yaxis.set_major_formatter(formatter)
ax2.set_ylabel('Recall Rate (%)', fontsize=20)
ax2.text(-0.1, 1.1, '(B)', transform=ax2.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

# plt.xlabel('Alignment Size', fontsize=20)
plt.xlabel('Minimum Alignment Size (R)', fontsize=20)

# Show the plot
plt.tight_layout()
plt.show()

plt.savefig(output_file)
