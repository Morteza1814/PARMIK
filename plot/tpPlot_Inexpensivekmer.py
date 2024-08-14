import matplotlib.pyplot as plt
import sys
import numpy as np
import matplotlib.ticker as ticker

# Initialize empty lists for x, y1, y2, y3, y4, and y5
x = [] # alignment sizes
y1 = [] # Baseline
y2 = [] # IKT0
y3 = [] # IKT0 + polishing
y4 = [] # IKT5000
y5 = [] # IKT5000 + polishing
y6 = [] # BLAST

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 9:
    print("Usage: python script.py input_file output_file first_ax_name second_ax_name third_ax_name fourth_ax_name 5th_ax_name 6th_ax_name")
    sys.exit(1)

# Extract input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
first_ax_name = sys.argv[3]
second_ax_name = sys.argv[4]
third_ax_name = sys.argv[5]
fourth_ax_name = sys.argv[6]
fifth_ax_name = sys.argv[7]
sixth_ax_name = sys.argv[8]

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
        y6.append(float(data[6]))  # Convert y6 value to float

# Create a figure and axes for both plots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

bar_width = 0.15
# Plot the first bar chart in the top subplot
# ax1.bar(x, y1, width=bar_width, label=first_ax_name, color='skyblue')
# ax1.bar([i + bar_width for i in x], y2, width=bar_width, label=second_ax_name, color='firebrick')  # Shift x values for Dataset 2
# ax1.bar([i + 2 * bar_width for i in x], y3, width=bar_width, label=third_ax_name, color='green')  # Shift x values for Dataset 3
# ax1.bar([i + 3 * bar_width for i in x], y4, width=bar_width, label=fourth_ax_name, color='purple')  # Shift x values for Dataset 4
# ax1.bar([i + 4 * bar_width for i in x], y5, width=bar_width, label=fifth_ax_name, color='orange')  # Shift x values for Dataset 5
# ax1.legend(fontsize=18)
# ax1.set_yscale('log')  # Set logarithmic scale on the y-axis
# ax1.tick_params(axis='both', which='major', labelsize=18)
# ax1.set_ylabel('No. of TP (log)', fontsize=20)
# ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

ax1.plot(x, y1, 'o-', label=first_ax_name, color='black', markersize=10)
ax1.plot(x, y2, 's-', label=second_ax_name, color='green', markersize=10)
ax1.plot(x, y4, '^-', label=fourth_ax_name, color='blue', markersize=10)
ax1.plot(x, y6, 'v-', label=sixth_ax_name, color='firebrick', markersize=10)

ax1.legend(fontsize=18)
ax1.set_yscale('log')  # Set logarithmic scale on the y-axis
ax1.tick_params(axis='both', which='major', labelsize=18)
ax1.set_ylabel('No. of TP (log)', fontsize=20)
ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')


# Specify the range for the zoomed-in plot
zoom_start = int(26 - x[0])
zoom_end = int(36 - x[0])

print("Zoom Start:", zoom_start, " Zoom End:", zoom_end)
# Plot the zoomed-in portion as a bar chart in the bottom subplot
ax2.bar(x[zoom_start:zoom_end], y1[zoom_start:zoom_end], label=first_ax_name, color='black', width=bar_width)
ax2.bar([i + bar_width for i in x[zoom_start:zoom_end]], y2[zoom_start:zoom_end], label=second_ax_name, color='green', width=bar_width)  # Shift x values for Dataset 2
ax2.bar([i + 2 * bar_width for i in x[zoom_start:zoom_end]], y3[zoom_start:zoom_end], label=third_ax_name, color='lightgreen', width=bar_width)  # Shift x values for Dataset 3
ax2.bar([i + 3 * bar_width for i in x[zoom_start:zoom_end]], y4[zoom_start:zoom_end], label=fourth_ax_name, color='blue', width=bar_width)  # Shift x values for Dataset 4
ax2.bar([i + 4 * bar_width for i in x[zoom_start:zoom_end]], y5[zoom_start:zoom_end], label=fifth_ax_name, color='skyblue', width=bar_width)  # Shift x values for Dataset 5
ax2.bar([i + 5 * bar_width for i in x[zoom_start:zoom_end]], y6[zoom_start:zoom_end], label=sixth_ax_name, color='firebrick', width=bar_width)  # Shift x values for Dataset 5
# ax2.legend(fontsize=18)
ax2.legend(fontsize=18, bbox_to_anchor=(.05, 1.1), loc='upper left', borderaxespad=0., ncol=2)
ax2.tick_params(axis='both', which='major', labelsize=18)

# Format the y-axis ticks
formatter = ticker.FuncFormatter(lambda x, pos: '0' if x == 0 else '{:.0f}'.format(x*1e-5))
ax2.yaxis.set_major_formatter(formatter)
scale_factor = 1e5 
ax2.set_ylabel(f'No. of TP (×10^{int(np.log10(scale_factor))})',  fontsize=20)

ax2.text(-0.1, 1.1, '(B)', transform=ax2.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

plt.xlabel('Alignment Size', fontsize=20)

# Show the plot
plt.tight_layout()
plt.show()

plt.savefig(output_file)
