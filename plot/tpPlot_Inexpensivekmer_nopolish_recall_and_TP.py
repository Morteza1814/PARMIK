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
if len(sys.argv) != 10:
    print("Usage: python script.py input_file output_file second_ax_name third_ax_name fourth_ax_name fifth_ax_name sixth_ax_name input_file2")
    sys.exit(1)

# Extract input and output file paths from command-line arguments
input_file1 = sys.argv[1]
output_file = sys.argv[2]
first_ax_name = sys.argv[3]
second_ax_name = sys.argv[4]
third_ax_name = sys.argv[5]
fourth_ax_name = sys.argv[6]
fifth_ax_name = sys.argv[7]
sixth_ax_name = sys.argv[8]
input_file2 = sys.argv[9]

# Read the data from the input file
with open(input_file1, 'r') as file:
    for line in file:
        data = line.strip().split('\t')
        x.append(float(data[0]))  # Convert x value to float
        y1.append(float(data[1]))  # Convert y1 value to float
        y2.append(float(data[2]))  # Convert y2 value to float
        y3.append(float(data[3]))  # Convert y3 value to float
        y4.append(float(data[4]))  # Convert y4 value to float
        y5.append(float(data[5]))  # Convert y5 value to float

# Initialize empty lists for x, y1, y2, y3, y4, and y5
tp_x = [] # alignment sizes
tp_y1 = [] # Baseline
tp_y2 = [] # IKT0
tp_y3 = [] # IKT0 + polishing
tp_y4 = [] # IKT5000
tp_y5 = [] # IKT5000 + polishing
tp_y6 = [] # BLAST

# Read the data from the input file
with open(input_file2, 'r') as file:
    for line in file:
        data = line.strip().split('\t')
        tp_x.append(float(data[0]))  # Convert x value to float
        tp_y1.append(float(data[1]))  # Convert y1 value to float
        tp_y2.append(float(data[2]))  # Convert y2 value to float
        tp_y3.append(float(data[3]))  # Convert y3 value to float
        tp_y4.append(float(data[4]))  # Convert y4 value to float
        tp_y5.append(float(data[5]))  # Convert y5 value to float
        tp_y6.append(float(data[6]))  # Convert y6 value to float

# Create a figure and axes for both plots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 15))

bar_width = 0.15

# Plot using dots
ax1.plot(x, y1, 's-', label=second_ax_name, color='navy', markersize=10)
ax1.plot(x, y2, '^-', label=third_ax_name, color='darkseagreen', markersize=10)
ax1.plot(x, y5, 'v-', label=sixth_ax_name, color='firebrick', markersize=10)

ax1.legend(fontsize=16, loc='lower right')
ax1.tick_params(axis='both', which='major', labelsize=18)
ax1.set_ylabel('Recall Rate (%)', fontsize=20)
ax1.set_xlabel('Minimum Alignment Size (R)', fontsize=20)
ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

# Specify the range for the zoomed-in plot
zoom_start = int(26 - x[0])
zoom_end = int(36 - x[0])

print("Zoom Start:", zoom_start, " Zoom End:", zoom_end)

# Plot the zoomed-in portion as a bar chart in the bottom subplot
ax2.bar(x[zoom_start:zoom_end], y1[zoom_start:zoom_end], label=second_ax_name, color='navy', width=bar_width)
ax2.bar([i + bar_width for i in x[zoom_start:zoom_end]], y3[zoom_start:zoom_end], label=fourth_ax_name, color='skyblue', width=bar_width)  # Shift x values for Dataset 3
ax2.bar([i + 2 * bar_width for i in x[zoom_start:zoom_end]], y2[zoom_start:zoom_end], label=third_ax_name, color='green', width=bar_width)  # Shift x values for Dataset 4
ax2.bar([i + 3 * bar_width for i in x[zoom_start:zoom_end]], y4[zoom_start:zoom_end], label=fifth_ax_name, color='darkseagreen', width=bar_width)  # Shift x values for Dataset 4
ax2.bar([i + 4 * bar_width for i in x[zoom_start:zoom_end]], y5[zoom_start:zoom_end], label=sixth_ax_name, color='firebrick', width=bar_width)  # Shift x values for Dataset 4

ax2.legend(fontsize=16, loc='lower right')
ax2.tick_params(axis='both', which='major', labelsize=18)

# Set y-axis limits to 50-100 for the zoomed-in plot
ax2.set_ylim(80, 100)

# Format the y-axis ticks
# formatter = ticker.FuncFormatter(lambda x, pos: '0' if x == 0 else '{:.0f}k'.format(x*1e-3))
# ax2.yaxis.set_major_formatter(formatter)
ax2.set_ylabel('Recall Rate (%)', fontsize=20)
ax2.set_xlabel('Minimum Alignment Size (R)', fontsize=20)

ax2.text(-0.1, 1.1, '(B)', transform=ax2.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

#ax3
ax3.bar(x[zoom_start:zoom_end], tp_y1[zoom_start:zoom_end], label=first_ax_name, color='grey', width=bar_width)
ax3.bar([i + bar_width for i in tp_x[zoom_start:zoom_end]], tp_y2[zoom_start:zoom_end], label=second_ax_name, color='navy', width=bar_width)  # Shift x values for Dataset 2
ax3.bar([i + 2 * bar_width for i in tp_x[zoom_start:zoom_end]], tp_y3[zoom_start:zoom_end], label=fourth_ax_name, color='skyblue', width=bar_width)  # Shift x values for Dataset 3
ax3.bar([i + 3 * bar_width for i in tp_x[zoom_start:zoom_end]], tp_y4[zoom_start:zoom_end], label=third_ax_name, color='green', width=bar_width)  # Shift x values for Dataset 4
ax3.bar([i + 4 * bar_width for i in tp_x[zoom_start:zoom_end]], tp_y5[zoom_start:zoom_end], label=fifth_ax_name, color='darkseagreen', width=bar_width)  # Shift x values for Dataset 5
ax3.bar([i + 5 * bar_width for i in tp_x[zoom_start:zoom_end]], tp_y6[zoom_start:zoom_end], label=sixth_ax_name, color='firebrick', width=bar_width)  # Shift x values for Dataset 5
# ax3.legend(fontsize=18)
# ax3.legend(fontsize=18, bbox_to_anchor=(.05, 1.1), loc='upper left', borderaxespad=0., ncol=2)
ax3.legend(fontsize=16, loc='upper right', ncols=2)
ax3.tick_params(axis='both', which='major', labelsize=18)

# Format the y-axis ticks
formatter = ticker.FuncFormatter(lambda x, pos: '0' if x == 0 else '{:.0f}'.format(x*1e-5))
ax3.yaxis.set_major_formatter(formatter)
scale_factor = 1e5 
ax3.set_ylabel(f'No. of TP (×10^{int(np.log10(scale_factor))})',  fontsize=20)

ax3.text(-0.1, 1.1, '(C)', transform=ax3.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')
ax3.set_xlabel('Alignment Size', fontsize=20)


# Show the plot
plt.tight_layout()
plt.show()

plt.savefig(output_file)
