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
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(18, 25))
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.95, hspace=0.4)

bar_width = 0.14

# Plot using dots
ax1.plot(x, y1, 's-', label=second_ax_name, color='navy', markersize=10)
ax1.plot(x, y2, '^-', label=third_ax_name, color='darkseagreen', markersize=10)
ax1.plot(x, y5, 'v-', label=sixth_ax_name, color='firebrick', markersize=10)

ax1.legend(fontsize=20, loc='lower right')
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_ylabel('Recall Rate (%)', fontsize=24)
ax1.set_xlabel('Minimum Alignment Size (R)', fontsize=24)
ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=26, fontweight='bold', va='top', ha='right', color='blue')

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

ax2.legend(fontsize=20, loc='lower right')
ax2.tick_params(axis='both', which='major', labelsize=20)

ax2.set_ylim(80, 100)
ax2.set_yticks([80, 85, 90, 95, 100])

# Format the y-axis ticks
# formatter = ticker.FuncFormatter(lambda x, pos: '0' if x == 0 else '{:.0f}k'.format(x*1e-3))
# ax2.yaxis.set_major_formatter(formatter)
ax2.set_ylabel('Recall Rate (%)', fontsize=24)
ax2.set_xlabel('Minimum Alignment Size (R)', fontsize=24)

ax2.text(-0.1, 1.1, '(B)', transform=ax2.transAxes, fontsize=26, fontweight='bold', va='top', ha='right', color='blue')

#ax3
bar_width = 0.12
ax3.bar(tp_x[zoom_start:zoom_end], tp_y1[zoom_start:zoom_end], label=first_ax_name, color='grey', width=bar_width)
ax3.bar([i + bar_width for i in tp_x[zoom_start:zoom_end]], tp_y2[zoom_start:zoom_end], label=second_ax_name, color='navy', width=bar_width)  # Shift x values for Dataset 2
ax3.bar([i + 2 * bar_width for i in tp_x[zoom_start:zoom_end]], tp_y3[zoom_start:zoom_end], label=fourth_ax_name, color='skyblue', width=bar_width)  # Shift x values for Dataset 3
ax3.bar([i + 3 * bar_width for i in tp_x[zoom_start:zoom_end]], tp_y4[zoom_start:zoom_end], label=third_ax_name, color='green', width=bar_width)  # Shift x values for Dataset 4
ax3.bar([i + 4 * bar_width for i in tp_x[zoom_start:zoom_end]], tp_y5[zoom_start:zoom_end], label=fifth_ax_name, color='darkseagreen', width=bar_width)  # Shift x values for Dataset 5
ax3.bar([i + 5 * bar_width for i in tp_x[zoom_start:zoom_end]], tp_y6[zoom_start:zoom_end], label=sixth_ax_name, color='firebrick', width=bar_width)  # Shift x values for Dataset 5

ax3.legend(fontsize=20, loc='lower right')
ax3.tick_params(axis='both', which='major', labelsize=20)

# Format the y-axis ticks
formatter = ticker.FuncFormatter(lambda x, pos: '0' if x == 0 else '{:.1f}'.format(x*1e-7))

ax3.yaxis.set_major_formatter(formatter)
scale_factor = 1e7 
ax3.set_ylabel(f'No. of TP (×10^{int(np.log10(scale_factor))})',  fontsize=24)

ax3.text(-0.1, 1.1, '(C)', transform=ax3.transAxes, fontsize=26, fontweight='bold', va='top', ha='right', color='blue')
ax3.set_xlabel('Minimum Alignment Size (R)', fontsize=24)

# ax4
# Data
x_values = np.array([132, 131, 130, 129, 128, 127, 126, 125, 124, 123, 122, 121, 120, 119, 
                     118, 117, 116, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106, 105, 
                     104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 
                     88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 
                     71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 
                     54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 
                     37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 
                     20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 
                     1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, 
                     -15, -16, -17, -18])

y_values = np.array([1, 2, 0, 0, 0, 0, 0, 0, 0, 2, 7, 3, 2, 2, 0, 0, 0, 1, 1, 1, 5, 0, 1, 
                     1, 0, 0, 0, 0, 10, 14, 1, 1, 1, 98, 0, 0, 0, 0, 1, 1, 39, 2, 60, 23, 
                     0, 0, 3, 5, 28, 4, 8, 14, 6, 25, 2, 1, 1, 1, 1, 1, 7, 8, 12, 14, 32, 
                     4, 3, 2, 114, 9, 16, 17, 9, 68, 28, 8, 7, 161, 14, 26, 69, 32, 29, 
                     154, 60, 16, 23, 63, 12, 351, 73, 56, 137, 307, 167, 406, 68, 41, 157, 
                     229, 184, 286, 658, 369, 365, 359, 579, 1393, 20083, 1729, 1559, 3121, 
                     2596, 2603, 5222, 4131, 3534, 18389, 9950, 10519, 22773, 16340, 20420, 
                     47410, 52987, 75389, 757787, 519113, 531663, 2347512, 5060046, 714383, 
                     36131627, 7754, 9050, 2842, 6259, 449, 2066, 875, 144, 8, 2, 1, 2, 0, 
                     0, 2, 2, 2, 104])

# Set colors based on x-values
colors = ['firebrick' if x < 0 else 'navy' if x > 0 else 'grey' for x in x_values]

# Labels
labels = ['PARMIK IKT0 < Baseline' if x < 0 else 'PARMIK IKT0 > Baseline' if x > 0 else 'PARMIK IKT0 = Baseline' for x in x_values]
ax4.bar(x_values, y_values, color=colors, label=labels)
ax4.set_yscale('log') 

# Adding custom legend
legend_labels = ['PARMIK IKT0+Polishing < Baseline', 'PARMIK IKT0+Polishing = Baseline', 'PARMIK IKT0+Polishing > Baseline']
legend_colors = ['firebrick', 'grey', 'navy']
handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in legend_colors]
ax4.legend(handles, legend_labels, fontsize=20)
ax4.tick_params(axis='both', which='major', labelsize=20)
ax4.set_xlabel('Number of Base Pairs', fontsize=24)
ax4.set_ylabel('Number of Alignments', fontsize=24)
ax4.text(-0.1, 1.1, '(D)', transform=ax4.transAxes, fontsize=26, fontweight='bold', va='top', ha='right', color='blue')


# Show the plot
plt.tight_layout()
plt.show()

plt.savefig(output_file)
