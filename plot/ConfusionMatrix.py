import matplotlib.pyplot as plt

# Confusion matrix values
# RmYn01
# PM
false_positive_tool1 = 0
false_negative_tool1 = 0
true_positive_tool1 = 1434191
true_negative_tool1 = 33
# BWA
false_positive_tool2 = 3885
false_negative_tool2 = 0
true_positive_tool2 = 17834
true_negative_tool2 = 0

# RmYn02
# # PM
# false_positive_tool1 = 0
# false_negative_tool1 = 0
# true_positive_tool1 = 9332694
# true_negative_tool1 = 0
# # BWA
# false_positive_tool2 = 1354
# false_negative_tool2 = 0
# true_positive_tool2 = 28168
# true_negative_tool2 = 0

# NC_045512.2
# PM
# false_positive_tool1 = 0
# false_negative_tool1 = 0
# true_positive_tool1 = 5593063
# true_negative_tool1 = 5789
# # BWA
# false_positive_tool2 = 19279
# false_negative_tool2 = 100
# true_positive_tool2 = 7716
# true_negative_tool2 = 2678

# Create a figure and axis
fig, ax = plt.subplots()

# Draw rectangles with edges
rect1 = plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black')
rect2 = plt.Rectangle((1, 0), 1, 1, facecolor='white', edgecolor='black')
rect3 = plt.Rectangle((0, 1), 1, 1, facecolor='white', edgecolor='black')
rect4 = plt.Rectangle((1, 1), 1, 1, facecolor='white', edgecolor='black')

# Add rectangles to the plot
ax.add_patch(rect1)
ax.add_patch(rect2)
ax.add_patch(rect3)
ax.add_patch(rect4)

# Add text to rectangles with specified colors
ax.text(0.5, 0.65, f'{"{:,}".format(false_positive_tool1)}', ha='center', va='center', color='blue', fontsize=18)
# ax.text(0.5, 0.5, f'/', ha='center', va='center', color='black', fontsize=20)
ax.text(0.5, 0.35, f'{"{:,}".format(false_positive_tool2)}', ha='center', va='center', color='red', fontsize=18)

ax.text(1.5, 0.65, f'{"{:,}".format(false_negative_tool1)}', ha='center', va='center', color='blue', fontsize=18)
# ax.text(1.5, 0.5, f'/', ha='center', va='center', color='black', fontsize=20)
ax.text(1.5, 0.35, f'{"{:,}".format(false_negative_tool2)}', ha='center', va='center', color='red', fontsize=18)

ax.text(0.5, 1.65, f'{"{:,}".format(true_positive_tool1)}', ha='center', va='center', color='blue', fontsize=18)
# ax.text(0.5, 1.5, f'/', ha='center', va='center', color='black', fontsize=20)
ax.text(0.5, 1.35, f'{"{:,}".format(true_positive_tool2)}', ha='center', va='center', color='red', fontsize=18)

ax.text(1.5, 1.65, f'{"{:,}".format(true_negative_tool1)}', ha='center', va='center', color='blue', fontsize=18)
# ax.text(1.5, 1.5, f'/', ha='center', va='center', color='black', fontsize=20)
ax.text(1.5, 1.35, f'{"{:,}".format(true_negative_tool2)}', ha='center', va='center', color='red', fontsize=18)

# Set axis properties
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_xticks([0.5, 1.5])
ax.set_xticklabels(['Positive', 'Negative'], fontsize=25)
ax.set_yticks([0.5, 1.5])
ax.set_yticklabels(['False', 'True'], fontsize=25)

# Hide grid lines
ax.grid(False)

# Hide axes ticks
ax.tick_params(axis='both', which='both', length=0)

# Set the aspect ratio to be equal
ax.set_aspect('equal')

# Remove axis labels
ax.set_xlabel('')
ax.set_ylabel('')

# Set the title
# ax.set_title('PM (blue) vs. BWA (red)', fontsize=25)

plt.show()

plt.savefig('confmat_RmYn01.png')
# plt.savefig('confmat_RmYn02.png')
# plt.savefig('confmat_NC_045512.2.png')

