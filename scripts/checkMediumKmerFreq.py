import re
import sys

# Function to check if there is any number over 5000 in the line
def has_number_over_m_below_n(line, m, n):
    # Extract all numbers within parentheses
    numbers = re.findall(r'\((\d+)\)', line)
    # Check if any of the numbers is greater than 5000
    for number in numbers:
        if int(number) > m and int(number) < n:
            return True
    return False

# Counter for lines with numbers over 5000
lines_with_numbers_in_range = 0

input_file = sys.argv[1]
m = int(sys.argv[2])
n = int(sys.argv[3])
# Read the file
with open(input_file, 'r') as file:
    for line in file:
        if has_number_over_m_below_n(line, m, n):
            lines_with_numbers_in_range += 1

print(f'Number of lines with numbers in range [{m}, {n}]: {lines_with_numbers_in_range}')