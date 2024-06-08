import sys

def process_file(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
        return

    # Dictionary to store counts
    string_counts = {}
    low_complexity_counts = {
        'all_same_char': 0,
        '0.85lowComplex': 0,
        '0.6lowComplex': 0,
        '0.4lowComplex': 0,
        '0.5lowComplex': 0
    }

    # Process each line
    for line in lines:
        first_col, value = line.split()
        value = int(value)
        
        # Check for low complexity strings (all same character)
        if len(set(first_col)) == 1:
            low_complexity_counts['all_same_char'] += 1
        # Check for mostly Gs
        if first_col.count('G') / len(first_col) > 0.85 or \
            first_col.count('C') / len(first_col) > 0.85 or \
            first_col.count('A') / len(first_col) > 0.85 or \
            first_col.count('T') / len(first_col) > 0.85:
            low_complexity_counts['0.85lowComplex'] += 1
        if first_col.count('G') / len(first_col) > 0.6 or \
            first_col.count('C') / len(first_col) > 0.6 or \
            first_col.count('A') / len(first_col) > 0.6 or \
            first_col.count('T') / len(first_col) > 0.6:
            low_complexity_counts['0.6lowComplex'] += 1
        if first_col.count('G') / len(first_col) > 0.5 or \
            first_col.count('C') / len(first_col) > 0.5 or \
            first_col.count('A') / len(first_col) > 0.5 or \
            first_col.count('T') / len(first_col) > 0.5:
            low_complexity_counts['0.5lowComplex'] += 1
        if first_col.count('G') / len(first_col) > 0.4 or \
            first_col.count('C') / len(first_col) > 0.4 or \
            first_col.count('A') / len(first_col) > 0.4 or \
            first_col.count('T') / len(first_col) > 0.4:
            low_complexity_counts['0.4lowComplex'] += 1

    print("Low complexity counts:", low_complexity_counts)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
    else:
        filename = sys.argv[1]
        process_file(filename)
