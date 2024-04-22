import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def read_data_from_file(filename):
    x = []
    y1 = []
    y2 = []
    with open(filename, 'r') as file:
        for line in file:
            data = line.strip().split(',')
            x.append(float(data[0]))
            y1.append(float(data[1]))
            y2.append(float(data[2]))
    return x, y1, y2

def main():
    parser = argparse.ArgumentParser(description='Plot data with two y-axes.')
    parser.add_argument('filename', help='Name of the file containing data')
    parser.add_argument('--xlabel', default='X Axis', help='Label for the x-axis')
    parser.add_argument('--ylabel1', default='Y1 Axis', help='Label for the first y-axis')
    parser.add_argument('--ylabel2', default='Y2 Axis', help='Label for the second y-axis')
    parser.add_argument('--output', default='output.png', help='Name of the output image file')
    args = parser.parse_args()

    try:
        x, y1, y2 = read_data_from_file(args.filename)
    except FileNotFoundError:
        print("File not found.")
        return
    except IndexError:
        print("File format is incorrect.")
        return

    if len(y1) != len(y2):
        print("The lengths of the two arrays are not the same.")
        return

    fig, ax1 = plt.subplots()

    line1, = ax1.plot(x, y1, 'b-', label=args.ylabel1)
    ax1.set_xlabel(args.xlabel, fontsize=14)
    # ax1.set_ylabel(args.ylabel1, color='b', fontsize=14)
    ax1.spines['left'].set_color('blue')  # Set color for y-axis line
    ax1.tick_params(axis='y', colors='blue')  # Set color for y-axis ticks

    ax2 = ax1.twinx()
    line2, = ax2.plot(x, y2, 'r-', label=args.ylabel2)
    # ax2.set_ylabel(args.ylabel2, color='r', fontsize=14)
    ax2.spines['right'].set_color('red')  # Set color for y-axis line
    ax2.tick_params(axis='y', colors='red')  # Set color for y-axis ticks

    # Adjust y-axis scales
    max_y1 = max(y1)
    max_y2 = max(y2)
    ax1.set_ylim(0, max_y1)
    ax2.set_ylim(0, max_y2)

    # Use FuncFormatter to shorten numbers for y-axis 1
    formatter1 = ticker.FuncFormatter(lambda x, pos: f'{x / 1e3:.0f}k' if x >= 1e3 else f'{x:.0f}')
    ax1.yaxis.set_major_formatter(formatter1)

    # Use FuncFormatter to shorten numbers for y-axis 2
    formatter2 = ticker.FuncFormatter(lambda x, pos: f'{x / 1e6:.0f}M' if x >= 1e6 else f'{x:.0f}')
    ax2.yaxis.set_major_formatter(formatter2)

    # Increase font size of tick labels
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)

    ax1.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax2.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    
    # Add legend
    lines = [line1, line2]
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, loc='upper center', fontsize=14)
    plt.show()

    plt.savefig(args.output)
    print(f"Figure saved as {args.output}")

if __name__ == "__main__":
    main()
