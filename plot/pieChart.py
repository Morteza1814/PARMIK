import matplotlib.pyplot as plt
import argparse

def read_data_from_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            label, percentages = parts[0], list(map(float, parts[1:]))
            data.append((label, percentages))
    print(data)
    return data

def create_pie_chart(ax, data, title, colors):
    # sizes = [value for _, value in data]
    ax.pie(data, colors=colors, autopct=lambda p: '{:.1f}%'.format(p), startangle=140, textprops={'fontsize': 18})
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax.set_title(title, fontsize=18, loc='center')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate two pie charts from a data file.")
    parser.add_argument("--input", help="Path to the data file.")
    parser.add_argument("--title1", help="Title of the first pie chart.")
    parser.add_argument("--title2", help="Title of the second pie chart.")
    parser.add_argument("--output", help="Output filename (without extension) for saving the charts.")
    args = parser.parse_args()

    data = read_data_from_file(args.input)

    # Extracting data from the array
    first_label, first_values = data[0]
    second_label, second_values = data[1]
    third_label, third_values = data[2]

    # Splitting data for two pie charts
    # sorted_labels = sorted([first_label, second_label, third_label])

    # Creating the sets with sorted labels
    # first_data = {(sorted_labels[0], first_values[0]), (sorted_labels[1], second_values[0]), (sorted_labels[2], third_values[0])}
    # second_data = {(sorted_labels[0], first_values[1]), (sorted_labels[1], second_values[1]), (sorted_labels[2], third_values[1])}

    first_data = [first_values[0], second_values[0], third_values[0]]
    second_data = [first_values[1], second_values[1], third_values[1]]

    print("1st data: ", first_data)
    print("2nd data:", second_data)

    labels = [first_label, second_label, third_label]
    colors = ['yellowgreen', 'skyblue', 'orange']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    create_pie_chart(ax1, first_data, args.title1, colors)
    create_pie_chart(ax2, second_data, args.title2, colors)
    
    plt.legend(labels, fontsize=11, loc='best', shadow=True, prop={'weight':'bold'})

    plt.tight_layout()  # Adjust layout to prevent overlap of titles and labels
    plt.savefig(args.output + ".png")
    plt.show()
