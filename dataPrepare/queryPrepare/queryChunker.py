import sys
import os

def split_file(input_file, output_dir, batch_size=1000):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    num_queries = int(len(lines) / 2)
    num_batches = (num_queries) // batch_size + (1 if num_queries % batch_size != 0 else 0)
    line_num = 0
    for i in range(num_batches):
        start = i * batch_size
        end = min((i + 1) * batch_size, len(lines))
        batch_query = 0
        output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_{start}-{end-1}")
        with open(output_file, 'w') as of:
            batch_query = 0
            while batch_query < batch_size:
                if line_num >= len(lines):
                    break
                line = lines[line_num]
                line_num += 1
                if "SARS" not in line:
                    batch_query+=1
                of.writelines(line)
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_directory>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    split_file(input_file, output_dir)
