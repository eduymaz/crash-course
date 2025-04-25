import sys

def main(input_file, output_file):
    with open(input_file) as inf, open(output_file, "w") as outf:
        data = inf.read().upper()
        outf.write(data)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2]) 