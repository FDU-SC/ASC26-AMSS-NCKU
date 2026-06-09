import numpy as np
import argparse

def judge(std_file, result_file):
    std = np.loadtxt(
        std_file,
        comments='#',
        delimiter=None,
        dtype=float
    )
    result = np.loadtxt(
        result_file,
        comments='#',
        delimiter=None,
        dtype=float
    )

    # print(std.shape, result.shape)

    count = 0
    sum = 0
    for i in range(min(std.shape[0],result.shape[0])):
        for j in [1,2,4,5]:
            maxn = max(abs(std[i][j]) , abs(result[i][j]))
            sum += ((std[i][j] - result[i][j]) / maxn) ** 2 if maxn > 0 else 0
            count += 1
    error = np.sqrt(sum/count)*100
    if error < 1:
        print("\033[32mAccepted\033[0m: {}%".format(error))
    else:
        print("\033[31mRejected\033[0m: {}%".format(error))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two files using the judge function.")
    parser.add_argument("--std_file", default="bssn_BH.dat")
    parser.add_argument("--result_file", default="binary_output/bssn_BH.dat")
    args = parser.parse_args()
    judge(args.std_file, args.result_file)