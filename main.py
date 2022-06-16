from testing import *
from clustering import *
from gen_input import *


def main():
    algo_clustering_to_file(index_size=6)
    # print(test_stats())
    # test_times(True, True, index_size=8)


if __name__ == "__main__":
    main()