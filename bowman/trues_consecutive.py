import itertools
from itertools import chain, dropwhile, takewhile

test_cases = (([1, 2, 3, 10, 11, 12, 4, 5, 6], [4, 5, 6, 1, 2, 3]),
              ([], []),
              ([1], [1]),
              ([1, 2, 3, 4], [1, 2, 3, 4]),
              ([10], []),
              ([10, 11, 12, 13], []),
              ([1, 10], [1]),
              ([1, 2, 3, 10, 11, 12], [1, 2, 3]),
              ([10, 11, 1, 2], [1, 2]),
              ([13, 14, 15, 1, 2, 3, 10, 11, 12], [1, 2, 3]),
              ([1, 2, 3, 10, 4, 5], [4, 5, 1, 2, 3]),
              ([10, 11, 1, 13, 14], [1]),
              ([1, 10, 2], [2, 1]),
              ([1, 2, 10, 11, 3], [3, 1, 2]),
              ([1, 10, 11, 2, 3], [2, 3, 1]))


def trues_consecutive(lst, pred):
    # TODO: Clean up?
    T1 = takewhile(pred, lst)
    T2 = filter(pred, dropwhile(pred, lst))
    return chain(T2, T1)

if __name__ == "__main__":
    for test_case, answer in test_cases:
        output = list(trues_consecutive(test_case, lambda x: x < 10))
        if output != answer:
            print(test_case, output)
    print("all test cases pass!")
