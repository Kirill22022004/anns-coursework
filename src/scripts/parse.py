import h5py
import numpy as np

with h5py.File("../../datasets/siftsmall-128-euclidean.hdf5", "r") as f:
    train = f["train"][:]
    test = f["test"][:]
    neighbors = f["neighbors"][:]

print(train.shape, test.shape, neighbors.shape)
np.savetxt("../../datasets/train.txt", train, fmt = "%.6f", delimiter = " ")
np.savetxt("../../datasets/queries.txt", test, fmt = "%.6f", delimiter = " ")
np.savetxt("../../datasets/answers.txt", neighbors, fmt = "%d", delimiter = " ")