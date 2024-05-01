# import the ZINC dataset

from khc import *
import pickle

# Load the ZINC dataset
from torch_geometric.datasets import ZINC


train = ZINC(root="data/ZINC", subset=True, split="train", pre_transform=process)
val = ZINC(root="data/ZINC", subset=True, split="val", pre_transform=process)
test = ZINC(root="data/ZINC", subset=True, split="test", pre_transform=process)


# Save the dataset
with open("data/ZINC/processed.pkl", "wb") as f:
    pickle.dump((train, val, test), f)
