import torch
from torch.utils.data import Dataset

class DNADataset(Dataset):
    def __init__(self, features, labels):
        self.features = features 
        self.labels = labels

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        feature_vector = list(self.features[idx].values())
        return torch.tensor(feature_vector, dtype=torch.float32), torch.tensor(self.labels[idx], dtype=torch.long)
