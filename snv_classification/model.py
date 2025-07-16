import torch
import torch.nn as nn

class simpleANN(nn.Module):
    def __init__(self, input_dim=48, output_dim=32):
        super(simpleANN, self).__init__()
        self.fc1 = nn.Linear(input_dim, output_dim)
        self.fc2 = nn.Linear(output_dim, 28)
        self.func = nn.ReLU()
        self.fc3 = nn.Linear(28, 2)

    def forward(self, x):
        x = self.fc1(x)
        x = self.func(self.fc2(x))
        x = self.fc3(x)
        return x
