import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from dataset import DNADataset
from model import simpleANN
from sklearn.model_selection import train_test_split
import pandas as pd

def train_model(max_epochs=40):
    df = pd.read_csv("snv_features.csv")
    drop_cols = [col for col in ['Chrom', 'Gene', 'Pos', 'Variant_Type'] if col in df.columns]
    feature_cols = [col for col in df.columns if col not in drop_cols + ['Label']]

    X = df[feature_cols].to_dict(orient='records')
    y = df['Label'].tolist()

    x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=32)
    train_loader = DataLoader(DNADataset(x_train, y_train), batch_size=32, shuffle=True)

    model = simpleANN(input_dim=len(feature_cols))
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    for epoch in range(max_epochs):
        model.train()
        total_loss = 0
        for inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        print(f"Epoch {epoch+1}, Train Loss: {total_loss:.4f}")

    return model, x_test, y_test
