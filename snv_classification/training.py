import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from dataset import DNADataset
from model import simpleANN
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np

def train_model(patience=5, max_epochs=100):
    df = pd.read_csv("snv_features.csv")
    
    # Drop unnecessary columns
    drop_cols = [col for col in ['Chrom', 'Gene', 'Pos', 'Variant_Type'] if col in df.columns]
    feature_cols = [col for col in df.columns if col not in drop_cols + ['Label']]

    X = df[feature_cols].to_dict(orient='records')
    y = df['Label'].tolist()

    x_train_full, x_test, y_train_full, y_test = train_test_split(X, y, test_size=0.2, random_state=32)
    x_train, x_val, y_train, y_val = train_test_split(x_train_full, y_train_full, test_size=0.2, random_state=32)

    train_loader = DataLoader(DNADataset(x_train, y_train), batch_size=32, shuffle=True)
    val_loader = DataLoader(DNADataset(x_val, y_val), batch_size=32)

    model = simpleANN(input_dim=len(feature_cols))
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    best_val_loss = np.inf
    best_model_state = model.state_dict()
    patience_counter = 0

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

        model.eval()
        val_loss = 0
        with torch.no_grad():
            for inputs, labels in val_loader:
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                val_loss += loss.item()

        avg_val_loss = val_loss / len(val_loader)
        print(f"Epoch {epoch+1}, Train Loss: {total_loss:.4f}, Val Loss: {avg_val_loss:.4f}")

        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            best_model_state = model.state_dict()
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print("Early stopping triggered!")
                break

    model.load_state_dict(best_model_state)
    return model, x_test, y_test
