from dataset import DNADataset
from torch.utils.data import DataLoader
import torch
from sklearn.metrics import accuracy_score

def model_evaluation(model, X_test, y_test):
    model.eval()
    test_data = DNADataset(X_test, y_test)
    test_loader = DataLoader(test_data, batch_size=1)

    predictions = []
    actual = []

    with torch.no_grad():
        for inputs, labels in test_loader:
            outputs = model(inputs)
            _, predicted = torch.max(outputs.data, 1)
            predictions.append(predicted.item())
            actual.append(labels.item())

    acc = accuracy_score(actual, predictions)
    print(f"Test Accuracy: {acc * 100:.2f}%")
