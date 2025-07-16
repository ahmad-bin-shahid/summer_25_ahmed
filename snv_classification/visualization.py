# visualization.py
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import torch
from torch.utils.data import DataLoader
from dataset import DNADataset
import os
import seaborn as sns
import pandas as pd

def plot_confusion_matrix(y_true, y_pred, save_path="plots/confusion_matrix_heatmap.png"):
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=["Neutral", "Deleterious"], yticklabels=["Neutral", "Deleterious"])
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Confusion Matrix (Heatmap)")
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_metric_bars(metrics_dict, save_path="plots/metric_bars.png"):
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.figure(figsize=(6, 4))
    keys = list(metrics_dict.keys())
    values = [metrics_dict[k] for k in keys]
    sns.barplot(x=keys, y=values, palette='viridis')
    plt.ylim(0, 1)
    plt.title("Evaluation Metrics")
    for i, val in enumerate(values):
        plt.text(i, val + 0.02, f"{val:.2f}", ha='center')
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def evaluate_and_visualize(model, X_test, y_test, save_dir="plots"):
    os.makedirs(save_dir, exist_ok=True)
    model.eval()
    test_data = DNADataset(X_test, y_test)
    test_loader = DataLoader(test_data, batch_size=1)

    predictions = []
    actual = []
    probabilities = []

    with torch.no_grad():
        for inputs, labels in test_loader:
            outputs = model(inputs)
            prob = torch.softmax(outputs, dim=1)
            _, predicted = torch.max(prob, 1)

            predictions.append(predicted.item())
            actual.append(labels.item())
            probabilities.append(prob[0][1].item())

    acc = accuracy_score(actual, predictions)
    prec = precision_score(actual, predictions)
    rec = recall_score(actual, predictions)
    f1 = f1_score(actual, predictions)

    print(f"Accuracy: {acc:.4f}, Precision: {prec:.4f}, Recall: {rec:.4f}, F1-Score: {f1:.4f}")

    # Plot bar graph of metrics
    metric_dict = {"Accuracy": acc, "Precision": prec, "Recall": rec, "F1-Score": f1}
    plot_metric_bars(metric_dict, save_path=os.path.join(save_dir, "metric_bars.png"))

    # Confusion matrix as heatmap
    plot_confusion_matrix(actual, predictions, save_path=os.path.join(save_dir, "confusion_matrix_heatmap.png"))