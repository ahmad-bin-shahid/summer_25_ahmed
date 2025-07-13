from extract_features import extract_features
from model import preprocess_and_train_model

def main():
    # Step 1: Extract features from CSV file
    print("Extracting features...")
    features_df = extract_features("promoter_activity_dataset.csv", "features.csv")
    print("Features extracted and saved to features.csv")
    
   
    print("Training model...")
    mse, r2 = preprocess_and_train_model("features.csv")
    print("Model trained!\n")
    print("Mean Squared Error:", mse)
    print("R² Score:", r2)
    
if __name__ == "__main__":
    main()