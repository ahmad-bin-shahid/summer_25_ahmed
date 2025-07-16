from training import train_model
from evaluate import model_evaluation
import pandas as pd
from extract_features import extract_features
from visualization import evaluate_and_visualize

if __name__ == "__main__":
    df = pd.read_csv("mapped_variants.csv")

    feature_dicts = extract_features(df)
    feature_df = pd.DataFrame(feature_dicts)
    feature_df.to_csv("snv_features.csv", index=False)

    print("Saved feature-enriched data to 'snv_features.csv'")
    print(feature_df.head())

    model, X_test, y_test = train_model()
    model_evaluation(model, X_test, y_test)
    evaluate_and_visualize(model, X_test, y_test)
