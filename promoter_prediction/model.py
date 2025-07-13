import pandas as pd
import numpy as np
from sklearn.preprocessing import RobustScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

def preprocess_and_train_model(features_csv):
    """Preprocess features, apply PCA if needed, and train a regression model."""
    # Load data
    df = pd.read_csv(features_csv)
    X = df.drop(['gene', 'activity'], axis=1)
    y = df['activity']
    
    # Scale features
    scaler = RobustScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Apply PCA if number of features is large (>100)
    if X_scaled.shape[1] > 100:
        pca = PCA(n_components=50)  # Reduce to 50 components
        X_scaled = pca.fit_transform(X_scaled)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)
    
    # Train Random Forest Regressor
    model = RandomForestRegressor(n_estimators=200, max_depth=10, random_state=42)
    model.fit(X_train, y_train)
    
    # Evaluate model
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    return mse, r2