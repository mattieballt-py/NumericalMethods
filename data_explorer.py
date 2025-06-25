import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(page_title="Data Explorer", layout="centered")
st.title("Data Explorer")

# Sidebar controls
st.sidebar.header("Settings")
points = st.sidebar.slider("Number of data points", min_value=10, max_value=500, value=100)
seed = st.sidebar.number_input("Random seed", value=42)
generate = st.sidebar.button("Generate Data")

# Generate data
if generate:
    np.random.seed(seed)
    data = {
        "x": np.linspace(0, 10, points),
        "y": np.sin(np.linspace(0, 10, points)) + np.random.normal(0, 0.5, points)
    }
    df = pd.DataFrame(data)

    st.success("Data generated!")
    st.line_chart(df)
    st.dataframe(df.style.highlight_max(axis=0))
else:
    st.info("Use the sidebar to generate data.")
