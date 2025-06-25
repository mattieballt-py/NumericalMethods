import streamlit as st 
import pandas as pd 
import numpy as np 

st.set_page_config(page_title="Data Explorer", layout="centered")
st.title("Data Explorer")

# sidebar controls
st.sidebar.header("settings")
points = st.sidebar.slider("no. data points")
seed = st.sidebar.number_input("random seed", value=42)
generate = st.sidebar.button("gen data")

# gen data
if generate:
    np.random.seed(seed)
    data = {
        "x": np.linspace(0,10,points)
        "y": np.sin(np.linspace(0,10,points)) + np.random.normal(0,0.5,points)
    }
    df = pd.DataFrame(data)

    st.success("data gened")
    st.line_chart(df)
    st.dataframe(df.style.highlight_max(axis=0))
else:
    st.info("use sidebar to gen data")