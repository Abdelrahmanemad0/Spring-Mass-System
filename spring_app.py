"""spring_app.py -- Interactive Streamlit demo for the two-mass spring
(train suspension) numerical simulation. Lets you tune the physical
parameters and compares Euler, Heun, and RK4 integration live.

Run:
    streamlit run spring_app.py
"""
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from spring_mass_sim import SimParams, run_all

st.set_page_config(page_title="Spring-Mass System Simulator", page_icon="\U0001F683", layout="wide")

st.title("Train Suspension Spring-Mass System")
st.caption(
    "Two-mass, two-spring numerical simulation comparing Euler, Heun, and "
    "4th-order Runge-Kutta integration. Adjust the parameters in the sidebar "
    "and re-run."
)

with st.sidebar:
    st.header("Parameters")
    k1 = st.slider("Spring constant k1 (N/m)", 0.5, 20.0, 5.0, 0.5)
    k2 = st.slider("Spring constant k2 (N/m)", 0.5, 20.0, 5.0, 0.5)
    m1 = st.slider("Mass 1 (kg)", 0.5, 20.0, 2.0, 0.5)
    m2 = st.slider("Mass 2 (kg)", 0.5, 20.0, 2.0, 0.5)
    L1 = st.slider("Unstretched length L1 (m)", 0.5, 10.0, 2.0, 0.5)
    L2 = st.slider("Unstretched length L2 (m)", 0.5, 10.0, 2.0, 0.5)
    w1 = st.slider("Width of mass 1 (m)", 0.5, 10.0, 5.0, 0.5)
    h = st.select_slider("Step size h (s)", options=[0.001, 0.005, 0.01, 0.05, 0.1], value=0.01)
    t_final = st.slider("Simulation time (s)", 5, 60, 20, 5)
    run_button = st.button("Run simulation", type="primary")

if "results" not in st.session_state or run_button:
    params = SimParams(k1=k1, k2=k2, m1=m1, m2=m2, L1=L1, L2=L2, w1=w1, h=h, t_final=float(t_final))
    st.session_state["results"] = run_all(params)

results = st.session_state["results"]

colors = {"euler": "#EF553B", "heun": "#00CC96", "rk4": "#636EFA"}
labels = {"euler": "Euler", "heun": "Heun", "rk4": "Runge-Kutta 4"}

tab_disp, tab_vel, tab_error = st.tabs(["Displacement", "Velocity", "Approximate Error"])

with tab_disp:
    col1, col2 = st.columns(2)
    for col, var, title in ((col1, "x1", "Displacement of Mass 1"), (col2, "x2", "Displacement of Mass 2")):
        fig = go.Figure()
        for method in ("euler", "heun", "rk4"):
            fig.add_trace(go.Scatter(
                x=results[method]["t"], y=results[method][var],
                mode="lines", name=labels[method], line=dict(color=colors[method]),
            ))
        fig.update_layout(title=title, xaxis_title="Time (s)", yaxis_title=f"{var} (m)", height=400)
        col.plotly_chart(fig, use_container_width=True)

with tab_vel:
    col1, col2 = st.columns(2)
    for col, var, title in ((col1, "v1", "Velocity of Mass 1"), (col2, "v2", "Velocity of Mass 2")):
        fig = go.Figure()
        for method in ("euler", "heun", "rk4"):
            fig.add_trace(go.Scatter(
                x=results[method]["t"], y=results[method][var],
                mode="lines", name=labels[method], line=dict(color=colors[method]),
            ))
        fig.update_layout(title=title, xaxis_title="Time (s)", yaxis_title=f"{var} (m/s)", height=400)
        col.plotly_chart(fig, use_container_width=True)

with tab_error:
    st.subheader("Final-step approximate percent error")
    rows = []
    for method in ("euler", "heun", "rk4"):
        r = results[method]
        rows.append({
            "Method": labels[method],
            "x1 (m)": r["x1"][-1], "x2 (m)": r["x2"][-1],
            "v1 (m/s)": r["v1"][-1], "v2 (m/s)": r["v2"][-1],
            "Ea x1 (%)": r["err_x1"][-1], "Ea x2 (%)": r["err_x2"][-1],
            "Ea v1 (%)": r["err_v1"][-1], "Ea v2 (%)": r["err_v2"][-1],
        })
    st.dataframe(pd.DataFrame(rows).set_index("Method"), use_container_width=True)
    st.caption(
        "Ea = |x(i+1) - x(i)| / x(i+1) * 100 — the approximate relative error "
        "between consecutive time steps, evaluated at the final step. RK4 "
        "converges fastest, so its curve typically overlaps Heun's closely "
        "while Euler drifts more visibly at larger step sizes."
    )
