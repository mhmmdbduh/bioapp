import streamlit as st

pages = {
    "Home": [
        st.Page("home.py", title="Home"),
        st.Page("Seqtools.py", title="Sequence Tools"),
    ]
}

pg = st.navigation(pages)
pg.run()