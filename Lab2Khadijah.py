# Name: Siti Khadijah binti Azhar
# Matric Number: A22EC0270

import streamlit as st
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
st.set_page_config(layout="wide")

# Function to retrieve PPI data from BioGRID
def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "9fa0658bc24af311ca0c4e0abaab38f7",
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,
        "searchbiogridids": True,
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    data = response.json()
    return pd.DataFrame.from_dict(data, orient='index')

# Function to retrieve PPI data from STRING
def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606
    }
    response = requests.get(string_url, params=params)
    data = response.json()
    return pd.json_normalize(data)

# Function that uses the PPI and creates the network using networkx library
def generate_network(dataframe):
    if "OFFICIAL_SYMBOL_A" in dataframe.columns and "OFFICIAL_SYMBOL_B" in dataframe.columns:
        network_graph = nx.from_pandas_edgelist(dataframe, "OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B")
    elif "preferredName_A" in dataframe.columns and "preferredName_B" in dataframe.columns:
        network_graph = nx.from_pandas_edgelist(dataframe, "preferredName_A", "preferredName_B")
    else:
        raise ValueError("Unexpected column names in DataFrame")
    
    return network_graph

# Function to retrieve all network centrality measures
def get_centralities(network_graph):
    degree_centrality = nx.degree_centrality(network_graph)
    closeness_centrality = nx.closeness_centrality(network_graph)
    betweenness_centrality = nx.betweenness_centrality(network_graph)
    eigenvector_centrality = "Calculation not performed"
    
    # Calculate eigenvector centrality 
    try:
        eigenvector_centrality = nx.eigenvector_centrality(network_graph, max_iter=500)
    except nx.PowerIterationFailedConvergence:
        eigenvector_centrality = "Calculation failed to converge"
    
    pagerank_centrality = nx.pagerank(network_graph)

    return {
        "Degree Centrality": degree_centrality,
        "Betweenness Centrality": betweenness_centrality,
        "Closeness Centrality": closeness_centrality,
        "Eigenvector Centrality": eigenvector_centrality,
        "PageRank Centrality": pagerank_centrality
    }

# Streamlit app
st.title("Protein-Protein Interaction (PPI) Network Viewer")

target_protein = st.text_input("Enter Protein ID")
database_choice = st.selectbox("Choose Database", ["BioGRID", "STRING"])

if st.button("Retrieve PPI Data"):
    if database_choice == "BioGRID":
        ppi_data = retrieve_ppi_biogrid(target_protein)
    else:
        ppi_data = retrieve_ppi_string(target_protein)

    if not ppi_data.empty:
        col1, col2 = st.columns(2)

        # Column 1: PPI Data Information
        with col1:
            st.subheader("PPI Data Information")
            # Display the PPI DataFrame
            st.dataframe(ppi_data.head(100))
            # Generate network graph
            network_graph = generate_network(ppi_data)
            st.write(f"Number of nodes: {network_graph.number_of_nodes()}")
            st.write(f"Number of edges: {network_graph.number_of_edges()}")
            # Plotting the network
            fig, ax = plt.subplots(figsize=(8, 6))
            pos = nx.spring_layout(network_graph, seed=42)
            nx.draw(network_graph, pos, with_labels=True, node_size=50, node_color="skyblue", font_size=8, ax=ax)
            st.pyplot(fig)

        # Column 2: Centrality Measures
        with col2:
            st.subheader("Centrality Measures")
            # Retrieve centrality measures
            centralities = get_centralities(network_graph)
            # Create and display a separate table for each centrality measure
            for centrality_name, centrality_data in centralities.items():
                st.write(centrality_name)
                centrality_df = pd.DataFrame(list(centrality_data.items()), columns=["Protein", centrality_name])
                st.dataframe(centrality_df)
