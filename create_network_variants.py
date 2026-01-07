"""
Create Network Variants for Mahkota Dewa - DN Analysis
"""
import pandas as pd
import requests
import json

# Load current network
edges = pd.read_csv(r'data\mahkota_dewa_dn\results\network_edges.csv')
centralities = pd.read_csv(r'data\mahkota_dewa_dn\results\network_centralities.csv')

print("=" * 50)
print("OPSI 2: Filter nodes with degree >= 4 (remove isolated)")
print("=" * 50)

# Filter nodes with degree >= 4
high_degree_genes = centralities[centralities['degree'] >= 4]['gene'].tolist()
print(f"Genes with degree >= 4: {len(high_degree_genes)}")
print(high_degree_genes)

# Filter edges
edges_filtered = edges[edges['source'].isin(high_degree_genes) & edges['target'].isin(high_degree_genes)]
print(f"Edges after filter: {len(edges_filtered)}")

# Save opsi 2
edges_filtered.to_csv(r'data\mahkota_dewa_dn\results\network_opsi2_filtered.csv', index=False)

# SIF format
sif_lines = []
for _, r in edges_filtered.iterrows():
    sif_lines.append(f"{r['source']}\tpp\t{r['target']}")
with open(r'data\mahkota_dewa_dn\results\network_opsi2_filtered.sif', 'w') as f:
    f.write('\n'.join(sif_lines))

print("Opsi 2 saved!")
print()

print("=" * 50)
print("OPSI 3: Add first neighbors from STRING")
print("=" * 50)

# Get intersection genes
with open(r'data\mahkota_dewa_dn\processed\intersection_genes.txt') as f:
    intersection_genes = [g.strip() for g in f.readlines() if g.strip()]

print(f"Intersection genes: {len(intersection_genes)}")

# Query STRING with add_nodes parameter for first neighbors
STRING_API = "https://string-db.org/api"
genes_str = "%0d".join(intersection_genes)

# Get network with first neighbors
params = {
    "identifiers": genes_str,
    "species": 9606,  # Human
    "add_nodes": 10,  # Add 10 first neighbors
    "required_score": 400,  # Medium confidence
}

try:
    response = requests.get(
        f"{STRING_API}/json/network",
        params=params,
        timeout=60
    )
    
    if response.status_code == 200:
        data = response.json()
        print(f"STRING returned {len(data)} interactions")
        
        # Build edges
        edges_opsi3 = []
        nodes_set = set()
        for item in data:
            source = item.get('preferredName_A', item.get('stringId_A', ''))
            target = item.get('preferredName_B', item.get('stringId_B', ''))
            score = item.get('score', 0)
            
            if source and target:
                edges_opsi3.append({
                    'source': source,
                    'target': target,
                    'weight': score
                })
                nodes_set.add(source)
                nodes_set.add(target)
        
        print(f"Nodes in network: {len(nodes_set)}")
        print(f"Edges in network: {len(edges_opsi3)}")
        
        # Save
        df_opsi3 = pd.DataFrame(edges_opsi3)
        df_opsi3.to_csv(r'data\mahkota_dewa_dn\results\network_opsi3_expanded.csv', index=False)
        
        # SIF format
        sif_lines3 = []
        for e in edges_opsi3:
            sif_lines3.append(f"{e['source']}\tpp\t{e['target']}")
        with open(r'data\mahkota_dewa_dn\results\network_opsi3_expanded.sif', 'w') as f:
            f.write('\n'.join(sif_lines3))
        
        # Save node types
        node_types = []
        for node in nodes_set:
            if node in intersection_genes:
                node_types.append({'gene': node, 'type': 'intersection'})
            else:
                node_types.append({'gene': node, 'type': 'first_neighbor'})
        
        pd.DataFrame(node_types).to_csv(r'data\mahkota_dewa_dn\results\network_opsi3_node_types.csv', index=False)
        
        print("Opsi 3 saved!")
    else:
        print(f"STRING API error: {response.status_code}")
        
except Exception as e:
    print(f"Error: {e}")

print()
print("=" * 50)
print("DONE! Files created:")
print("=" * 50)
print("Opsi 2: network_opsi2_filtered.csv, network_opsi2_filtered.sif")
print("Opsi 3: network_opsi3_expanded.csv, network_opsi3_expanded.sif, network_opsi3_node_types.csv")
