import types
from grnndata import GRNAnnData
import d3graph
from typing import List, Optional

import gseapy as gp
import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse
import seaborn as sns
import tqdm
from anndata import AnnData
from anndata import read_h5ad as anndata_read_h5ad
from d3graph import d3graph
from matplotlib import pyplot as plt
from pyvis import network as pnx
from sklearn.metrics.pairwise import cosine_similarity
from scprint2.tasks import GNInfer as _GNInfer

# Get the base seaborn color palette as hex list
base_color_palette = sns.color_palette().as_hex()


def customize_GRNAnnData(obj : GRNAnnData) -> GRNAnnData:

    def plot_subgraph(
        self,
        pathway_name: str,
        genes: List[str],
        gene_col: str = "symbol",
        palette: List[str] = base_color_palette,
        interactive: bool = False,
        do_enr: bool = False,
        # color_overlap: pd.DataFrame | None = None,
        color_edges : bool = True,
        node_size: int = 3000,
        font_size: int = 14,
        out_filepath: str = "",
        **kwargs: dict,
    ):
        """
        plot_subgraph plots a subgraph of the gene regulatory network (GRN) centered around a seed gene.

        Args:
            seed (str or list): The seed gene or list of genes from which the subgraph will be centered.
            gene_col (str, optional): The column name in the .var DataFrame that contains gene identifiers. Defaults to "symbol".
            max_genes (int, optional): The maximum number of genes to include in the subgraph. Defaults to 10.
            only (float, optional): The threshold for filtering connections. If less than 1, it is used as a minimum weight threshold. If 1 or greater, it is used as the number of top connections to retain. Defaults to 0.3.
            palette (list, optional): The color palette to use for plotting. Defaults to base_color_palette.
            interactive (bool, optional): Whether to create an interactive plot. Defaults to True.
            do_enr (bool, optional): Whether to perform enrichment analysis on the subgraph. Defaults to False.
            color_overlap (pd.DataFrame | None, optional): A DataFrame with geneA, geneB, and value columns to color edges based on overlap. Defaults to None.
            node_size (int, optional): Size of the nodes in the plot. Defaults to 3000.
            font_size (int, optional): Font size for the node labels. Defaults to 14.
            **kwargs: Additional keyword arguments to pass to the d3graph or networkx plotting functions.

        Returns:
            d3graph or None: The d3graph object if interactive is True, otherwise None.
        """
        
        # Gathering gene symbols
        rn = {k: v for k, v in self.var[gene_col].items()}

        available_genes = [gene for gene in genes if gene in rn.keys()]

        diff = list(set(genes).difference(set(available_genes)))

        if diff:
            print(f"Leaving behind {len(diff)} genes:")
            print(diff)
            genes = available_genes

        # Sous-matrice du GRN pour les gènes du pathway
        mat = self.grn.loc[genes, genes]

        # Rename with gene symbol instead of ensembl IDs
        mat = mat.rename(index=rn, columns=rn)

        mat = mat * 100
        color = [palette[0]] * len(mat)

        mat = mat.T

        if color_edges:
            import matplotlib.cm as cm
            import matplotlib.colors as mcolors

            # Create a mapping from (geneA, geneB) to value
            edge_colors = {}
            for geneA, row in mat.iterrows():
                for geneB, value in row.items():
                    edge_colors[(geneA, geneB)] = value if pd.notna(value) else 0

            # Get min and max values for normalization
            vmin, vmax = mat.min(axis=None), mat.max(axis=None)

            # Create viridis colormap and normalizer
            cmap = cm.get_cmap("viridis")
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # print(mat)
        # print(color)
        if interactive:
            d3 = d3graph()
            d3.graph(mat, color=None)
            d3.set_node_properties(color=color, fontcolor="#000000", **kwargs)
            d3.set_edge_properties(directed=True)
            for i, gene1 in enumerate(mat.index):
                for j, gene2 in enumerate(mat.columns):
                    if mat.iloc[i, j] != 0 and gene1 != gene2:
                        if (gene1, gene2) in edge_colors:
                            rgba = cmap(norm(edge_colors[(gene1, gene2)]))
                            hex_color = mcolors.rgb2hex(rgba)
                            d3.edge_properties[(gene1, gene2)]["color"] = hex_color
                        else:
                            d3.edge_properties[(gene1, gene2)][
                                "color"
                            ] = "#808080"  # Default gray
            return d3
        else:
            # Create a graph from the DataFrame
            G = nx.from_pandas_adjacency(mat, create_using=nx.DiGraph())
            # Draw the graph
            plt.figure(figsize=(15, 15))  # Increase the size of the plot
            # Use spring layout with adjusted parameters for better spacing
            pos = nx.spring_layout(G, k=2.0, iterations=50, seed=42)

            # Or try kamada_kawai layout for better node distribution
            # pos = nx.kamada_kawai_layout(G)

            # Draw nodes with larger size
            nx.draw_networkx_nodes(
                G, pos, node_size=node_size, node_color=color, alpha=0.6
            )

            if color_edges:
                # For networkx plotting
                edge_colors_nx = []
                for edge in G.edges():
                    if edge in edge_colors:
                        edge_colors_nx.append(cmap(norm(edge_colors[edge])))
                    elif (edge[1], edge[0]) in edge_colors:  # Check reverse direction
                        edge_colors_nx.append(
                            cmap(norm(edge_colors[(edge[1], edge[0])]))
                        )
                    else:
                        edge_colors_nx.append("#808080")  # Default gray
                # Draw edges with transparency
                nx.draw_networkx_edges(
                    G,
                    pos,
                    alpha=0.7,
                    width=3.5,
                    edge_color=edge_colors_nx,
                    arrows=True,
                    edge_cmap=cmap,
                    node_size=node_size,
                    connectionstyle="arc3,rad=0.1",
                    min_source_margin=15,
                    min_target_margin=15,
                )
                # Add colorbar
                sm = cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                plt.colorbar(sm, ax=plt.gca(), label="Overlap Value")
            else:
                nx.draw_networkx_edges(
                    G,
                    pos,
                    alpha=0.7,
                    width=2.0,
                    arrows=True,
                    node_size=node_size,
                    connectionstyle="arc3,rad=0.1",
                    min_source_margin=15,
                    min_target_margin=15,
                )
            # Draw labels with better formatting
            nx.draw_networkx_labels(
                G,
                pos,
                font_size=font_size,
                font_weight="bold",
                font_family="sans-serif",
                bbox=dict(
                    boxstyle="round,pad=0.3",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.7,
                ),
            )
            plt.axis("off")
            plt.title(pathway_name)
            plt.tight_layout()

            plt.savefig(out_filepath)
            
        if do_enr:
            enr = gp.enrichr(
                gene_list=list(G.nodes),
                gene_sets=[
                    "KEGG_2021_Human",
                    "MSigDB_Hallmark_2020",
                    "Reactome_2022",
                    "Tabula_Sapiens",
                    "WikiPathway_2023_Human",
                    "TF_Perturbations_Followed_by_Expression",
                    "Reactome",
                    "PPI_Hub_Proteins",
                    "OMIM_Disease",
                    "GO_Molecular_Function_2023",
                ],
                organism="Human",  # change accordingly
                # description='pathway',
                # cutoff=0.08, # test dataset, use lower value for real case
                background=self.var.symbol.tolist(),
            )
            print(enr.res2d.head(20))
        return G

    # Assigning custom method here
    obj.plot_subgraph = types.MethodType(plot_subgraph,
            obj)
    
    return obj

class GNInfer(_GNInfer):

    def __init__(self, *args, **kwargs):
        self.genelist = kwargs["genelist"]
        super().__init__(*args, **kwargs)
    