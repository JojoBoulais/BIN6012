import json
from posixpath import sep
import re
import pandas as pd

from utils_filepath import *
from grnndata import GRNAnnData
from customize_libs import customize_GRNAnnData
from anndata import AnnData
import bionty as bt
import sqlite3

import lamindb as ln


class PathwayLoader:
    _instance = None

    def __new__(cls, version=3):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._version = version
            cls._instance._filepath = get_filepath(
                "data", f"pathways_v{version}.json")
            cls._instance._pathways = None
        return cls._instance

    def __init__(self, version=3):
        if self._pathways is None:
            with open(self._filepath, "r") as f_pathways:
                self._pathways = json.load(f_pathways)

        self.pathway_names = list(self._pathways.keys())

    def get_pathways(self):
        return self.pathway_names

    def load_pathway_by_name(self, name):
        assert name in self._pathways
        return Pathway(name, self._pathways[name])

    def load_all_pathways(self):
        return [Pathway(name, self._pathways[name]) for name in self.pathway_names]


class Pathway:

    _pathways = None
    _re_human = re.compile(r"^ENSG")
    _re_mus = re.compile(r"^ENSMUS")

    def __init__(self, name, json_pathway):
        self.name = name
        self.pathway = json_pathway
        self.genes = list(self.pathway.keys())

    def get_genes_ensembl(self, organism="human"):
        genes = []

        for gene, ensembl_ids in self.pathway.items():
            genes.extend(ensembl_ids)

        return genes

    def get_name_as_path(self):
        return self.name.replace(" ", "_").replace("/", "_")

    def get_name(self):
        return self.name


if __name__ == "__main__":
    # ========================== PATHWAYS ==========================
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8913803/

    # PPI (ensembl)
    # https://string-db.org/cgi/download?sessionId=bZwq8wTEfbPY&species_text=Homo+sapiens

    # inflammation response
    NLRP3_caspase1 = {"NLRP3": ["ENSG00000162711"],
                      "caspase-1": ["ENSG00000137752"],
                      "IL1B": ["ENSG00000125538"],
                      "IL18": ["ENSG00000150782"]
                      }

    TLR4_MyD88_NFkB = {
        "TLR2": ["ENSG00000137462"],
        "TLR4": ["ENSG00000136869"],
        "MyD88": ["ENSG00000172936"],
        "NFKB1": ["ENSG00000109320"],
        "NFKB2": ["ENSG00000077150"],
        "IL1A": ["ENSG00000115008"],
        "IL1R1": ["ENSG00000115594"],
        "IRAK4": ["ENSG00000198001"],
        "TRAF6": ["ENSG00000175104"],
        "S100A1": ["ENSG00000160678"],
        "S100A8": ["ENSG00000143546"],
        "LGALS3": ["ENSG00000131981"],
        "S100B": ["ENSG00000160307"]
    }

    # oxidative stress and apoptosis
    Notch = {
        "NOTCH1": ["ENSG00000148400"],
        "NFKB1": ["ENSG00000109320"],
        "NOTCH2": ["ENSG00000134250"],
        "Delta-like 1": ["ENSG00000275555", "ENSG00000198719"],
        "GSK3B": ["ENSG00000082701"],
        "JAG1": ["ENSG00000101384"],
        "JAG2": ["ENSG00000184916"],
        "NOTCH3": ["ENSG00000074181"],
        "NOTCH4": ["ENSG00000234876"],
    }

    Hippo = {
        "DCHS1": ["ENSG00000166341"],
        "DCHS2": ["ENSG00000197410"],
        "STK4": ["ENSG00000101109"],  # Hippo
        "STK3": ["ENSG00000104375"],
        "SAV1": ["ENSG00000151748"],
        "Yap": ["ENSG00000137693"],
        "WWTR1": ["ENSG00000018408"],
        "TEAD1": ["ENSG00000187079"],
        "TEAD2": ["ENSG00000074219"],
        "TEAD3": ["ENSG00000007866"],
        "TEAD4": ["ENSG00000197905"]
    }

    Nrf2_HO_1 = {
        "NFE2L2": ["ENSG00000116044"],
        "NRF2": ["ENSG00000116044"],
        "KEAP1": ["ENSG00000079999"],
        "HMOX1": ["ENSG00000100292"],  # HO-1
        "HMOX2": ["ENSG00000103415"],
        "GSTA1": ["ENSG00000243955"],
        "GSTP1": ["ENSG00000084207"],
        "NQO1": ["ENSG00000181019"]
    }

    RhoA_ROCK = {
        "ROCK1": ["ENSG00000067900"],
        "ROCK2": ["ENSG00000134318"],
        "RHOA": ["ENSG00000067560"],
        "ISL1": ["ENSG00000016082"],
        "CFL1": ["ENSG00000172757"],
        "PTEN": ["ENSG00000171862", "ENSG00000284792"],
        "DPYSL2": ["ENSG00000092964"]  # CRMP-2
    }

    MAPK = {
        "MAPK1": ["ENSG00000100030"],
        "MAPK3": ["ENSG00000102882"],
        "MAPK7": ["ENSG00000166484"],  # ERK5
        "MAPK8": ["ENSG00000107643"],  # JNK1
        "MAPK9": ["ENSG00000050748"],  # JNK2
        "MAPK10": ["ENSG00000109339"],  # JNK3
        "MAPK11": ["ENSG00000185386"],
        "ACKR3": ["ENSG00000144476"],
        "HSP90AA1": ["ENSG00000080824"],
        "ADRA1A": ["ENSG00000120907"],
        "MST1": ["ENSG00000173531"],
        "RAF1": ["ENSG00000132155"],
        "RGS5": ["ENSG00000143248"],
        "ANO1": ["ENSG00000131620"],
        "RGS5": ["ENSG00000143248"],
        "PTGS2": ["ENSG00000073756"]
    }

    # main regulator of angiogenesis
    PI3K_Akt = {
        # PI3K https://en.wikipedia.org/wiki/P110%CE%B1
        "PIK3CA": ["ENSG00000121879"],
        "PIK3CB": ["ENSG00000051382"],
        "AKT1": ["ENSG00000142208"],
        "AKT2": ["ENSG00000105221"],
        "AKT3": ["ENSG00000117020", "ENSG00000275199"],
        "VEGFA": ["ENSG00000112715"],
        "NOS3": ["ENSG00000164867"],  # eNOS
        "MTOR": ["ENSG00000198793"],
        "GSK3B": ["ENSG00000082701"],
        "FOXO1": ["ENSG00000150907"],
        "MST1": ["ENSG00000173531"],
        "EPO": ["ENSG00000130427"],
        "NOS2": ["ENSG00000007171"],
        "VEGFA": ["ENSG00000112715"]
    }

    JAK_STAT = {
        "STAT1": ["ENSG00000115415"],
        "STAT3": ["ENSG00000168610"],
        "JAK1": ["ENSG00000162434"],
        "FOS": ["ENSG00000170345"],
        "NFKB1": ["ENSG00000109320"],
        "NFKB2": ["ENSG00000077150"],
        "CASP1": ["ENSG00000137752"],
        "BCL2": ["ENSG00000171791"],
        "BAX": ["ENSG00000087088"]
    }

    TGF_SMADs = {  # TGF-β/SMADs
        "TGFB1": ["ENSG00000105329"],
        "TGFBR1": ["ENSG00000106799"],
        "TGFBR2": ["ENSG00000163513"],
        "SMAD2": ["ENSG00000175387"],
        "SMAD3": ["ENSG00000166949"],
        "SMAD4": ["ENSG00000141646"],
        "SMAD7": ["ENSG00000101665"],
        "KLF5": ["ENSG00000102554"],
        "CYTL1": ["ENSG00000170891"],
        "ANO1": ["ENSG00000131620"]
    }

    Wnt_catenin = {  # Wnt/β-catenin
        "WNT1": ["ENSG00000125084"],
        "WNT2": ["ENSG00000105989"],
        "WNT3": ["ENSG00000108379"],
        "WNT4": ["ENSG00000162552"],
        "WNT5A": ["ENSG00000114251"],
        "IL1B": ["ENSG00000125538"],
        "IL6": ["ENSG00000136244"],
        "WIF1": ["ENSG00000156076"],
        "WNT10B": ["ENSG00000169884"],
        "SFRP1": ["ENSG00000104332"],
        "MIR145": ["ENSG00000276365"],
        "WNT11": ["ENSG00000085741"],
        "LEF1": ["ENSG00000138795"],
        "LRP5": ["ENSG00000162337"],
        "LRP6": ["ENSG00000070018"],
        "CTNNB1": ["ENSG00000168036"],
        "APC": ["ENSG00000134982"],
        "GSK3B": ["ENSG00000082701"],
        "DKK2": ["ENSG00000155011"]
    }

    #  !!!!! Sonic hedgehog not done !!!!!
    sonic_hedgehog = {
        "SHH": ["ENSG00000164690"],
        "PTCH1": ["ENSG00000185920"],
        "SMO": ["ENSG00000128602"],
        "GLI1": ["ENSG00000111087"],
        "GLI2": ["ENSG00000074047"],
        "GLI3": ["ENSG00000106571"]
    }

    pathways = {"NLRP3/caspase1": NLRP3_caspase1,
                "TLR4/MyD88/NFkB": TLR4_MyD88_NFkB,
                "Nrf2/HO-1": Nrf2_HO_1,
                "Notch": Notch,
                "Hippo": Hippo,
                "RhoA/ROCK": RhoA_ROCK,
                "MAPK": MAPK,
                "PI3K/Akt": PI3K_Akt,
                "JAK/STAT": JAK_STAT,
                "TGF-β/SMADs": TGF_SMADs,
                "Wnt/β-catenin": Wnt_catenin,
                "Sonic/Hedgehog": sonic_hedgehog
                }

    def generate_json_file(pathways, version):
        fp_pathways = get_filepath("data", f"pathways_v{str(version)}.json")

        with open(fp_pathways, "w") as f_pathways:
            json.dump(pathways, f_pathways)

    # generate_json_file(pathways, 3)

    import gseapy as gp

    reactome = gp.get_library(name="Reactome_2022")

    fp = get_filepath("tmp", "pathways.txt")

    with open(fp, "w") as f:

        hearts = [k for k in reactome.keys() if "cardiac" in k.lower()
                  or "heart" in k.lower()]

        print(hearts, file=f)

        for pathway, genes in reactome.items():
            print(pathway, len(genes), file=f)

    reactome_pathways = {
        "TLR4/MyD88/NFkB_reactome": "MyD88-independent TLR4 Cascade R-HSA-166166",
        "Notch1_reactome": "NOTCH1 Intracellular Domain Regulates Transcription R-HSA-2122947",
        "Notch2_1_reactome": "NOTCH2 Activation And Transmission Of Signal To Nucleus R-HSA-2979096",
        "Notch2_2_reactome": "NOTCH2 Intracellular Domain Regulates Transcription R-HSA-2197563",
        "Notch3_1_reactome": "NOTCH3 Activation And Transmission Of Signal To Nucleus R-HSA-9013507",
        "Notch3_2_reactome": "NOTCH3 Intracellular Domain Regulates Transcription R-HSA-9013508",
        "Notch4_1_reactome": "NOTCH4 Activation And Transmission Of Signal To Nucleus R-HSA-9013700",
        "Notch4_2_reactome": "NOTCH4 Intracellular Domain Regulates Transcription R-HSA-9013695",
        "Hippo_reactome": "Signaling By Hippo R-HSA-2028269",
        "RhoA/ROCK_reactome": "RHO GTPases Activate ROCKs R-HSA-5627117",
        "MAPK1_reactome": "MAPK Targets/ Nuclear Events Mediated By MAP Kinases R-HSA-450282",
        "MAPK2_reactome": "MAPK Family Signaling Cascades R-HSA-5683057",
        "PI3K/Akt_reactome": "Negative Regulation Of PI3K/AKT Network R-HSA-199418",
        "JAK/STAT_reactome": "Gene And Protein Expression By JAK-STAT Signaling After Interleukin-12 Stimulation R-HSA-8950505",
        "TGF-β/SMADs_reactome": "TGF-beta Receptor Signaling Activates SMADs R-HSA-2173789",
        "Wnt/β-catenin_reactome": "Beta-catenin Independent WNT Signaling R-HSA-3858494"
    }

    # conn = sqlite3.connect("/home/jordboul/scratch/scPrint2/BIN6012/testdb/.lamindb/lamin.db")
    # cur = conn.cursor()

    # cur.execute("SELECT * FROM bionty_gene;")
    # genes = cur.fetchall()

    # conn.close()

    genes_fp = get_filepath("data", "genes-v48-hg38.bed")

    genes_df = pd.read_csv(genes_fp,
                           sep="\t",
                           header=None)

    genes_df.columns = ["chr", "start", "end", "symbol", "ensembl"]

    gene_map = genes_df.set_index("symbol")["ensembl"].to_dict()

    out_pathway = {}

    for pathway, pathway_reactome in reactome_pathways.items():
        out_pathway[pathway] = {}

        pt_reactome_genes = reactome[pathway_reactome]

        for gene in pt_reactome_genes:
            if gene in gene_map.keys():
                ensmbl = gene_map[gene]
                if ensmbl and Pathway.re_human.match(ensmbl):
                    out_pathway[pathway][gene] = [ensmbl]

    generate_json_file(out_pathway, 4)

    # print(len(sample_rows))
    # print(sample_rows[0])

    # 2, 4

    # source = bt.Source.get(entity="bionty.Gene", source="cl", version="2022-08-16")
    # default_sources = bt.Source.filter(entity="bionty.CellType", currently_used=True).to_dataframe()
    # genes = bt.Gene.import_source("ensembl", organism="human")

    # source = bt.Source.get(
    #     entity="bionty.Gene",
    #     name="ensembl",
    #     organism="human"
    # )

    # genes = bt.Gene.import_source(source)

    # source = bt.Source.get(
    # entity="bionty.Gene",
    # name="ensembl",
    # organism="human",
    # )

    # genes = bt.Gene.public(source=source)
    # print(type(genes))

    # print(list(genes))

    # bionty_gene
    # bionty_organism
    # bionty_organism_parents
    # bionty_pathway
    # bionty_pathway_genes
    # bionty_pathway_parents

    # print(genes)

    # print(reactome[reactome_pathways[0]])

    # anndata = customize_GRNAnnData(GRNAnnData())

    # print(anndata.rn)

    # obj.gene_col = "symbol"
    # obj.available_genes = []
    # obj.rn = {k: v for k, v in obj.var[obj.gene_col].items()}

    # print(AnnData().var)

    # dataset_fp = get_filepath("data", "myocardial_infarction_preprocessed_4.h5ad")
    # model_fp = get_filepath("model", "small-v2.ckpt")

    # model = load_model_with_cuda_if_avail(model_fp)
    # adata = scanpy.read_h5ad(dataset_fp)

    # for name, pathway in pathways.items():
    #     print(f"=========== {name} ===========")
    #     for gene, ls_ensembl in pathway.items():
    #         print(gene, ls_ensembl[0], ls_ensembl[0] in model.genes)
