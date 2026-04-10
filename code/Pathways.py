import json
import re

from utils_filepath import *

class Pathway():

    pathways = ["NLRP3/caspase1",
                "TLR4/MyD88/NFkB",
                "Nrf2/HO-1",
                "Notch",
                "Hippo",
                "RhoA/ROCK",
                "MAPK",
                "PI3K/Akt",
                "JAK/STAT",
                "TGF-β/SMADs",
                "Wnt/β-catenin"]

    re_human = re.compile(r"^ENSG")
    re_mus = re.compile(r"^ENSMUS")

    def __init__(self, pathway_name, version=2):
        assert pathway_name in self.pathways

        fp_pathways = get_filepath("data", f"pathways_v{version}.json")

        with open(fp_pathways, "r") as f_pathways:
            _pathways = json.load(f_pathways)

        self.pathway = _pathways[pathway_name]
        self.genes = list(self.pathway.keys())

    def get_genes_ensembl(self, organism="human"):
        genes = []

        if organism == "human":
            pattern = self.re_human
        elif organism in ["mouse", "mus", "murine"]:
            pattern = self.re_mus
        else:
            raise ValueError("organism must be 'human' or 'mouse'")

        for gene, ensembl_ids in self.pathway.items():
            genes.extend([ensembl for ensembl in ensembl_ids if pattern.match(ensembl)])

        return genes

    @classmethod
    def load_all_pathways(cls):
        return [Pathway(pathway) for pathway in cls.pathways]

    # @classmethod
    # def extract_pathway_from_filepath(cls, filepath):

    #     for pathway in cls.pathways():
    #         if pathway in filepath:
    #             return pathway
            
    #     print("No pathway found")
    #     return ""


if __name__ == "__main__":
    # ========================== PATHWAYS ==========================
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8913803/

    # PPI (ensembl)
    # https://string-db.org/cgi/download?sessionId=bZwq8wTEfbPY&species_text=Homo+sapiens

    # inflammation response
    NLRP3_caspase1 = {  "NLRP3" : ["ENSG00000162711", "ENSMUSG00000032691"],
                        "caspase-1" : ["ENSG00000137752", "ENSMUSG00000025888"],
                        "IL1B" : ["ENSG00000125538"],
                        "IL18" : ["ENSG00000150782"]
                    }

    TLR4_MyD88_NFkB = {
        "TLR2" : ["ENSG00000137462"],
        "TLR4" : ["ENSG00000136869", "ENSMUSG00000039005"],
        "MyD88" : ["ENSG00000172936", "ENSMUSG00000032508"],
        "NFKB1" : ["ENSG00000109320", "ENSMUSG00000028163"],
        "NFKB2" : ["ENSG00000077150", "ENSMUSG00000025225"],
        "IL1A" : ["ENSG00000115008"],
        "IL1R1" : ["ENSG00000115594"],
        "IRAK4" : ["ENSG00000198001"],
        "TRAF6" : ["ENSG00000175104"],
        "S100A1" : ["ENSG00000160678"],
        "S100A8" : ["ENSG00000143546"],
        "LGALS3" : ["ENSG00000131981"],
        "S100B" : ["ENSG00000160307"]
    }

    # oxidative stress and apoptosis
    Notch = {
        "NOTCH1" : ["ENSG00000148400", "ENSMUSG00000026923"],
        "NFKB1" : ["ENSG00000109320", "ENSMUSG00000028163"],
        "NOTCH2" : ["ENSG00000134250", "ENSMUSG00000027878"],
        "Delta-like 1": ["ENSG00000275555", "ENSG00000198719", "ENSMUSG00000014773"],
        "GSK3B" : ["ENSG00000082701", "ENSMUSG00000022812"],
        "JAG1" : ["ENSG00000101384", "ENSMUSG00000027276"],
        "JAG2" : ["ENSG00000184916", "ENSMUSG00000002799"],
        "NOTCH3" : ["ENSG00000074181", "ENSMUSG00000038146"],
        "NOTCH4" : ["ENSG00000234876", "ENSMUSG00000015468"],
    }

    Hippo = {
        "DCHS1" : ["ENSG00000166341", "ENSMUSG00000036862"],
        "DCHS2" : ["ENSG00000197410", "ENSG00000284227"],
        "STK4" : ["ENSG00000101109", "ENSMUSG00000018209"], # Hippo
        "STK3" : ["ENSG00000104375", "ENSMUSG00000022329"],
        "SAV1" : ["ENSG00000151748", "ENSMUSG00000021067"],
        "Yap" : ["ENSG00000137693", "ENSMUSG00000053110"],
        "WWTR1" : ["ENSG00000018408", "ENSMUSG00000027803"],
        "TEAD1" : ["ENSG00000187079", "ENSMUSG00000055320"],
        "TEAD2" : ["ENSG00000074219", "ENSMUSG00000030796"],
        "TEAD3" : ["ENSG00000007866", "ENSMUSG00000002249"],
        "TEAD4" : ["ENSG00000197905", "ENSMUSG00000030353"]
    }

    Nrf2_HO_1 = {
        "NFE2L2" : ["ENSG00000116044"],
        "NRF2" : ["ENSG00000116044", "ENSMUSG00000015839"],
        "KEAP1" : ["ENSG00000079999", "ENSMUSG00000003308"],
        "HMOX1" : ["ENSG00000100292", "ENSMUSG00000005413"], #HO-1
        "HMOX2" : ["ENSG00000103415"],
        "GSTA1" : ["ENSG00000243955", "ENSMUSG00000057933"],
        "GSTP1" : ["ENSG00000084207"],
        "NQO1" : ["ENSG00000181019", "ENSMUSG00000003849"]
    }

    RhoA_ROCK = {
        "ROCK1" : ["ENSG00000067900", "ENSMUSG00000024290"],
        "ROCK2" : ["ENSG00000134318", "ENSMUSG00000020580"],
        "RHOA" : ["ENSG00000067560", "ENSMUSG00000007815"],
        "ISL1" : ["ENSG00000016082", "ENSMUSG00000042258"],
        "CFL1" : ["ENSG00000172757", "ENSMUSG00000056201"],
        "PTEN" : ["ENSG00000171862", "ENSG00000284792", "ENSMUSG00000013663"],
        "DPYSL2" : ["ENSG00000092964", "ENSMUSG00000022048"] # CRMP-2
    }

    MAPK = {
        "MAPK1" : ["ENSG00000100030", "ENSMUSG00000063358"],
        "MAPK3" : ["ENSG00000102882", "ENSMUSG00000063065"],
        "MAPK7" : ["ENSG00000166484", "ENSMUSG00000001034"], # ERK5
        "MAPK8" : ["ENSG00000107643", "ENSMUSG00000021936"], # JNK1
        "MAPK9" : ["ENSG00000050748", "ENSMUSG00000020366"],  # JNK2
        "MAPK10" : ["ENSG00000109339", "ENSMUSG00000046709"],  # JNK3
        "MAPK11" : ["ENSG00000185386", "ENSMUSG00000053137"],
        "ACKR3" : ["ENSG00000144476", "ENSMUSG00000044337"],
        "HSP90AA1" : ["ENSG00000080824"],
        "ADRA1A" : ["ENSG00000120907"],
        "MST1" : ["ENSG00000173531"],
        "RAF1" : ["ENSG00000132155"],
        "RGS5" : ["ENSG00000143248"],
        "ANO1" : ["ENSG00000131620"],
        "RGS5" : ["ENSG00000143248"],
        "PTGS2" : ["ENSG00000073756"]
    }

    # main regulator of angiogenesis
    PI3K_Akt = {
        "PIK3CA" : ["ENSG00000121879", "ENSMUSG00000027665"], # PI3K https://en.wikipedia.org/wiki/P110%CE%B1
        "PIK3CB" : ["ENSG00000051382", "ENSMUSG00000032462"],
        "AKT1" : ["ENSG00000142208", "ENSMUSG00000001729"],
        "AKT2" : ["ENSG00000105221", "ENSMUSG00000004056"],
        "AKT3" : ["ENSG00000117020", "ENSG00000275199", "ENSMUSG00000019699"],
        "VEGFA" : ["ENSG00000112715", "ENSMUSG00000023951"],
        "NOS3" : ["ENSG00000164867", "ENSMUSG00000028978"], # eNOS
        "MTOR" : ["ENSG00000198793", "ENSMUSG00000028991"],
        "GSK3B" : ["ENSG00000082701", "ENSMUSG00000022812"],
        "FOXO1" : ["ENSG00000150907"],
        "MST1" : ["ENSG00000173531", "ENSMUSG00000032591"],
        "EPO" : ["ENSG00000130427", "ENSMUSG00000029711"],
        "NOS2" : ["ENSG00000007171"],
        "VEGFA" : ["ENSG00000112715"]
    }

    JAK_STAT = {
        "STAT1" : ["ENSG00000115415", "ENSMUSG00000026104"],
        "STAT3" : ["ENSG00000168610", "ENSMUSG00000004040"],
        "JAK1" : ["ENSG00000162434", "ENSMUSG00000028530"],
        "FOS" : ["ENSG00000170345", "ENSMUSG00000021250"],
        "NFKB1": ["ENSG00000109320", "ENSMUSG00000028163"],
        "NFKB2" : ["ENSG00000077150", "ENSMUSG00000025225"],
        "CASP1" : ["ENSG00000137752", "ENSMUSG00000025888"],
        "BCL2" : ["ENSG00000171791", "ENSMUSG00000057329"],
        "BAX" : ["ENSG00000087088", "ENSMUSG00000003873"]
    }

    TGF_SMADs = { # TGF-β/SMADs
        "TGFB1" : ["ENSG00000105329"],
        "TGFBR1" : ["ENSG00000106799"],
        "TGFBR2" : ["ENSG00000163513"],
        "SMAD2" : ["ENSG00000175387"],
        "SMAD3" : ["ENSG00000166949"],
        "SMAD4" : ["ENSG00000141646"],
        "SMAD7" : ["ENSG00000101665"],
        "KLF5" : ["ENSG00000102554"],
        "CYTL1" : ["ENSG00000170891"],
        "ANO1" : ["ENSG00000131620"]
    }

    Wnt_catenin = { # Wnt/β-catenin
        "WNT1" : ["ENSG00000125084"],
        "WNT2" : ["ENSG00000105989"],
        "WNT3" : ["ENSG00000108379"],
        "WNT4" : ["ENSG00000162552"],
        "WNT5A" : ["ENSG00000114251"],
        "IL1B" : ["ENSG00000125538"],
        "IL6" : ["ENSG00000136244"],
        "WIF1" : ["ENSG00000156076"],
        "WNT10B" : ["ENSG00000169884"],
        "SFRP1" : ["ENSG00000104332"],
        "MIR145" : ["ENSG00000276365"],
        "WNT11" : ["ENSG00000085741"],
        "LEF1" : ["ENSG00000138795"],
        "LRP5" : ["ENSG00000162337"],
        "LRP6" : ["ENSG00000070018"],
        "CTNNB1" : ["ENSG00000168036"],
        "APC" : ["ENSG00000134982"],
        "GSK3B" : ["ENSG00000082701"],
        "DKK2" : ["ENSG00000155011"]
    }

    #  !!!!! Sonic hedgehog not done !!!!!

    pathways = {"NLRP3/caspase1" : NLRP3_caspase1, 
                "TLR4/MyD88/NFkB" : TLR4_MyD88_NFkB,
                "Nrf2/HO-1" : Nrf2_HO_1,
                "Notch" : Notch,
                "Hippo" : Hippo,
                "RhoA/ROCK" : RhoA_ROCK,
                "MAPK" : MAPK,
                "PI3K/Akt" : PI3K_Akt,
                "JAK/STAT" : JAK_STAT,
                "TGF-β/SMADs" : TGF_SMADs,
                "Wnt/β-catenin" : Wnt_catenin
                }

    def generate_json_file(version):
        fp_pathways = get_filepath("data", f"pathways_v{str(version)}.json")

        with open(fp_pathways, "w") as f_pathways:
            json.dump(pathways, f_pathways)


    # generate_json_file(2)


    import gseapy as gp

    reactome = gp.get_library(name="Reactome_2022")

    fp = get_filepath("tmp", "pathways.txt")

    with open(fp, "w") as f:

        hearts= [k for k in reactome.keys() if "cardiac" in k.lower() or "heart" in k.lower()]

        print(hearts, file=f)

        for pathway, genes in reactome.items():
            print(pathway, len(genes), file=f)

    [
        "MyD88-independent TLR4 Cascade R-HSA-166166",
        "NOTCH1 Intracellular Domain Regulates Transcription R-HSA-2122947",
        "NOTCH2 Activation And Transmission Of Signal To Nucleus R-HSA-2979096",
        "NOTCH2 Intracellular Domain Regulates Transcription R-HSA-2197563",
        "NOTCH3 Activation And Transmission Of Signal To Nucleus R-HSA-9013507",
        "NOTCH3 Intracellular Domain Regulates Transcription R-HSA-9013508",
        "NOTCH4 Activation And Transmission Of Signal To Nucleus R-HSA-9013700",
        "NOTCH4 Intracellular Domain Regulates Transcription R-HSA-9013695",
        "Signaling By Hippo R-HSA-2028269",
        "RHO GTPases Activate ROCKs R-HSA-5627117",
        "MAPK Targets/ Nuclear Events Mediated By MAP Kinases R-HSA-450282",
        "MAPK Family Signaling Cascades R-HSA-5683057",
        "Negative Regulation Of PI3K/AKT Network R-HSA-199418",
        ""
    ]


    # dataset_fp = get_filepath("data", "myocardial_infarction_preprocessed_4.h5ad")
    # model_fp = get_filepath("model", "small-v2.ckpt")

    # model = load_model_with_cuda_if_avail(model_fp)
    # adata = scanpy.read_h5ad(dataset_fp)

    # for name, pathway in pathways.items():
    #     print(f"=========== {name} ===========")
    #     for gene, ls_ensembl in pathway.items():
    #         print(gene, ls_ensembl[0], ls_ensembl[0] in model.genes)