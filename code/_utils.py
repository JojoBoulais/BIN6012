from scprint2 import scPRINT2
import torch

def load_model_with_cuda_if_avail(model_fp):
    """
    Load model from filepath, enable cuda if available
    """
    model = scPRINT2.load_from_checkpoint(
        model_fp,
        precpt_gene_emb=None,
        gene_pos_file=None,
    )
    if not torch.cuda.is_available():
        model = model.to(torch.float32)

    model = model.to("cuda" if torch.cuda.is_available() else "cpu")

    return model