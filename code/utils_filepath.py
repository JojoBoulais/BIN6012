from pathlib import Path
import os

__all__ = ["get_data_filepath", "get_model_filepath"]

CURRENT_DIR = Path(__file__).parent

DATA_DIR = os.path.join(CURRENT_DIR, "..", "data")
MODEL_DIR = os.path.join(CURRENT_DIR, "..", "model_ckpt")

def get_data_filepath(data_file):
    return os.path.join(DATA_DIR, data_file)

def get_model_filepath(model_file):
    return os.path.join(MODEL_DIR, model_file)
