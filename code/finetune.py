import os
from typing import Any, Dict, List, Optional

import bionty as bt
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import scanpy as sc
import seaborn as sns
import torch
import torch.nn.functional as F
from torch.masked import masked_tensor
from anndata import AnnData, concat
from scdataloader import Collator, Preprocessor
from scdataloader.data import SimpleAnnDataset
from scdataloader.utils import get_descendants, random_str
from scib_metrics.benchmark import Benchmarker
from scipy.stats import spearmanr
from simpler_flash import FlashTransformer
from sklearn.metrics import f1_score
from torch.utils.data import DataLoader
from tqdm import tqdm

from scprint2.model import loss
from scprint2.model import utils

FILE_LOC = os.path.dirname(os.path.realpath(__file__))


class FinetuneGRN:
    def __init__(
        self,
        batch_key: str = "batch",
        max_len: int = 5000,
        learn_batches_on: Optional[str] = None,
        num_workers: int = 8,
        batch_size: int = 16,
        num_epochs: int = 8,
        do_mmd_on: Optional[str] = None,
        lr: float = 0.0002,
        ft_mode: str = "xpressor",
        frac_train: float = 0.8,
        loss_scalers: dict = {},
        use_knn: bool = True,
        mask_frac: float = 0.3
    ) -> torch.nn.Module:
        """
        Embedder a class to embed and annotate cells using a model

        Args:
            batch_key (str, optional): The key in adata.obs that indicates the batch information. Defaults to "batch".
            learn_batches_on (str, optional): The key in adata.obs to learn batch embeddings on. Defaults to None.
                if none, will not learn the batch embeddings.
                the goal is e.g. when having a new species, to learn an embedding for it during finetuning and replace
                the "learn_batches_on" embedding in the model with it, in this case it should be "organism_ontology_term_id".
                batch correction might indeed be better learnt with this additional argument in some cases.
            do_mmd_on (str, optional):The key in adata.obs to learn batch embeddings on. Defaults to None.
                this embedding should have less batch information in it, after finetuning.
            batch_size (int, optional): The size of the batches to be used in the DataLoader. Defaults to 64.
            num_workers (int, optional): The number of worker processes to use for data loading. Defaults to 8.
            max_len (int, optional): The maximum length of the sequences to be processed. Defaults to 5000.
            lr (float, optional): The learning rate for the optimizer. Defaults to 0.0002.
            num_epochs (int, optional): The number of epochs to train the model. Defaults to 8.
            ft_mode (str, optional): The fine-tuning mode, either "xpressor" or "full". Defaults to "xpressor".
            frac_train (float, optional): The fraction of data to be used for training. Defaults to 0.8.
            loss_scalers (dict, optional): A dictionary specifying the scaling factors for different loss components. Defaults to {}.
                expr, class, mmd, kl
            use_knn (bool, optional): Whether to use k-nearest neighbors information. Defaults to True.
            mask_frac (float, optional): Fraction of the input to mask. Defaults to 0.3.
        """
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.batch_key = batch_key
        self.learn_batches_on = learn_batches_on
        self.max_len = max_len
        self.lr = lr
        self.num_epochs = num_epochs
        self.ft_mode = ft_mode
        self.frac_train = frac_train
        self.batch_emb = None
        self.batch_encoder = {}
        self.do_mmd_on = do_mmd_on
        self.loss_scalers = loss_scalers
        self.use_knn = use_knn
        self.mask_ratio = mask_frac

    def __call__(
        self,
        model: torch.nn.Module,
        adata: AnnData = None,
        train_data: AnnData = None,
        val_data: AnnData = None,
    ) -> torch.nn.Module:
        """
        __call__ function to call the embedding

        Args:
            model (torch.nn.Module): The scPRINT model to be used for embedding and annotation.
            adata (AnnData): The annotated data matrix of shape n_obs x n_vars. Rows correspond to cells and columns to genes.
                Defaults to None.
                if provided, it will be split into training and validation sets.
            train_data (AnnData, optional): The training data. Defaults to None.
                if adata is provided, this will be ignored.
            val_data (AnnData, optional): The validation data. Defaults to None.
                if adata is provided, this will be ignored.

        Raises:
            ValueError: If the model does not have a logger attribute.
            ValueError: If the model does not have a global_step attribute.

        Returns:
            torch.nn.Module: the fine-tuned model
        """
        # one of "all" "sample" "none"
        model.predict_mode = "none"
        if self.ft_mode == "xpressor":
            for val in model.parameters():
                val.requires_grad = False
                # setting all to TRUE

            for val in model.cell_transformer.parameters():
                val.requires_grad = True
            for val in model.transformer.blocks[-1].parameters():
                val.requires_grad = True
            for i in model.transformer.blocks:
                i.cross_attn.requires_grad = True
            for val in model.compressor.parameters():
                val.requires_grad = True
        elif self.ft_mode == "full":
            for val in model.parameters():
                val.requires_grad = True
        else:
            raise ValueError("ft_mode must be one of 'xpressor' or 'full'")

        # PREPARING THE DATA
        if adata is not None:
            n_train = int(self.frac_train * len(adata))
            train_idx = np.random.choice(len(adata), n_train, replace=False)
            val_idx = np.setdiff1d(np.arange(len(adata)), train_idx)

            train_data = adata[train_idx].copy()
            val_data = adata[val_idx].copy()

            print(f"Training data: {train_data.shape}")
            print(f"Validation data: {val_data.shape}")

        mencoders = {}
        for k, v in model.label_decoders.items():
            mencoders[k] = {va: ke for ke, va in v.items()}
        # this needs to remain its original name as it is expect like that by collator, otherwise need to send org_to_id as params

        # create datasets
        self.batch_encoder = {
            i: n
            for n, i in enumerate(
                train_data.obs[self.batch_key].astype("category").cat.categories
            )
        }
        mencoders[self.batch_key] = self.batch_encoder
        train_dataset = SimpleAnnDataset(
            train_data,
            obs_to_output=[self.batch_key],
            get_knn_cells=model.expr_emb_style == "metacell" and self.use_knn,
            encoder=mencoders,
        )
        if val_data is not None:
            self.batch_encoder.update(
                {
                    i: n + len(self.batch_encoder)
                    for n, i in enumerate(
                        val_data.obs[self.batch_key].astype("category").cat.categories
                    )
                    if i not in self.batch_encoder
                }
            )
            mencoders[self.batch_key] = self.batch_encoder
            val_dataset = SimpleAnnDataset(
                val_data,
                obs_to_output=[self.batch_key],
                get_knn_cells=model.expr_emb_style == "metacell" and self.use_knn,
                encoder=mencoders,
            )

        # Create collator
        collator = Collator(
            organisms=model.organisms,
            valid_genes=model.genes,
            class_names=[self.batch_key],
            how="random expr",  # or "all expr" for full expression
            max_len=self.max_len,
            org_to_id=mencoders.get("organism_ontology_term_id", {}),
        )

        # Create data loaders
        train_loader = DataLoader(
            train_dataset,
            collate_fn=collator,
            batch_size=self.batch_size,  # Adjust based on GPU memory
            num_workers=self.num_workers,
            shuffle=True,
        )
        if val_data is not None:
            val_loader = DataLoader(
                val_dataset,
                collate_fn=collator,
                batch_size=self.batch_size,
                num_workers=self.num_workers,
                shuffle=False,
            )

        if self.learn_batches_on is not None:
            if val_data is not None:
                print(
                    "all batch key values in val_data should also be present in train_adata!!!"
                )
            self.batch_emb = torch.nn.Embedding(
                num_embeddings=train_data.obs[self.batch_key].nunique(),
                embedding_dim=(
                    model.compressor[self.learn_batches_on].fc_mu.weight.shape[0]
                    if hasattr(model, "compressor")
                    else model.d_model
                ),
            )

        ## PREPARING THE OPTIM
        all_params = (
            list(model.parameters())
            + (
                list(self.batch_emb.parameters())
                if self.learn_batches_on is not None
                else []
            )
        )

        # Setup optimizer
        optimizer = torch.optim.AdamW(
            all_params,
            lr=self.lr,
            weight_decay=0.01,
            betas=(0.9, 0.999),
            eps=1e-8,
        )

        # Setup scheduler
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode="min", factor=0.5, patience=2
        )

        # Setup automatic mixed precision
        scaler = torch.cuda.amp.GradScaler() if torch.cuda.is_available() else None

        for k, i in model.mat_labels_hierarchy.items():
            model.mat_labels_hierarchy[k] = i.to(model.device)

        ## train
        for epoch in range(self.num_epochs):
            print(f"\nEpoch {epoch + 1}/{self.num_epochs}")
            print(f"Current learning rate: {optimizer.param_groups[0]['lr']:.2e}")

            # Training phase
            train_loss = 0.0
            train_steps = 0
            avg_expr = 0
            avg_cls = 0
            avg_mmd = 0

            pbar = tqdm(train_loader, desc="Training")
            model.train()
            for batch_idx, batch in enumerate(pbar):
                optimizer.zero_grad()
                total_loss, loss_expr = self.batch_corr_pass(
                    batch,
                    model,
                    self.mask_ratio
                )
                # Backward pass
                scaler.scale(total_loss).backward()
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                scaler.step(optimizer)
                scaler.update()

                train_loss += total_loss.item()
                train_steps += 1
                avg_cls += cls_loss.item()
                avg_expr += loss_expr.item()
                avg_mmd += mmd
                # Update progress bar
                pbar.set_postfix(
                    {
                        "loss": f"{total_loss.item():.4f}",
                        "avg_loss": f"{train_loss / train_steps:.4f}",
                        "cls_loss": f"{cls_loss.item():.4f}",
                        "mmd_loss": f"{mmd:.4f}",
                        "expr_loss": f"{loss_expr.item():.4f}",
                    }
                )

            # Validation phase
            if val_data is not None:
                model.eval()
                val_loss = 0.0
                val_steps = 0
                val_loss_expr = 0.0
                val_mmd = 0.0
                val_cls = 0.0
                val_loss_to_prt = 0.0

                with torch.no_grad():
                    for batch in val_loader:  # tqdm(val_loader, desc="Validation"):
                        loss_val, loss_expr = self.batch_corr_pass(
                            batch,
                            model,
                            self.mask_ratio
                        )
                        val_loss_to_prt += loss_val.item()
                        val_loss += loss_val.item()
                        val_steps += 1
                        val_loss_expr += loss_expr.item()
                        val_mmd += mmd
                        val_cls += cls_loss.item()
                try:
                    avg_val_loss = val_loss_to_prt / val_steps
                    avg_train_loss = train_loss / train_steps
                except ZeroDivisionError:
                    print(
                        "Error: Division by zero occurred while calculating average losses."
                    )
                    avg_train_loss = 0
                print(
                    "cls_loss: {:.4f}, mmd_loss: {:.4f}, expr_loss: {:.4f}".format(
                        val_cls / val_steps,
                        val_mmd / val_steps,
                        val_loss_expr / val_steps,
                    )
                )
                print(f"Train Loss: {avg_train_loss:.4f}, Val Loss: {avg_val_loss:.4f}")

                # Store LR before scheduler step for comparison
                lr_before = optimizer.param_groups[0]["lr"]

                # Update learning rate
                scheduler.step(avg_val_loss)

                # Check if LR was reduced
                lr_after = optimizer.param_groups[0]["lr"]
                if lr_after < lr_before:
                    print(
                        f"🔻 Learning rate reduced from {lr_before:.2e} to {lr_after:.2e} (factor: {lr_after / lr_before:.3f})"
                    )
                else:
                    print(f"✅ Learning rate unchanged: {lr_after:.2e}")

                # Early stopping check (simple implementation)
                if epoch > 3 and val_loss / val_steps > 1.3 * avg_train_loss:
                    print("Early stopping due to overfitting")
                    break

        print("Manual fine-tuning completed!")
        model.eval()
        return model

    def batch_corr_pass(self, batch, model, mask_frac):
        genes = batch["genes"].to(model.device)
        expression = batch["x"].to(model.device)
        depth = batch["depth"].to(model.device)
        class_elem = batch["class"].long().to(model.device)
        total_loss = 0
        mask_ratio = mask_frac
        # Forward pass with automatic mixed precisio^n
        with torch.cuda.amp.autocast():
            # Forward pass
            output = model.forward(
                genes,
                expression,
                req_depth=depth,
                depth_mult=expression.sum(1),
                do_class=True,
                metacell_token=torch.zeros_like(depth),
                mask=simple_masker([batch, seq_length], mask_ratio=mask_ratio)
            )

            if self.learn_batches_on is not None:
                batch_pos = model.classes.index(self.learn_batches_on) + 1
                output["output_cell_embs"][:, batch_pos, :] = self.batch_emb(
                    class_elem[:, -1]
                )

            output_gen = model._generate(
                cell_embs=output["output_cell_embs"],
                gene_pos=gene_pos,
                depth_mult=expression.sum(1),
                req_depth=depth,
            )
            ## generate expr loss
            if "zero_logits" in output_gen:
                loss_expr = loss.zinb(
                    theta=output_gen["disp"],
                    pi=output_gen["zero_logits"],
                    mu=output_gen["mean"],
                    target=expression,
                )
                if model.zinb_and_mse:
                    loss_expr += (
                        loss.mse(
                            input=torch.log(output_gen["mean"] + 1)
                            * (1 - torch.sigmoid(output_gen["zero_logits"])),
                            target=torch.log(expression + 1),
                        )
                        / 10  # scale to make it more similar to the zinb
                    )
            else:
                loss_expr = loss.mse(
                    input=torch.log(output_gen["mean"] + 1),
                    target=torch.log(expression + 1),
                )
            # Add expression loss to total
            total_loss += loss_expr * self.loss_scalers.get("expr", 0.5)

        return total_loss, loss_expr
