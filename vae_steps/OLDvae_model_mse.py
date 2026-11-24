# VAE version "mse" (chosen in vae_model.py; used by steps 3 and 4).
#
# Input: a SCimilarity embedding (128 numbers).
# Encoder squeezes it to latent_dim numbers (h); decoder rebuilds the 128 numbers from h.
# Trained on reference cells of ONE cell type, it learns what that type normally looks like.
#
# Standard VAE:
#   loss = reconstruction loss (summed squared error, i.e. MSE with reduction="sum")
#        + KL divergence between q(h|x) = N(mu, sigma^2) and the prior p(h) = N(0, I)
#   Minimizing this loss = maximizing the ELBO.
# score(x) = elbo_score(x) = -loss for each cell (higher = cell looks more like the training type).
#
# Result on SCimilarity embeddings (2026-11-17): collapsed, KL = 0 for every type; the VAE
# rebuilds every cell as the average cell of its type, so score = -(squared distance to
# that average).

import numpy as np
import torch
from torch import nn
import torch.nn.functional as F

latent_dim = 128     # size of the summary h (mentor: 128, same as the input; was 10).
                     # VAEs and results are saved per version (3_vaes_<run_name>/, <run_name>/)
hidden = 256
train_steps = 20000  # 4000 was too few: validation curves were still rising in a straight line
batch_size = 256
check_every = 200    # steps between checks on the validation cells
run_name = f"mse_latent{latent_dim}"  # names the results folders for this version
version_description = ("textbook VAE: loss = reconstruction (MSE, summed over 128 numbers) + KL to N(0, I); "
                       "score = -(reconstruction + KL), with h = mu")
score_label = "VAE score = -(reconstruction + KL)"
recon_label = "reconstruction loss (MSE)"


def write_run_info(out_dir, query_cell, script):
    # Adds a line to <out_dir>/run_info.txt saying which VAE version and settings made the results.
    import datetime
    import os
    with open(os.path.join(out_dir, "run_info.txt"), "a") as f:
        f.write(f"{datetime.datetime.now():%Y-%m-%d %H:%M}  {script}  query cell {query_cell}\n"
                f"    VAE version: {run_name} ({os.path.abspath(__file__)})\n"
                f"    {version_description}\n"
                f"    latent_dim={latent_dim}, hidden={hidden}, train_steps={train_steps}, "
                f"batch_size={batch_size}, check_every={check_every}, learning rate=1e-3\n")


class VAE(nn.Module):
    def __init__(self, n_in=128):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(n_in, hidden), nn.ReLU(),
                                     nn.Linear(hidden, hidden), nn.ReLU(),
                                     nn.Linear(hidden, 2 * latent_dim))
        self.decoder = nn.Sequential(nn.Linear(latent_dim, hidden), nn.ReLU(),
                                     nn.Linear(hidden, hidden), nn.ReLU(),
                                     nn.Linear(hidden, n_in))

    def encode(self, x):
        mu, logvar = self.encoder(x).chunk(2, dim=-1)
        return mu, logvar

    def loss_per_cell(self, x, sample=True):
        # sample=True during training (reparameterization trick);
        # sample=False when scoring, so the score is the same every time (uses h = mu)
        mu, logvar = self.encode(x)
        h = mu + torch.randn_like(mu) * (0.5 * logvar).exp() if sample else mu
        recon = F.mse_loss(self.decoder(h), x, reduction="none").sum(-1)
        kl = 0.5 * (mu ** 2 + logvar.exp() - 1 - logvar).sum(-1)
        return recon, kl

    @torch.no_grad()
    def elbo_score(self, x):
        # -(reconstruction + KL) per cell; higher = more typical of the training type
        recon, kl = self.loss_per_cell(torch.as_tensor(x, dtype=torch.float32), sample=False)
        return (-(recon + kl)).numpy()

    def score(self, x):
        # the score steps 3-4 use (same name in every VAE version)
        return self.elbo_score(x)


def train_vae(x_train, rng):
    # 10% of the training cells are validation cells: checked every check_every steps,
    # and the model with the lowest validation loss is kept (stops it memorising the training cells).
    # Returns the VAE and its validation curve (step, mean recon, mean KL).
    perm = rng.permutation(len(x_train))
    n_val = len(x_train) // 10
    x_val = torch.tensor(x_train[perm[:n_val]], dtype=torch.float32)
    x_fit = torch.tensor(x_train[perm[n_val:]], dtype=torch.float32)

    vae = VAE()
    opt = torch.optim.Adam(vae.parameters(), lr=1e-3)
    best_val, best_state, curve = np.inf, None, []
    for step in range(train_steps):
        batch = x_fit[torch.randint(len(x_fit), (batch_size,))]
        recon, kl = vae.loss_per_cell(batch)
        loss = (recon + kl).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
        if step % check_every == 0 or step == train_steps - 1:
            with torch.no_grad():
                r, k = vae.loss_per_cell(x_val, sample=False)
                val = (r + k).mean().item()
            curve.append((step, r.mean().item(), k.mean().item()))
            if val < best_val:
                best_val, best_state = val, {n: p.clone() for n, p in vae.state_dict().items()}
    vae.load_state_dict(best_state)
    vae.eval()
    return vae, curve
