# VAE version "learned_variance" (chosen in vae_model.py; used by steps 3 and 4).
#
# Input: a SCimilarity embedding (128 numbers).
# Encoder squeezes it to latent_dim numbers (h); decoder rebuilds the 128 numbers from h.
# Trained on reference cells of ONE cell type, it learns what that type normally looks like.
#
# Difference from the textbook version (vae_model_mse.py): the reconstruction loss divides
# each squared error by a spread the VAE LEARNS for each of the 128 numbers (dec_logvar):
#   recon = 1/2 * sum over 128 numbers of [ (x - rebuilt)^2 / spread^2 + log spread^2 + log 2pi ]
#         = Gaussian negative log-likelihood with learned variance (a weighted MSE).
#   loss  = recon + KL to N(0, I)   (minimizing = maximizing the ELBO)
# score(x) = log p(x), estimated by importance sampling over n_importance summaries h
#            (higher = cell looks more like the training type).

import numpy as np
import torch
from torch import nn

latent_dim = 128     # size of the summary h (mentor: 128, same as the input; first runs used 10)
hidden = 256
train_steps = 20000  # 4000 was too few: validation curves were still rising in a straight line
batch_size = 256
check_every = 200    # steps between checks on the validation cells
n_importance = 1000  # samples per cell for the log-likelihood estimate (200 was noisier, and
                     # too low for the query cell by 20-130; papers often use 5000)
run_name = f"learnedvar_latent{latent_dim}"  # names the results folders for this version
version_description = ("VAE with learned spread per number: loss = Gaussian reconstruction with "
                       "learned variance (weighted MSE) + KL to N(0, I); score = log-likelihood by "
                       f"importance sampling ({n_importance} samples per cell)")
score_label = "log-likelihood under that type's VAE"
recon_label = "reconstruction loss\n(Gaussian, learned spread)"


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
        self.dec_logvar = nn.Parameter(torch.full((n_in,), np.log(0.03 ** 2)))  # learned spread per number

    def encode(self, x):
        mu, logvar = self.encoder(x).chunk(2, dim=-1)
        return mu, logvar

    def log_p_x_given_h(self, x, h):
        # Gaussian log-likelihood of the cell given summary h, with the learned spread
        mean = self.decoder(h)
        return -0.5 * (((x - mean) ** 2) / self.dec_logvar.exp()
                       + self.dec_logvar + np.log(2 * np.pi)).sum(-1)

    def loss_per_cell(self, x, sample=True):
        # sample=True during training (reparameterization trick); sample=False uses h = mu
        mu, logvar = self.encode(x)
        h = mu + torch.randn_like(mu) * (0.5 * logvar).exp() if sample else mu
        recon = -self.log_p_x_given_h(x, h)
        kl = 0.5 * (mu ** 2 + logvar.exp() - 1 - logvar).sum(-1)
        return recon, kl

    def elbo(self, x, sample=True):
        # ELBO = log p(x|h) - KL = -(recon + KL); training maximizes this (loss = -ELBO)
        recon, kl = self.loss_per_cell(x, sample)
        return -(recon + kl)

    @torch.no_grad()
    def score(self, x):
        # log p(x) by importance sampling: average p(x|h) p(h) / q(h|x) over n_importance h's
        x = torch.as_tensor(x, dtype=torch.float32)
        mu, logvar = self.encode(x)
        std = (0.5 * logvar).exp()
        out = []
        for _ in range(n_importance):
            h = mu + torch.randn_like(mu) * std
            log_q = (-0.5 * (((h - mu) / std) ** 2 + logvar + np.log(2 * np.pi))).sum(-1)
            log_prior = (-0.5 * (h ** 2 + np.log(2 * np.pi))).sum(-1)
            out.append(self.log_p_x_given_h(x, h) + log_prior - log_q)
        return (torch.logsumexp(torch.stack(out), dim=0) - np.log(n_importance)).numpy()


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
        loss = -vae.elbo(batch).mean()  # = mean of (recon + KL)
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
