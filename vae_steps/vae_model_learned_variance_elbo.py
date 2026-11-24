# VAE version "learned_variance_elbo" (chosen in vae_model.py; used by steps 3 and 4).
#
# SAME model and training as vae_model_learned_variance.py; ONLY the score is different:
#   learned_variance       score = log p(x) by importance sampling (logsumexp over 1000 h's)
#   learned_variance_elbo  score = -(reconstruction + KL) with h = mu, an approximation of the
#                          ELBO (same score code as fixed_variance)
# The ELBO is a lower bound: log p(x) = ELBO + gap, and the gap can be bigger for unusual
# cells. Comparing the two versions shows how much the choice of score matters.
# It reuses the trained learned_variance VAEs (vae_folder below): only run step 4 for it.
#
# Input: a SCimilarity embedding (128 numbers).
# Encoder squeezes it to latent_dim numbers (h); decoder rebuilds the 128 numbers from h.
# Trained on reference cells of ONE cell type, it learns what that type normally looks like.
#
# Difference from the textbook version (vae_model_fixed_variance.py): the reconstruction loss divides
# each squared error by a spread the VAE LEARNS for each of the 128 numbers (dec_logvar):
#   recon = 1/2 * sum over 128 numbers of [ (x - rebuilt)^2 / spread^2 + log spread^2 + log 2pi ]
#         = Gaussian negative log-likelihood with learned variance (a weighted MSE).
#   loss  = recon + KL to N(0, I)   (minimizing = maximizing the ELBO)
# score(x) = elbo_score(x) = -loss for each cell, with h = mu
#            (higher = cell looks more like the training type).

import numpy as np
import torch
from torch import nn

latent_dim = 10      # size of the summary h (10; also run with 128 -> results in *_learnedvar_latent128)
hidden = 256
train_steps = 20000  # 4000 was too few: validation curves were still rising in a straight line
batch_size = 256
check_every = 200    # steps between checks on the validation cells
run_name = f"learnedvar_elbo_mu_latent{latent_dim}"  # names the results folders for this version
vae_folder = f"3_vaes_learnedvar_latent{latent_dim}"  # reuses the learned_variance VAEs
                     # (same model and training); only the score is new
version_description = ("VAE with learned spread per number: loss = Gaussian reconstruction with "
                       "learned variance (weighted MSE) + KL to N(0, I); "
                       "score = -(reconstruction + KL), with h = mu")
score_label = "ELBO (h = mu) under that type's VAE"
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
    # The training loss of every step is kept in vae.train_log (step, mean recon, mean KL).
    perm = rng.permutation(len(x_train))
    n_val = len(x_train) // 10
    x_val = torch.tensor(x_train[perm[:n_val]], dtype=torch.float32)
    x_fit = torch.tensor(x_train[perm[n_val:]], dtype=torch.float32)

    vae = VAE()
    opt = torch.optim.Adam(vae.parameters(), lr=1e-3)
    best_val, best_state, curve, train_log = np.inf, None, [], []
    for step in range(train_steps):
        batch = x_fit[torch.randint(len(x_fit), (batch_size,))]
        recon, kl = vae.loss_per_cell(batch)
        loss = (recon + kl).mean()  # = -ELBO, same number as -vae.elbo(batch).mean()
        train_log.append((step, recon.mean().item(), kl.mean().item()))
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
    vae.train_log = train_log
    return vae, curve
