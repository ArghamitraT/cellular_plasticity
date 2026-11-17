# Chooses which VAE version steps 3 and 4 use. Change `version` and rerun steps 3-4.
#
#   "learned_variance"  vae_model_learned_variance.py
#                       reconstruction = squared error divided by a LEARNED spread per number
#                       (Gaussian with learned variance); score = log-likelihood estimated by
#                       importance sampling. Learns the shape of each cell type.
#   "mse"               vae_model_mse.py
#                       textbook VAE: reconstruction = plain MSE; score = -(MSE + KL).
#                       On SCimilarity embeddings it collapsed (KL = 0): it only learned the
#                       average cell of each type, so its score = distance to that average.
#
# Each version saves to its own folders (3_vaes_<run_name>/ and <query cell>/<run_name>/),
# so switching never overwrites the other version's results.

version = "learned_variance"

if version == "learned_variance":
    from vae_model_learned_variance import *  # noqa: F401,F403
elif version == "mse":
    from vae_model_mse import *  # noqa: F401,F403
else:
    raise ValueError(f"version must be 'learned_variance' or 'mse', not {version!r}")
