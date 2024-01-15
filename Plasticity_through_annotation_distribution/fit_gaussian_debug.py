import os
import csv
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import multivariate_normal
import utils_AT

# Define the path to the embeddings folder and the output CSV file
embeddings_folder = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/cell_type_embeddings')
Gaussian_folder = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')
output_csv = os.path.join(os.getcwd(), '../../files/scimilarity_trainingset_cellnum_AuthorLable.csv')

# Prepare to write to CSV file

# Loop over each file in the embeddings directory
for file_name in os.listdir(embeddings_folder):
    if file_name.endswith("_embeddings.npy"):
        cell_type = file_name.replace('_embeddings.npy', '')
        file_path = os.path.join(embeddings_folder, file_name)

        # Load data for the cell type
        data = np.load(file_path)

        # Perform PCA on the data to reduce it to 3 dimensions
        pca = PCA(n_components=3)
        reduced_data = pca.fit_transform(data)

        # Fit a Gaussian Mixture Model to the reduced data
        gmm = GaussianMixture(n_components=1, covariance_type='full')
        gmm.fit(reduced_data)

        # The mean and covariance of the fitted Gaussian distribution
        mean = gmm.means_[0]
        covariance = gmm.covariances_[0]

        # Create a grid for the 3D plot
        grid_x, grid_y = np.meshgrid(
            np.linspace(reduced_data[:, 0].min(), reduced_data[:, 0].max(), num=100),
            np.linspace(reduced_data[:, 1].min(), reduced_data[:, 1].max(), num=100)
        )
        grid_z = np.zeros_like(grid_x)

        # Flatten the grid for PDF evaluation
        grid_flat = np.stack((grid_x.ravel(), grid_y.ravel(), grid_z.ravel()), axis=-1)

        # Evaluate the PDF of the Gaussian distribution on the grid
        rv = multivariate_normal(mean, covariance)
        pdf = rv.pdf(grid_flat).reshape(grid_x.shape)

        # Plot the surface plot
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(grid_x, grid_y, pdf, cmap='viridis')

        ax.set_xlabel('PCA Dimension 1')
        ax.set_ylabel('PCA Dimension 2')
        ax.set_zlabel('Probability Density')

        plt.title(f'3D Surface Plot of {cell_type}')
        plt.savefig(os.path.join('../figures/fitted_gaussian_author_lables_3d', utils_AT.create_image_name(cell_type + "_3d_gaussian_")))

        #plt.show()

        # # Evaluate the PDF on the reduced data
        # rv = multivariate_normal(mean, covariance)
        # pdf = rv.pdf(reduced_data)
        #
        # # Normalize the PDF for better visualization
        # pdf_normalized = (pdf - pdf.min()) / (pdf.max() - pdf.min())
        #
        # # Plot the 3D scatter plot
        # fig = plt.figure(figsize=(10, 7))
        # ax = fig.add_subplot(111, projection='3d')
        #
        # # Scatter plot with the normalized PDF values used for the color map
        # scatter = ax.scatter(reduced_data[:, 0], reduced_data[:, 1], reduced_data[:, 2], c=pdf_normalized,
        #                      cmap='viridis')
        #
        # # Create a color bar
        # cbar = fig.colorbar(scatter, shrink=0.5, aspect=5)
        # cbar.set_label('Probability Density')
        #
        # # Set labels
        # ax.set_xlabel('PCA Dimension 1')
        # ax.set_ylabel('PCA Dimension 2')
        # ax.set_zlabel('PCA Dimension 3')
        #
        # plt.title('3D Scatter Plot of Gaussian Mixture Model')
        # plt.show()

        # Save the mean and covariance
        # np.save(os.path.join(Gaussian_folder, f'{cell_type}_mean.npy'), gmm.means_)
        # np.save(os.path.join(Gaussian_folder, f'{cell_type}_cov.npy'), gmm.covariances_)

        # # Optionally perform PCA and plot
        # pca = PCA(n_components=2)
        # reduced_data = pca.fit_transform(data)
        #
        # plt.figure(figsize=(10, 8))
        # sns.scatterplot(x=reduced_data[:, 0], y=reduced_data[:, 1], color='blue', alpha=0.5)
        # # Create grid and multivariate normal
        # x = np.linspace(min(reduced_data[:, 0]), max(reduced_data[:, 0]), num=100)
        # y = np.linspace(min(reduced_data[:, 1]), max(reduced_data[:, 1]), num=100)
        # X, Y = np.meshgrid(x, y)
        #
        # # Fit a multivariate normal distribution to the PCA-reduced data
        # # Calculate the mean of the reduced data for the multivariate normal
        # mean_reduced = np.mean(reduced_data, axis=0)
        # # Use the PCA object's covariance matrix, which is already in the reduced form
        # cov_reduced = np.cov(reduced_data, rowvar=False)
        #
        # # Create the multivariate normal model
        # model = multivariate_normal(mean=mean_reduced, cov=cov_reduced)
        #
        # # Evaluate the PDF of the model on the grid
        # pos = np.empty(X.shape + (2,))
        # pos[:, :, 0] = X
        # pos[:, :, 1] = Y
        # Z = model.pdf(pos)
        #
        # # Plot the PDF contour
        # plt.contour(X, Y, Z, levels=10, cmap="Reds")
        #
        # plt.title(f'PCA of {cell_type} embeddings')
        #
        # # Optional: Save the plot
        # plt.savefig(os.path.join('../figures/', utils_AT.create_image_name(cell_type + "_2d_gaussian_")))
        #
        # plt.close()  # Close the plot to free memory

print("Gaussian fitting completed for all cell types and shapes saved to CSV.")




# import os
# from sklearn.mixture import GaussianMixture
# from sklearn.decomposition import PCA
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import multivariate_normal
# import utils_AT
#
# # Define the path to the embeddings folder
# embeddings_folder = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/cell_type_embeddings')
# Gaussian_folder = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')
#
# # Loop over each file in the embeddings directory
# for file_name in os.listdir(embeddings_folder):
#     if file_name.endswith("_embeddings.npy"):
#         cell_type = file_name.replace('_embeddings.npy', '')
#         file_path = os.path.join(embeddings_folder, file_name)
#
#         # Load data for the cell type
#         data = np.load(file_path)
#         print(f'{cell_type}: {data.shape[0]}')
#
#         # Fit Gaussian Mixture Model
#         gmm = GaussianMixture(n_components=1, covariance_type='full')
#         gmm.fit(data)
#
#         # Save the mean and covariance
#         np.save(os.path.join(Gaussian_folder, f'{cell_type}_mean.npy'), gmm.means_)
#         np.save(os.path.join(Gaussian_folder, f'{cell_type}_cov.npy'), gmm.covariances_)
#
#         # Optionally perform PCA and plot
#         pca = PCA(n_components=2)
#         reduced_data = pca.fit_transform(data)
#
#         plt.figure(figsize=(10, 8))
#         sns.scatterplot(x=reduced_data[:, 0], y=reduced_data[:, 1], color='blue', alpha=0.5)
#         # Create grid and multivariate normal
#         x = np.linspace(min(reduced_data[:, 0]), max(reduced_data[:, 0]), num=100)
#         y = np.linspace(min(reduced_data[:, 1]), max(reduced_data[:, 1]), num=100)
#         X, Y = np.meshgrid(x, y)
#
#         # Fit a multivariate normal distribution to the PCA-reduced data
#         # Calculate the mean of the reduced data for the multivariate normal
#         mean_reduced = np.mean(reduced_data, axis=0)
#         # Use the PCA object's covariance matrix, which is already in the reduced form
#         cov_reduced = np.cov(reduced_data, rowvar=False)
#
#         # Create the multivariate normal model
#         model = multivariate_normal(mean=mean_reduced, cov=cov_reduced)
#
#         # Evaluate the PDF of the model on the grid
#         pos = np.empty(X.shape + (2,))
#         pos[:, :, 0] = X
#         pos[:, :, 1] = Y
#         Z = model.pdf(pos)
#
#         # Plot the PDF contour
#         plt.contour(X, Y, Z, levels=10, cmap="Reds")
#
#         plt.title(f'PCA of {cell_type} embeddings')
#
#         # Optional: Save the plot
#         plt.savefig(os.path.join('../figures/', utils_AT.create_image_name(cell_type+"_2d_gaussian_")))
#
#         plt.close()  # Close the plot to free memory
#
# print("Gaussian fitting completed for all cell types.")
#


# """ This function fits Gaussian to the existing cell types distribution """
#
# import os
#
#
# from sklearn.mixture import GaussianMixture
# import numpy as np
# from sklearn.decomposition import PCA
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import multivariate_normal
# import utils_AT as util
#
#
# #data = np.random.randn(10000, 128)  # Replace this with your actual data
# embedding_path = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/cell_type_embeddings/capillary_endothelial_cell_embeddings.npy')
# data = np.load(embedding_path)
#
# # Calculate the empirical mean and covariance matrix
# mean_vector = np.mean(data, axis=0)
# covariance_matrix = np.cov(data, rowvar=False)
#
# # Assuming `data` is your dataset as described above
# gmm = GaussianMixture(n_components=1, covariance_type='full')
#
# # Fit the model
# gmm.fit(data)
#
# # The mean and covariance of the fitted Gaussian distribution
# mean_vector = gmm.means_
# covariance_matrix = gmm.covariances_
#
#
# # ... assume 'data' is your original data and 'mean_vector' is the original mean before PCA ...
#
# # Perform PCA on the data
# pca = PCA(n_components=2)
# reduced_data = pca.fit_transform(data)
#
# # Plot the reduced data
# plt.figure(figsize=(10, 8))
# sns.scatterplot(x=reduced_data[:, 0], y=reduced_data[:, 1], color='blue', alpha=0.5)
#
# # Create grid and multivariate normal
# x = np.linspace(min(reduced_data[:, 0]), max(reduced_data[:, 0]), num=100)
# y = np.linspace(min(reduced_data[:, 1]), max(reduced_data[:, 1]), num=100)
# X, Y = np.meshgrid(x, y)
#
# # Fit a multivariate normal distribution to the PCA-reduced data
# # Calculate the mean of the reduced data for the multivariate normal
# mean_reduced = np.mean(reduced_data, axis=0)
# # Use the PCA object's covariance matrix, which is already in the reduced form
# cov_reduced = np.cov(reduced_data, rowvar=False)
#
# # Create the multivariate normal model
# model = multivariate_normal(mean=mean_reduced, cov=cov_reduced)
#
# # Evaluate the PDF of the model on the grid
# pos = np.empty(X.shape + (2,))
# pos[:, :, 0] = X
# pos[:, :, 1] = Y
# Z = model.pdf(pos)
#
# # Plot the PDF contour
# plt.contour(X, Y, Z, levels=10, cmap="Reds")
# plt.savefig()
# plt.show('../figures/' + util.create_image_name("2D_gaussian_"))
#
#
