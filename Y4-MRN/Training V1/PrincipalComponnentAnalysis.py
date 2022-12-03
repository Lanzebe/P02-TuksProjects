from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt



InputData = np.load('TrainingSetInputs.npy', allow_pickle=True)
OutputData = np.load('TrainingSetOutputs.npy', allow_pickle=True)
AdditionalData = np.load('TrainingSetAdditional.npy', allow_pickle=True)


print('Input shape',np.shape(InputData))
print('Output shape',np.shape(OutputData))
print('Output shape',np.shape(AdditionalData))

x = np.array([InputData[:,0], InputData[:,1],InputData[:,2],AdditionalData[:,0],AdditionalData[:,1],AdditionalData[:,2],OutputData[:,0]]).T


pc = PCA(n_components=6)
pc = pc.fit(x)

print("components:", pc.components_)
print("mean:      ", pc.mean_)
print("explained_variance_:      ", pc.explained_variance_)
print("covariance:", pc.get_covariance()) 

pca = PCA(n_components=6)
pca = pca.fit_transform(x)

print(np.shape(pca))
plt.scatter(x[:,0],x[:,6], label = "Original")
plt.scatter(pca[:,4], pca[:,5], label = "Projected")
plt.legend(loc="best", fancybox=True, shadow=True)
plt.grid(True)
plt.title('Original Correlation vs transformed correlation for the first natural frequency.')
plt.xlabel('Principal Componnent 1')
plt.ylabel('First Natural Frequency')
plt.show()