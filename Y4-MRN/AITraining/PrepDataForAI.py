import numpy as np
import matplotlib.pyplot as plt


RawDataSet1 = np.load('Dataset1.npy', allow_pickle=True)
RawDataSet2 = np.load('Dataset2.npy', allow_pickle=True)
RawDataSet3 = np.load('Dataset3.npy', allow_pickle=True)
RawDataSet4 = np.load('Dataset4.npy', allow_pickle=True)
RawDataSet5 = np.load('Dataset5.npy', allow_pickle=True)

CompleteDataSet = np.append(RawDataSet1,RawDataSet2,axis=0)
CompleteDataSet = np.append(CompleteDataSet,RawDataSet3,axis=0)
CompleteDataSet = np.append(CompleteDataSet,RawDataSet4,axis=0)
CompleteDataSet = np.append(CompleteDataSet,RawDataSet5,axis=0)

CompleteDataSet[:,6] = np.log10(CompleteDataSet[:,6]*CompleteDataSet[:,0]*np.sqrt(CompleteDataSet[:,4]/CompleteDataSet[:,3]))
CompleteDataSet[:,7] = np.log10(CompleteDataSet[:,7]*CompleteDataSet[:,0]*np.sqrt(CompleteDataSet[:,4]/CompleteDataSet[:,3]))
CompleteDataSet[:,8] = np.log10(CompleteDataSet[:,8]*CompleteDataSet[:,0]*np.sqrt(CompleteDataSet[:,4]/CompleteDataSet[:,3]))
CompleteDataSet[:,9] = np.log10(CompleteDataSet[:,9]*CompleteDataSet[:,0]*np.sqrt(CompleteDataSet[:,4]/CompleteDataSet[:,3]))
CompleteDataSet[:,10] = np.log10(CompleteDataSet[:,10]*CompleteDataSet[:,0]*np.sqrt(CompleteDataSet[:,4]/CompleteDataSet[:,3]))
CompleteDataSet[:,11] = np.log10(CompleteDataSet[:,11]*CompleteDataSet[:,0]*np.sqrt(CompleteDataSet[:,4]/CompleteDataSet[:,3]))
CompleteDataSet[:,12] = np.log10(CompleteDataSet[:,12]*CompleteDataSet[:,0]*np.sqrt(CompleteDataSet[:,4]/CompleteDataSet[:,3]))
CompleteDataSet[:,13] = np.log10(CompleteDataSet[:,13]*CompleteDataSet[:,0]*np.sqrt(CompleteDataSet[:,4]/CompleteDataSet[:,3]))

CompleteDataSetMean = np.mean(CompleteDataSet,axis=0)
CompleteDataSetStd = np.std(CompleteDataSet,axis=0)

IdxSetToDel =[]

print(CompleteDataSet[0])
#print(CompleteDataSetStd)

for i in range(len(CompleteDataSet)):
    if (CompleteDataSet[i,6:]<(CompleteDataSetMean[6:] - 2*CompleteDataSetStd[6:])).any():
        IdxSetToDel.append(i)



CleanedDataSet = np.delete(CompleteDataSet,IdxSetToDel,axis=0)



CleanedDataSet[:,1] = CleanedDataSet[:,0]/CleanedDataSet[:,1]
CleanedDataSet[:,2] = CleanedDataSet[:,2]/CleanedDataSet[:,1]

DataSetMin = np.min(CleanedDataSet,axis=0)
DataSetMax = np.max(CleanedDataSet,axis=0)

#DataSetMin[0],DataSetMin[1],DataSetMin[2],DataSetMin[3],DataSetMin[4],DataSetMin[5] = 0.0005,5  ,5  ,50E9 ,800 ,0.15
#DataSetMax[0],DataSetMax[1],DataSetMax[2],DataSetMax[4],DataSetMax[4],DataSetMax[5] = 0.01  ,100,100,500E9,10E3,0.4


#DataSetMin[0],DataSetMin[1],DataSetMin[2],DataSetMin[3],DataSetMin[4],DataSetMin[5] = 0.0005, 0.0005/(0.01*100) ,(0.0005*5)/(0.01*100)    ,50E9 ,800 ,0.15
#DataSetMax[0],DataSetMax[1],DataSetMax[2],DataSetMax[4],DataSetMax[4],DataSetMax[5] = 0.01  ,0.01/(0.0005*5)    ,1    ,500E9,10E3,0.4

#
#Gradient = 1/(DataSetMax-DataSetMin)
#C_Intercept = -(DataSetMin*Gradient)
#CleanedDataSet = Gradient*CleanedDataSet + C_Intercept
#
#np.save('Gradients.npy',Gradient)
#np.save('C_Intercept.npy',C_Intercept)

np.save('CleanedDataSet.npy',CleanedDataSet)


TrainingSet = CleanedDataSet[0:38400]
TestSet = CleanedDataSet[38400:48000]
ValidationSet = CleanedDataSet[48000:]



print('CompleteDataSet:', len(CompleteDataSet))
print('CleanedDataSet:', len(CleanedDataSet))
print('TrainingSet:', len(TrainingSet))
print('TestSet:', len(TestSet))
print('ValidationSet:', len(ValidationSet))

Plotting = False

if Plotting:


    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
    
    # Create Plot

    scatter_plot = ax.scatter3D(CleanedDataSet[:,1],CleanedDataSet[:,2],CleanedDataSet[:,6], c= CleanedDataSet[:,5])
    

    plt.colorbar(scatter_plot)
    # Show plot

    plt.show()

""" 

np.save('TrainingSetInputs',TrainingSet[:,[1,2,5]])
np.save('TrainingSetOutputs',TrainingSet[:,[6,7,8,9,10,11,12,13]])
np.save('TrainingSetAdditional',TrainingSet[:,[0,3,4]])

np.save('TestSetInputs',TestSet[:,[1,2,5]])
np.save('TestSetOutputs',TestSet[:,[6,7,8,9,10,11,12,13]])
np.save('TestSetAdditional',TestSet[:,[0,3,4]])

np.save('ValidationSetInputs',ValidationSet[:,[1,2,5]])
np.save('ValidationSetOutputs',ValidationSet[:,[6,7,8,9,10,11,12,13]])
np.save('ValidationSetAdditional',ValidationSet[:,[0,3,4]]) """