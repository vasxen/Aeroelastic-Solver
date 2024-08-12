from Preprocessor import * 
from Solver import *

Wing = WingPlanform(Wingspan = 12, Croot = 6, Ctip = 6, Sweepback= np.deg2rad(5))
Nspan = 15
Nchord = 10
NodeMatrix, PanelNodes = CreateWingNodes(Wing, Nspan,Nchord)
Elematrix = ConnectStructuredGrid(Nspan, Nchord)
Aeropanels = CreateCaeroPanels(Nspan , Nchord + 1, list(range(Nspan*Nchord, 2*Nspan*Nchord + Nspan +1 )) )
FixedDoFs = ClampLeftEdge(Nspan, Nchord)
Mat1 = IsotropicMaterial(1, 1600, 9E10, 0.33)
Prop1 = PSHELL(1, Mat1, 0.001)
Caero = CAERO(PanelNodes, Aeropanels)
PropertyAssignment = AssignProperty(Elematrix[:,0].tolist(), Prop1)


ShellElements = Q4Elements(Elematrix, NodeMatrix, PropertyAssignment) #type: ignore
K = ShellElements.StifnessMatrix
M = ShellElements.MassMatrix
K = (K.T + K ) / 2 
M = (M.T + M ) / 2

W, V = eig(K, M) # type: ignore
Wr = W.real
sorted_indices = np.argsort(Wr)
sorted_eigenvalues = Wr[sorted_indices]
# print(sorted_eigenvalues)
V = V[:, sorted_indices]

Displacements = np.zeros_like(NodeMatrix)
for i in range(NodeMatrix.shape[0]):
    Displacements[i,:] = NodeMatrix[i,:] +  V[[6*i, 6*i+1, 6*i +2], 13]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(Displacements[:,0],Displacements[:,1],Displacements[:,2])

for element in Elematrix:
    # Extract x, y, z coordinates for each point in the element
    e = np.append(element[1:], element[1])
    x_coords = Displacements[e, 0]
    y_coords = Displacements[e, 1]
    z_coords = Displacements[e, 2]
    
    # Plot the line connecting the four points
    ax.plot(x_coords, y_coords, z_coords, color = 'black')
    
    


plt.xlabel("X")
plt.ylabel("Y")
plt.show()



