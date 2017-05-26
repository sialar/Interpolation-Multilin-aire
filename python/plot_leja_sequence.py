from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import sys

leja_sequence_file = open( "../data/leja_sequence.txt", "r")
lines = leja_sequence_file.readlines()
leja_sequence_file.close()

leja_sequence_x, leja_sequence_y, leja_sequence_z = [], [], []
values = []

if len(sys.argv)!=3:
    print("Invalid number of arguments")
    sys.exit(0)

nb_pts = int(sys.argv[2])
dim = int(sys.argv[1])

for i in range(0,nb_pts):
    values.append(float(lines[i]))

if dim == 3:
    for v1 in values:
        for v2 in values:
            for v3 in values:
                leja_sequence_x.append(v1)
                leja_sequence_y.append(v2)
                leja_sequence_z.append(v3)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(leja_sequence_x, leja_sequence_y, leja_sequence_z, c='r')
    plt.title('Nuage de points de Leja en 3D')
    plt.savefig('../images/leja_sequence_3D.png')

else:
    for v1 in values:
        for v2 in values:
            leja_sequence_x.append(v1)
            leja_sequence_y.append(v2)
    plt.scatter(leja_sequence_x, leja_sequence_y, c='r')
    plt.title('Nuage de points de Leja en 2D')
    plt.savefig('../images/leja_sequence_2D.png')
plt.show()
