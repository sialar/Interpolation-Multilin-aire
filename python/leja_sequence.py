from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt

leja_sequence_file = open( "leja_sequence.txt", "r")
lines = leja_sequence_file.readlines()
leja_sequence_file.close()

leja_sequence_x, leja_sequence_y, leja_sequence_z = [], [], []
values_x, values_y, values_z = [], [], []


nb_leja_points_x = int(lines[0])
for i in range(1,nb_leja_points_x+1):
    values_x.append(float(lines[i]))
offset = nb_leja_points_x+1
nb_leja_points_y = int(lines[offset])
for i in range(1,nb_leja_points_y+1):
    values_y.append(float(lines[i+offset]))
offset = nb_leja_points_x + nb_leja_points_y + 2

if len(lines) > offset:
    nb_leja_points_z = int(lines[offset])
    for i in range(1,nb_leja_points_z+1):
        values_z.append(float(lines[i+offset]))
    for v1 in values_x:
        for v2 in values_y:
            for v3 in values_z:
                leja_sequence_x.append(v1)
                leja_sequence_y.append(v2)
                leja_sequence_z.append(v3)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(leja_sequence_x, leja_sequence_y, leja_sequence_z, c='r')
    plt.title('Nuage de points de Leja en 3D')

else:
    for v1 in values_x:
        for v2 in values_y:
            leja_sequence_x.append(v1)
            leja_sequence_y.append(v2)
    plt.scatter(leja_sequence_x, leja_sequence_y, c='r')
    plt.title('Nuage de points de Leja en 2D')

plt.savefig('leja_sequence.png')
plt.show()
