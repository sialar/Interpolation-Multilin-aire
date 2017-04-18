import matplotlib.pyplot as plt

leja_sequence_file = open( "leja_sequence.txt", "r")
lines = leja_sequence_file.readlines()
leja_sequence_file.close()

leja_sequence_x = []
leja_sequence_y = []
values_x = []
values_y = []

nb_leja_points_x = int(lines[0])
for i in range(1,nb_leja_points_x+1):
    values_x.append(float(lines[i]))

nb_leja_points_y = int(lines[nb_leja_points_x+1])
for i in range(1,nb_leja_points_y+1):
    values_y.append(float(lines[i+nb_leja_points_x+1]))

for v1 in values_x:
    for v2 in values_y:
        leja_sequence_x.append(v1)
        leja_sequence_y.append(v2)

plt.scatter(leja_sequence_x,leja_sequence_y,s=1)
plt.title('Nuage de points de Leja')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('lejaSequence_2dDistribution.png')
plt.show()
